#!/bin/bash
# Run_Water.sh - Setup and Run 10-mer Water Simulation (100ns)
# Adapted from Run_Phenol.sh for consistency and GPU support.

# Exit on error
set -e

echo "----------------------------------------------------------------"
echo "Starting Water Simulation Setup (Semaglutide 10-mer)..."
echo "----------------------------------------------------------------"

# Directories
DIR_MAIN="Simulation_Water"
mkdir -p "$DIR_MAIN"

# DEBUG: Print environment info
echo "=== Debugging Environment ==="
hostname
echo "Modules Loaded:"
module list 2>&1 || echo "Module command failed"
echo "============================="

# GROMACS Command Detection & Launcher Setup
# Priority:
# 1. gmx (Thread-MPI) -> No launcher needed
# 2. gmx_mpi / gmx-mpi -> Needs mpirun/mpiexec/srun

GMX_CMD=""
LAUNCHER=""
EXEC_SERIAL=""
EXEC_MD=""

if command -v gmx &> /dev/null; then
    GMX_CMD="gmx"
    echo "Found 'gmx' (Thread-MPI). No external launcher required."
elif command -v gmx_mpi &> /dev/null; then
    GMX_CMD="gmx_mpi"
    echo "Found 'gmx_mpi' (MPI-enabled)."
elif command -v gmx-mpi &> /dev/null; then
    GMX_CMD="gmx-mpi"
    echo "Found 'gmx-mpi' (MPI-enabled)."
else
    echo "Error: No GROMACS binary (gmx, gmx_mpi, gmx-mpi) found in PATH."
    echo "Please ensure 'module load gromacs...' is in your submission script."
    exit 1
fi

# Detect Launcher if using MPI version
if [[ "$GMX_CMD" == *"mpi"* ]]; then
    if command -v mpirun &> /dev/null; then
        LAUNCHER="mpirun"
    elif command -v mpiexec &> /dev/null; then
        LAUNCHER="mpiexec"
    elif command -v srun &> /dev/null; then
        LAUNCHER="srun"
    else
        echo "Error: MPI GROMACS found ($GMX_CMD) but no launcher (mpirun, mpiexec, srun) found."
        exit 1
    fi
    
    echo "Using MPI Launcher: $LAUNCHER"
    
    if [[ "$LAUNCHER" == "srun" ]]; then
        FLAGS="--mpi=pmi2"
        EXEC_SERIAL="srun $FLAGS -n 1"
        EXEC_MD="srun $FLAGS"
    else
        EXEC_SERIAL="$LAUNCHER -n 1"
        EXEC_MD="$LAUNCHER"
    fi
else
    # Thread-MPI (gmx)
    EXEC_SERIAL=""
    EXEC_MD=""
fi

echo "Using GROMACS execution command: $GMX_CMD"
echo "Serial Launcher: '$EXEC_SERIAL'"
echo "MD Launcher:     '$EXEC_MD'"

# --- Files & Paths ---
# We expect the initial structure and topology to be in Simulation_Water
# existing files: system_10mer.gro, topol_10mer.top
cd "$DIR_MAIN"

INPUT_GRO="system_10mer.gro"
TOPOL="topol_10mer.top"

if [ ! -f "$INPUT_GRO" ]; then
    echo "Error: $INPUT_GRO not found in $DIR_MAIN"
    exit 1
fi

# --- 1. Solvation ---
SOLV_GRO="01_solvation/system_10mer_solv.gro"
mkdir -p 01_solvation

if [ ! -f "$SOLV_GRO" ]; then
    echo "[1/5] Solvating system..."
    $EXEC_SERIAL $GMX_CMD solvate -cp "$INPUT_GRO" -cs spc216.gro -o "$SOLV_GRO" -p "$TOPOL"
else
    echo "[1/5] Solvated system already exists (skipping)"
fi

# --- 2. Ions ---
IONS_DIR="02_ions"
ION_GRO="${IONS_DIR}/system_10mer_final.gro"
IONS_TPR="${IONS_DIR}/ions.tpr"
mkdir -p "$IONS_DIR"

if [ ! -f "$ION_GRO" ]; then
    echo "[2/5] Adding Ions..."
    
    # Fix topology if needed (calling the python script from root if it exists)
    if [ -f "../fix_topology.py" ]; then
        python3 ../fix_topology.py "$TOPOL"
    fi

    # We need ions.mdp. Check if it exists here, else cp from root
    if [ ! -f "ions.mdp" ] && [ -f "../ions.mdp" ]; then
        cp "../ions.mdp" .
    fi

    # Generate TPR
    $EXEC_SERIAL $GMX_CMD grompp -f ions.mdp -c "$SOLV_GRO" -p "$TOPOL" -o "$IONS_TPR" -maxwarn 10
    
    # Genion
    echo "SOL" | $EXEC_SERIAL $GMX_CMD genion -s "$IONS_TPR" -o "$ION_GRO" -p "$TOPOL" -pname POT -nname CLA -neutral -conc 0.15
else
    echo "[2/5] Ionized system already exists (skipping)"
fi

# --- 3. Minimization ---
MIN_DIR="03_minimization"
MIN_PREFIX="step4.0_minimization_10mer"
MIN_GRO="${MIN_DIR}/${MIN_PREFIX}.gro"
MIN_TPR="${MIN_DIR}/${MIN_PREFIX}.tpr"
mkdir -p "$MIN_DIR"

# Ensure mdp exists
if [ ! -f "step4.0_minimization.mdp" ] && [ -f "../step4.0_minimization.mdp" ]; then
    cp "../step4.0_minimization.mdp" .
fi

if [ ! -f "$MIN_GRO" ]; then
    echo "[3/5] Running Minimization..."
    $EXEC_SERIAL $GMX_CMD grompp -f step4.0_minimization.mdp -o "$MIN_TPR" -c "$ION_GRO" -r "$ION_GRO" -p "$TOPOL" -maxwarn 10
    
    # Run
    $EXEC_MD $GMX_CMD mdrun -v -deffnm "${MIN_DIR}/${MIN_PREFIX}"
else
    echo "[3/5] Minimization already complete (skipping)"
fi

# --- 4. Equilibration ---
EQUI_DIR="04_equilibration"
EQUI_PREFIX="step4.1_equilibration_10mer"
EQUI_GRO="${EQUI_DIR}/${EQUI_PREFIX}.gro"
EQUI_TPR="${EQUI_DIR}/${EQUI_PREFIX}.tpr"
mkdir -p "$EQUI_DIR"

# Ensure mdp exists
if [ ! -f "step4.1_equilibration.mdp" ] && [ -f "../step4.1_equilibration.mdp" ]; then
    cp "../step4.1_equilibration.mdp" .
fi

if [ ! -f "$EQUI_GRO" ]; then
    echo "[4/5] Running Equilibration..."
    
    # Index file needed?
    if [ ! -f "index.ndx" ]; then
        echo "Generating Index file..."
        echo -e "\"Protein\" | \"Other\"\nname 17 SOLU\n! \"SOLU\"\nname 18 SOLV\nq" | $EXEC_SERIAL $GMX_CMD make_ndx -f "$ION_GRO" -o index.ndx
    fi

    $EXEC_SERIAL $GMX_CMD grompp -f step4.1_equilibration.mdp -o "$EQUI_TPR" -c "$MIN_GRO" -r "$ION_GRO" -p "$TOPOL" -n index.ndx
    $EXEC_MD $GMX_CMD mdrun -v -deffnm "${EQUI_DIR}/${EQUI_PREFIX}"
else
    echo "[4/5] Equilibration already complete (skipping)"
fi

# --- 5. Production ---
PROD_DIR="05_production"
PROD_PREFIX="step5_production_10mer"
PROD_GRO="${PROD_DIR}/${PROD_PREFIX}.gro"
PROD_TPR="${PROD_DIR}/${PROD_PREFIX}.tpr"
PROD_CPT="${PROD_DIR}/${PROD_PREFIX}.cpt"
mkdir -p "$PROD_DIR"

# Ensure mdp exists
if [ ! -f "step5_production.mdp" ] && [ -f "../step5_production.mdp" ]; then
    cp "../step5_production.mdp" .
fi
# Note: We likely already updated step5_production.mdp in Simulation_Water/ during the previous step

echo "----------------------------------------------------------------"
echo "Starting Production Pipeline (100ns)"
echo "----------------------------------------------------------------"

# Always run grompp to ensure we have the latest MDP settings (100ns)
# Only if TPR doesn't exist or we strictly want to update? 
# Robustness: existing TPR might be old 1ns version. 
# Safe approach: Backup old TPR if it exists and nsteps differs? 
# For now, we will regenerate TPR if the CPT doesn't exist, OR if we force it.
# Actually, the user might be extending. 
# If CPT exists, we prioritize appending.

# Always generate Production TPR to ensure latest MDP settings (100ns) are used
echo "Generating Production TPR..."
$EXEC_SERIAL $GMX_CMD grompp -f step5_production.mdp -o "$PROD_TPR" -c "$EQUI_GRO" -p "$TOPOL" -n index.ndx

if [ -f "$PROD_CPT" ]; then
    echo "Resuming from checkpoint..."
    $EXEC_MD $GMX_CMD mdrun -v -deffnm "${PROD_DIR}/${PROD_PREFIX}" -cpi "${PROD_DIR}/${PROD_PREFIX}.cpt" -append
else
    if [ ! -f "$PROD_GRO" ]; then
        echo "Starting fresh production run..."
        $EXEC_MD $GMX_CMD mdrun -v -deffnm "${PROD_DIR}/${PROD_PREFIX}"
    else
        echo "Production run completed (GRO file exists)."
    fi
fi

# --- 6. Analysis (Rg) ---
RG_OUT="analysis/gyrate.xvg"
mkdir -p "analysis"

if [ -f "${PROD_DIR}/${PROD_PREFIX}.tpr" ]; then
    echo "Calculating Radius of Gyration..."
    # Select Group 1 (Protein) for Rg calculation
    echo "1" | $EXEC_SERIAL $GMX_CMD gyrate -s "${PROD_DIR}/${PROD_PREFIX}.tpr" -f "${PROD_DIR}/${PROD_PREFIX}.gro" -o "$RG_OUT"
    echo "Rg analysis saved to $RG_OUT"
fi

echo "----------------------------------------------------------------"
echo "Water Simulation Pipeline Complete!"
echo "----------------------------------------------------------------"
