#!/bin/bash
# Run_Phenol.sh - Setup and Run 10-mer + Phenol Simulation

# Exit on error
set -e

echo "----------------------------------------------------------------"
echo "Starting Phenol Simulation Setup..."
echo "----------------------------------------------------------------"

# Directories
DIR_MAIN="Simulation_Phenol"
mkdir -p "$DIR_MAIN"

# DEBUG: Print environment info
echo "=== Debugging Environment ==="
hostname
echo "Modules Loaded:"
module list 2>&1 || echo "Module command failed"
echo "============================="

# Files (Assumed to be in root GROMACS_Setup)
PHENOL_PDB="phenol_charmm.pdb"
PHENOL_ITP="phenol_charmm.itp"
PHENOL_NMOL=50

# Check required files
if [ ! -f "$PHENOL_PDB" ] || [ ! -f "$PHENOL_ITP" ]; then
    echo "Error: Missing phenol files ($PHENOL_PDB or $PHENOL_ITP)"
    exit 1
fi

# GROMACS Command Detection & Launcher Setup
# We prefer 'mpirun' for MPI binaries to avoid PMPI_Init failures,
# unless 'gmx' (Thread-MPI) is used.

# GROMACS Command & Launcher Detection
# Priority:
# 1. gmx (Thread-MPI) -> No launcher needed (Simpler, less prone to MPI errors)
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
    
    # Set Exec Variables
    # Note: For srun, we might need '--exclusive' or similar if "credential expired" persists,
    # but 'mpiexec' is usually safer if available.
    echo "Using MPI Launcher: $LAUNCHER"
    
    if [[ "$LAUNCHER" == "srun" ]]; then
        # Try --mpi=pmi2 to avoid "Job credential expired" or PMI errors
        # Also --exclusive sometimes helps prevent step conflicts
        # Using --mpi=pmi2 as primary fix attempt.
        FLAGS="--mpi=pmi2"
        EXEC_SERIAL="srun $FLAGS -n 1"
        EXEC_MD="srun $FLAGS"
    else
        EXEC_SERIAL="$LAUNCHER -n 1"
        EXEC_MD="$LAUNCHER"
    fi
fi

echo "Using GROMACS execution command: $GMX_CMD"
echo "Serial Launcher: '$EXEC_SERIAL'"
echo "MD Launcher:     '$EXEC_MD'"

# Fix Atom Name Mismatch (OH1 vs OH) in PDB if present
# ITP uses "OH", PDB might have "OH1" from earlier rename.
sed -i 's/OH1 /OH  /g' "$PHENOL_PDB"

# 1. Prepare Base System (Protein Only)
# We look for the Water simulation output in Simulation_Water
INPUT_GRO="Simulation_Water/05_production/step5_production_10mer.gro"
TPR_REF="Simulation_Water/05_production/step5_production_10mer.tpr"

if [ ! -f "$INPUT_GRO" ]; then
    # Fallback to equilibration if production files missing
    INPUT_GRO="Simulation_Water/04_equilibration/step4.1_equilibration_10mer.gro"
    TPR_REF="Simulation_Water/04_equilibration/step4.1_equilibration_10mer.tpr"  # Rough approx for PBC
fi

echo "Using input structure: $INPUT_GRO"

if [ ! -f "$DIR_MAIN/protein_only.gro" ]; then
    # Extract Protein (Group 1 usually)
    # CRITICAL FIX: PROH topology includes KWF1 (SNAC?), but 'Protein' group excludes it.
    # We must merge Protein (1) and KWF1 (13) to match the PROH molecule definition.
    echo "Creating merged index (Protein + KWF1)..."
    # "1 | 13" creates new group. "name 19 PROH_Complex" renames it. "q" to quit.
    # We assume groups 0-18 exist (based on previous logs). New group is 19.
    echo -e "1 | 13 \n name 19 PROH_Complex \n q" | $EXEC_SERIAL $GMX_CMD make_ndx -f "$TPR_REF" -o "$DIR_MAIN/index_merged.ndx"

    echo "Extracting Protein + KWF1..."
    # Center on Protein (1), Output PROH_Complex (19)
    echo "1 19" | $EXEC_SERIAL $GMX_CMD trjconv -f "$INPUT_GRO" -s "$TPR_REF" -n "$DIR_MAIN/index_merged.ndx" -o "$DIR_MAIN/protein_only.gro" -pbc mol -center
else
    echo "Step 1: Protein Extraction already done. Skipping."
fi

# 2. Insert Phenol
if [ ! -f "$DIR_MAIN/box_phenol.gro" ]; then
    echo "Inserting $PHENOL_NMOL Phenol molecules..."
    $EXEC_SERIAL $GMX_CMD insert-molecules -f "$DIR_MAIN/protein_only.gro" -ci "$PHENOL_PDB" -nmol "$PHENOL_NMOL" -o "$DIR_MAIN/box_phenol.gro"
else
    echo "Step 2: Phenol Insertion already done. Skipping."
fi

# 3. Create Topology (Fresh Write)
if [ ! -f "$DIR_MAIN/topol_phenol.top" ]; then
    # Simulation_Phenol/topol_phenol.top
    cat <<EOF > "$DIR_MAIN/topol_phenol.top"
;;
;; Generated by Run_Phenol.sh for Phenol Simulation
;;

; Include forcefield parameters
#include "../toppar/forcefield.itp"

; Ligand Topology
#include "../phenol_charmm.itp"

; Protein and Water Topologies
#include "../toppar/PROH.itp"
#include "../toppar/SOL.itp"
#include "../toppar/POT.itp"
#include "../toppar/CLA.itp"

[ system ]
; Name
Semaglutide 10-mer with Phenol

[ molecules ]
; Compound        #mols
PROH              10
PHEN              $PHENOL_NMOL
EOF
else
    echo "Step 3: Topology already exists. Skipping."
fi

# 4. Solvate
if [ ! -f "$DIR_MAIN/box_solvated.gro" ]; then
    echo "Solvating..."
    $EXEC_SERIAL $GMX_CMD solvate -cp "$DIR_MAIN/box_phenol.gro" -cs spc216.gro -o "$DIR_MAIN/box_solvated.gro" -p "$DIR_MAIN/topol_phenol.top"
else
    echo "Step 4: Solvation already done. Skipping."
fi

# 5. Add Ions
if [ ! -f "$DIR_MAIN/box_ionized.gro" ]; then
    echo "Adding Ions..."
    # Create TPR for genion
    $EXEC_SERIAL $GMX_CMD grompp -f ions.mdp -c "$DIR_MAIN/box_solvated.gro" -p "$DIR_MAIN/topol_phenol.top" -o "$DIR_MAIN/ions.tpr" -maxwarn 1
    
    # Run genion
    echo "SOL" | $EXEC_SERIAL $GMX_CMD genion -s "$DIR_MAIN/ions.tpr" -o "$DIR_MAIN/box_ionized.gro" -p "$DIR_MAIN/topol_phenol.top" -pname POT -nname CLA -neutral -conc 0.15
else
    echo "Step 5: Ionization already done. Skipping."
fi

# 6. Setup Simulation Steps (Min -> Equil -> Prod)
# Copy MDP files (Always update in case of changes)
cp step4.0_minimization.mdp "$DIR_MAIN/"
cp step4.1_equilibration.mdp "$DIR_MAIN/"
cp step5_production.mdp "$DIR_MAIN/"

cd "$DIR_MAIN"

# Minimization
if [ ! -f "step4.0_minimization_phenol.gro" ]; then
    echo "Starting Minimization..."
    $EXEC_SERIAL $GMX_CMD grompp -f step4.0_minimization.mdp -c box_ionized.gro -r box_ionized.gro -p topol_phenol.top -o step4.0_minimization_phenol.tpr
    $EXEC_MD $GMX_CMD mdrun -v -deffnm step4.0_minimization_phenol
else
    echo "Step 6a: Minimization already done. Skipping."
fi

# Equilibration
if [ ! -f "step4.1_equilibration_phenol.gro" ]; then
    echo "Creating Index Groups for T-Coupling (SOLU/SOLV)..."
    # Define SOLU (Protein + KWF1) and SOLV (Everything else)
    # Using Group IDs from log: 1=Protein, 13=KWF1 (adjust if needed)
    echo -e "1 | 13 \n name 20 SOLU \n ! 20 \n name 21 SOLV \n q" | $EXEC_SERIAL $GMX_CMD make_ndx -f step4.0_minimization_phenol.gro -o index_phenol.ndx

    echo "Starting Equilibration..."
    $EXEC_SERIAL $GMX_CMD grompp -f step4.1_equilibration.mdp -c step4.0_minimization_phenol.gro -r step4.0_minimization_phenol.gro -p topol_phenol.top -n index_phenol.ndx -o step4.1_equilibration_phenol.tpr
    $EXEC_MD $GMX_CMD mdrun -v -deffnm step4.1_equilibration_phenol
else
    echo "Step 6b: Equilibration already done. Skipping."
fi

# Production (Resume capable)
echo "Starting/Continuing Production (Target: 10 ns)..."

# Always grompp to allow updating parameters (like nsteps)
$EXEC_SERIAL $GMX_CMD grompp -f step5_production.mdp -c step4.1_equilibration_phenol.gro -p topol_phenol.top -n index_phenol.ndx -o step5_production_phenol.tpr

# Check if we should resume from a checkpoint
if [ -f "step5_production_phenol.cpt" ]; then
    echo "Resuming from checkpoint..."
    $EXEC_MD $GMX_CMD mdrun -v -deffnm step5_production_phenol -cpi step5_production_phenol.cpt -append
else
    echo "Starting fresh production run..."
    $GMX_CMD mdrun -v -deffnm step5_production_phenol
fi

# --- 7. Analysis (Rg) ---
# Analysis logic for Phenol Simulation
if [ -f "step5_production_phenol.tpr" ]; then
    echo "Calculating Radius of Gyration (Protein)..."
    # Select Group 1 (Protein) for Rg calculation
    # Note: Ensure 'Protein' group exists or use index_phenol.ndx if needed. 
    # Usually Group 1 is Protein in standard index files.
    echo "1" | $EXEC_SERIAL $GMX_CMD gyrate -s step5_production_phenol.tpr -f step5_production_phenol.gro -o gyrate_phenol.xvg
    echo "Rg analysis saved to gyrate_phenol.xvg"
fi

echo "----------------------------------------------------------------"
echo "Run Complete check"
echo "----------------------------------------------------------------"
