#!/bin/bash
# Run_Single_Phenol.sh - Setup and Run Single Semaglutide Chain + 5 Phenols

# Exit on error
set -e

echo "----------------------------------------------------------------"
echo "Starting Single Chain Phenol Simulation Setup..."
echo "----------------------------------------------------------------"

# Directories
DIR_MAIN="Simulation_Single_Phenol"
mkdir -p "$DIR_MAIN"

# GROMACS Command Detection
if command -v gmx &> /dev/null; then
    GMX_CMD="gmx"
elif command -v gmx_mpi &> /dev/null; then
    GMX_CMD="gmx_mpi"
else
    echo "Error: No GROMACS binary found."
    exit 1
fi

# Detect Launcher (for cluster compatibility)
if [[ "$GMX_CMD" == *"mpi"* ]]; then
    if command -v srun &> /dev/null; then
        LAUNCHER="srun --mpi=pmi2" # Common SLURM flag
    elif command -v mpirun &> /dev/null; then
        LAUNCHER="mpirun"
    else
        LAUNCHER=""
    fi
else
    LAUNCHER=""
fi

EXEC_SERIAL="$LAUNCHER -n 1"
EXEC_MD="$LAUNCHER"

# Avoid double launcher for serial commands if unnecessary
if [[ "$GMX_CMD" == "gmx" ]]; then
    EXEC_SERIAL=""
    EXEC_MD=""
fi

echo "Using GROMACS: $GMX_CMD"
echo "Launcher: $EXEC_MD"

# Files
PHENOL_PDB="phenol_charmm.pdb"
PHENOL_NMOL=5 # 5 Phenols for 1 Chain (5:1 ratio)

# 1. Prepare Base System (Single Chain)
# Extract Chain 1 from the 10-mer production run (or system_10mer.gro)
INPUT_GRO="Simulation_Water/05_production/step5_production_10mer.gro"
TPR_REF="Simulation_Water/05_production/step5_production_10mer.tpr"

if [ ! -f "$INPUT_GRO" ]; then
    echo "Error: Input 10-mer simulation not found ($INPUT_GRO)"
    exit 1
fi

echo "Extracting Single Chain from 10-mer..."
if [ ! -f "$DIR_MAIN/single_chain.gro" ]; then
    # Create index for Molecule 1 (Assuming semaglutide is ~50 atoms/res, 30 res -> ~500 atoms)
    # Better approach: Use 'gmx trjconv' with index group.
    # Group 1 is usually 'Protein'. We need a sub-group.
    # Let's try to extract residue 1-31 (Chain A).
    # Using 'make_ndx' to create a group for just the first molecule.
    
    # "keep 1" -> Keep group 1 (Protein)
    # "r 1-31" -> Select residues 1 to 31 (Assuming chain 1 is res 1-31)
    # Wait, if 10-mer is one moleculetype with 10 mols, they might have unique atom/res IDs or be repetitive.
    # Safe bet: Select by atom number if we know the size, or just "chain 1".
    # Assuming standard GROMACS numbering: Chains are sequential.
    
    echo -e "ri 1-31 \n name 20 Chain_A \n q" | $EXEC_SERIAL $GMX_CMD make_ndx -f "$TPR_REF" -o "$DIR_MAIN/index_chain1.ndx"
    
    echo "Chain_A" | $EXEC_SERIAL $GMX_CMD trjconv -f "$INPUT_GRO" -s "$TPR_REF" -n "$DIR_MAIN/index_chain1.ndx" -o "$DIR_MAIN/single_chain_raw.gro" -pbc mol -center
    
    # 2. Define Box (Small cubic box, e.g., 6nm)
    $EXEC_SERIAL $GMX_CMD editconf -f "$DIR_MAIN/single_chain_raw.gro" -o "$DIR_MAIN/single_chain_box.gro" -c -d 1.2 -bt cubic
else
    echo "Single chain extraction already done."
fi

# 3. Insert Phenol
if [ ! -f "$DIR_MAIN/box_phenol.gro" ]; then
    echo "Inserting $PHENOL_NMOL Phenol molecules..."
    $EXEC_SERIAL $GMX_CMD insert-molecules -f "$DIR_MAIN/single_chain_box.gro" -ci "$PHENOL_PDB" -nmol "$PHENOL_NMOL" -o "$DIR_MAIN/box_phenol.gro"
else
    echo "Phenol insertion already done."
fi

# 4. Create Topology
if [ ! -f "$DIR_MAIN/topol_single.top" ]; then
    cat <<EOF > "$DIR_MAIN/topol_single.top"
;;
;; Single Chain Phenol Topology
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
Single Semaglutide Chain with Phenol

[ molecules ]
; Compound        #mols
PROH              1
PHEN              $PHENOL_NMOL
EOF
fi

# 5. Solvate
if [ ! -f "$DIR_MAIN/box_solvated.gro" ]; then
    echo "Solvating..."
    $EXEC_SERIAL $GMX_CMD solvate -cp "$DIR_MAIN/box_phenol.gro" -cs spc216.gro -o "$DIR_MAIN/box_solvated.gro" -p "$DIR_MAIN/topol_single.top"
fi

# 6. Add Ions
if [ ! -f "$DIR_MAIN/box_ionized.gro" ]; then
    echo "Adding Ions..."
    $EXEC_SERIAL $GMX_CMD grompp -f ions.mdp -c "$DIR_MAIN/box_solvated.gro" -p "$DIR_MAIN/topol_single.top" -o "$DIR_MAIN/ions.tpr" -maxwarn 1
    echo "SOL" | $EXEC_SERIAL $GMX_CMD genion -s "$DIR_MAIN/ions.tpr" -o "$DIR_MAIN/box_ionized.gro" -p "$DIR_MAIN/topol_single.top" -pname POT -nname CLA -neutral -conc 0.15
fi

# 7. Simulation Steps
cp step4.0_minimization.mdp "$DIR_MAIN/"
cp step4.1_equilibration.mdp "$DIR_MAIN/"
cp step5_production.mdp "$DIR_MAIN/"

cd "$DIR_MAIN"

# Min
if [ ! -f "min.gro" ]; then
    echo "Running Minimization..."
    $EXEC_SERIAL $GMX_CMD grompp -f step4.0_minimization.mdp -c box_ionized.gro -p topol_single.top -o min.tpr
    $EXEC_MD $GMX_CMD mdrun -v -deffnm min
fi

# Equil (NVT/NPT combined or just NVT for now)
if [ ! -f "equil.gro" ]; then
    echo "Running Equilibration..."
    # Make index for T-coupling (Protein+Phenol vs Water+Ions)
    # Simplified: Protein (1) and Non-Protein. Or just system.
    # Let's create SOLU (Protein) and SOLV (Rest)
    echo -e "1 | 13 \n name 19 SOLU \n ! 19 \n name 20 SOLV \n q" | $EXEC_SERIAL $GMX_CMD make_ndx -f min.gro -o index.ndx
    
    $EXEC_SERIAL $GMX_CMD grompp -f step4.1_equilibration.mdp -c min.gro -r min.gro -p topol_single.top -n index.ndx -o equil.tpr
    $EXEC_MD $GMX_CMD mdrun -v -deffnm equil
fi

# Prod
if [ ! -f "prod.gro" ]; then
    echo "Running Production..."
    $EXEC_SERIAL $GMX_CMD grompp -f step5_production.mdp -c equil.gro -p topol_single.top -n index.ndx -o prod.tpr
    $EXEC_MD $GMX_CMD mdrun -v -deffnm prod
fi

echo "Simulation Complete!"
