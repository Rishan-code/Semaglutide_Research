#!/bin/bash

# Stop on erro
set -e

echo "----------------------------------------------------------------"
echo "Setting up 10-mer Semaglutide Simulation"
echo "----------------------------------------------------------------"

# 1. Solvate
# Note: Input is system_10mer.gro (the box with 10 proteins)
if [ ! -f "system_10mer_solv.gro" ]; then
    echo "[1/4] Solvating system..."
    gmx solvate -cp system_10mer.gro -cs spc216.gro -o system_10mer_solv.gro -p topol_10mer.top
else
    echo "[1/4] Solvated system already exists (skipping)"
fi

# 2. Add Ions
if [ ! -f "system_10mer_final.gro" ]; then
    echo "[2/4] Adding Ions (Neutralizing + 0.15M)..."
    
    # FIX: Sanitize topology before running genion/grompp
    python3 fix_topology.py topol_10mer.top

    # Generate TPR for genion
    gmx grompp -f ions.mdp -c system_10mer_solv.gro -p topol_10mer.top -o ions.tpr -maxwarn 10
    
    # Replace SOL with ions (POT and CLA)
    # Using 'SOL' as the group to replace (usually group 13 or similar, we use -pname/-nname to specify ion names)
    # echo 13 | ... selects the group 'SOL' automatically if possible, or we let user select interactively? 
    # Better to force it if we know the group number, but group numbers change. 
    # We will use 'SOL' in the echo input.
    echo "SOL" | gmx genion -s ions.tpr -o system_10mer_final.gro -p topol_10mer.top -pname POT -nname CLA -neutral -conc 0.15
else
    echo "[2/4] Ionized system already exists (skipping)"
fi

# Configuration for Run
INIT="system_10mer_final"
MINI_PREFIX="step4.0_minimization_10mer"
EQUI_PREFIX="step4.1_equilibration_10mer"
PROD_PREFIX="step5_production_10mer"

echo "----------------------------------------------------------------"
echo "Starting Simulation Pipeline"
echo "----------------------------------------------------------------"

# 2b. Generate Index File
if [ ! -f "index.ndx" ]; then
    echo "[2.5/4] Generating Index File..."
    # Create SOLU (Protein + KWF1/Other) and SOLV (Everything else) groups
    # We rely on "Protein" and "Other" group names being present.
    # Group 17 will be created (SOLU), then Group 18 (SOLV).
    echo -e "\"Protein\" | \"Other\"\nname 17 SOLU\n! \"SOLU\"\nname 18 SOLV\nq" | gmx make_ndx -f system_10mer_final.gro -o index.ndx
fi

# 3. Minimization
if [ ! -f "${MINI_PREFIX}.gro" ]; then
    echo "[3/4] Running Energy Minimization..."
    # We use step4.0_minimization.mdp (reusing the one from single sim)
    gmx grompp -f step4.0_minimization.mdp -o "${MINI_PREFIX}.tpr" -c "${INIT}.gro" -r "${INIT}.gro" -p topol_10mer.top -maxwarn 10
    gmx mdrun -v -deffnm "${MINI_PREFIX}" -ntmpi 1
else
    echo "[3/4] Minimization already complete (skipping)"
fi

# 4. Equilibration
if [ ! -f "${EQUI_PREFIX}.gro" ]; then
    echo "[4/4] Running Equilibration..."
    # We use step4.1_equilibration.mdp
    gmx grompp -f step4.1_equilibration.mdp -o "${EQUI_PREFIX}.tpr" -c "${MINI_PREFIX}.gro" -r "${INIT}.gro" -p topol_10mer.top -n index.ndx
    gmx mdrun -v -deffnm "${EQUI_PREFIX}" -ntmpi 1
else
    echo "[4/4] Equilibration already complete (skipping)"
fi

# 5. Production
if [ ! -f "${PROD_PREFIX}.gro" ]; then
    echo "[5/5] Running Production Simulation..."
    gmx grompp -f step5_production.mdp -o "${PROD_PREFIX}.tpr" -c "${EQUI_PREFIX}.gro" -p topol_10mer.top -n index.ndx
    gmx mdrun -v -deffnm "${PROD_PREFIX}" -ntmpi 1
else
    echo "[5/5] Production already complete (skipping)"
fi

echo "----------------------------------------------------------------"
echo "Simulation Pipeline Complete!"
echo "----------------------------------------------------------------"
