#!/bin/bash

# Stop on error
set -e

# Configuration
INIT="step3_input"
MINI_PREFIX="step4.0_minimization"
EQUI_PREFIX="step4.1_equilibration"
PROD_PREFIX="step5_production"

echo "----------------------------------------------------------------"
echo "Starting Semaglutide Simulation (WSL/Linux Mode)"
echo "----------------------------------------------------------------"

# 1. Minimization
echo "[1/3] Running Energy Minimization..."
gmx grompp -f "${MINI_PREFIX}.mdp" -o "${MINI_PREFIX}.tpr" -c "${INIT}.gro" -r "${INIT}.gro" -p topol.top -n index.ndx -maxwarn 10
gmx mdrun -v -deffnm "${MINI_PREFIX}" -ntmpi 1

# 2. Equilibration
echo "[2/3] Running Equilibration..."
gmx grompp -f "${EQUI_PREFIX}.mdp" -o "${EQUI_PREFIX}.tpr" -c "${MINI_PREFIX}.gro" -r "${INIT}.gro" -p topol.top -n index.ndx
gmx mdrun -v -deffnm "${EQUI_PREFIX}" -ntmpi 1

# 3. Production
echo "[3/3] Running Production (Part 1)..."
# Using _1 suffix for the first part of production
PROD_STEP="${PROD_PREFIX}_1"

gmx grompp -f "${PROD_PREFIX}.mdp" -o "${PROD_STEP}.tpr" -c "${EQUI_PREFIX}.gro" -p topol.top -n index.ndx
gmx mdrun -v -deffnm "${PROD_STEP}" -ntmpi 1

echo "----------------------------------------------------------------"
echo "Simulation Complete!"
echo "----------------------------------------------------------------"
