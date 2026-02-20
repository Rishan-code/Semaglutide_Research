#!/bin/bash
# =============================================================================
# Run_Single_Water.sh
# Single Semaglutide Chain in WATER ONLY (no phenol) — 150 ns Production
# Control experiment for comparison with phenol simulation
# =============================================================================
set -e

echo "================================================================"
echo "  Single Semaglutide Chain — Water Only (Control)"
echo "================================================================"

# --- Configuration ---
DIR_MAIN="Simulation_Single_Water"
PROD_NSTEPS=75000000   # 150 ns at dt=0.002

mkdir -p "$DIR_MAIN"

# --- GROMACS Detection ---
if command -v gmx &> /dev/null; then
    GMX="gmx"
elif command -v gmx_mpi &> /dev/null; then
    GMX="gmx_mpi"
else
    echo "Error: GROMACS not found." && exit 1
fi

# Launcher for MPI builds
LAUNCHER=""
if [[ "$GMX" == *"mpi"* ]]; then
    if command -v srun &> /dev/null; then
        LAUNCHER="srun --mpi=pmi2"
    elif command -v mpirun &> /dev/null; then
        LAUNCHER="mpirun"
    fi
fi

SERIAL="$LAUNCHER"
MDRUN="$LAUNCHER"
if [[ "$GMX" == "gmx" ]]; then
    SERIAL=""
    MDRUN=""
fi

echo "GROMACS: $GMX | Launcher: ${MDRUN:-none}"

# --- Input Files ---
INPUT_GRO="Simulation_Water/05_production/step5_production_10mer.gro"
TPR_REF="Simulation_Water/05_production/step5_production_10mer.tpr"

if [ ! -f "$INPUT_GRO" ] || [ ! -f "$TPR_REF" ]; then
    echo "Error: 10-mer simulation files not found."
    echo "  Need: $INPUT_GRO"
    echo "  Need: $TPR_REF"
    exit 1
fi

# =============================================================================
# STEP 1: Extract Single Chain from 10-mer
# =============================================================================
if [ ! -f "$DIR_MAIN/single_chain.gro" ]; then
    echo ""
    echo "--- Step 1: Extracting Single Chain ---"

    # Select residue indices 1-30 (one PROH molecule = 30 residues)
    echo -e "ri 1-30\nq" | $SERIAL $GMX make_ndx -f "$TPR_REF" -o "$DIR_MAIN/index_extract.ndx"

    # Extract the last created group
    LAST_GROUP=$(grep "^\[" "$DIR_MAIN/index_extract.ndx" | wc -l)
    LAST_GROUP=$((LAST_GROUP - 1))

    echo "$LAST_GROUP" | $SERIAL $GMX trjconv \
        -f "$INPUT_GRO" \
        -s "$TPR_REF" \
        -n "$DIR_MAIN/index_extract.ndx" \
        -o "$DIR_MAIN/single_chain_raw.gro" \
        -pbc mol

    # Define box (cubic, 1.2 nm padding)
    $SERIAL $GMX editconf \
        -f "$DIR_MAIN/single_chain_raw.gro" \
        -o "$DIR_MAIN/single_chain.gro" \
        -c -d 1.2 -bt cubic

    echo "  Chain extracted and boxed."
else
    echo "Step 1: Single chain already extracted. Skipping."
fi

# =============================================================================
# STEP 2: Create Topology (NO phenol — water only)
# =============================================================================
if [ ! -f "$DIR_MAIN/topol_single_water.top" ]; then
    echo ""
    echo "--- Step 2: Creating Topology ---"

    cat <<EOF > "$DIR_MAIN/topol_single_water.top"
;;
;; Single Semaglutide in Water — Control Topology (CHARMM36)
;;

; Force field
#include "../toppar/forcefield.itp"

; Protein (Semaglutide)
#include "../toppar/PROH.itp"

; Water & Ions
#include "../toppar/SOL.itp"
#include "../toppar/POT.itp"
#include "../toppar/CLA.itp"

[ system ]
Single Semaglutide Chain in Water

[ molecules ]
; Compound        #mols
PROH              1
EOF

    echo "  Topology created (water only, no phenol)."
else
    echo "Step 2: Topology already exists. Skipping."
fi

# =============================================================================
# STEP 3: Solvate
# =============================================================================
if [ ! -f "$DIR_MAIN/box_solvated.gro" ]; then
    echo ""
    echo "--- Step 3: Solvating ---"

    $SERIAL $GMX solvate \
        -cp "$DIR_MAIN/single_chain.gro" \
        -cs spc216.gro \
        -o "$DIR_MAIN/box_solvated.gro" \
        -p "$DIR_MAIN/topol_single_water.top"

    echo "  System solvated."
else
    echo "Step 3: Already solvated. Skipping."
fi

# =============================================================================
# STEP 4: Add Ions (neutralize + 0.15M KCl)
# =============================================================================
if [ ! -f "$DIR_MAIN/box_ionized.gro" ]; then
    echo ""
    echo "--- Step 4: Adding Ions ---"

    $SERIAL $GMX grompp \
        -f ions.mdp \
        -c "$DIR_MAIN/box_solvated.gro" \
        -p "$DIR_MAIN/topol_single_water.top" \
        -o "$DIR_MAIN/ions.tpr" \
        -maxwarn 1

    echo "SOL" | $SERIAL $GMX genion \
        -s "$DIR_MAIN/ions.tpr" \
        -o "$DIR_MAIN/box_ionized.gro" \
        -p "$DIR_MAIN/topol_single_water.top" \
        -pname POT -nname CLA -neutral -conc 0.15

    echo "  Ions added."
else
    echo "Step 4: Already ionized. Skipping."
fi

# =============================================================================
# STEP 5: Copy MDP files & create production MDP for 150 ns
# =============================================================================
cd "$DIR_MAIN"

cp ../step4.0_minimization.mdp .
cp ../step4.1_equilibration.mdp .

# Create production MDP with 150 ns
cat <<EOF > step5_production.mdp
integrator              = md
dt                      = 0.002
nsteps                  = $PROD_NSTEPS  ; 150 ns
nstxout-compressed      = 50000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
rlist                   = 1.2
rcoulomb                = 1.2
coulombtype             = PME
;
tcoupl                  = v-rescale
tc_grps                 = SOLU SOLV
tau_t                   = 1.0 1.0
ref_t                   = 303.15 303.15
;
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 5.0
compressibility         = 4.5e-5
ref_p                   = 1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SOLU SOLV
;
EOF

echo "  MDP files prepared (production = 150 ns)."

# =============================================================================
# STEP 6: Energy Minimization
# =============================================================================
if [ ! -f "min.gro" ]; then
    echo ""
    echo "--- Step 6: Energy Minimization ---"

    $SERIAL $GMX grompp \
        -f step4.0_minimization.mdp \
        -c box_ionized.gro \
        -r box_ionized.gro \
        -p topol_single_water.top \
        -o min.tpr

    $MDRUN $GMX mdrun -v -deffnm min

    echo "  Minimization complete."
else
    echo "Step 6: Minimization already done. Skipping."
fi

# =============================================================================
# STEP 7: Create Index Groups (SOLU / SOLV)
# =============================================================================
if [ ! -f "index.ndx" ]; then
    echo ""
    echo "--- Step 7: Creating Index Groups ---"

    # For water-only system: SOLU = Protein, SOLV = everything else
    # Robust: don't hardcode group numbers — they vary across systems
    # Step 1: Create complement of Protein (group 1)
    echo -e "! 1\nq" | $SERIAL $GMX make_ndx -f min.gro -o index_raw.ndx

    # Step 2: Find the new group number (last group)
    NGRP=$(grep -c '^\[' index_raw.ndx)
    LAST=$((NGRP - 1))

    # Step 3: Rename Protein→SOLU and complement→SOLV
    echo -e "name 1 SOLU\nname $LAST SOLV\nq" | \
        $SERIAL $GMX make_ndx -f min.gro -n index_raw.ndx -o index.ndx

    rm -f index_raw.ndx

    echo "  Index groups created."
else
    echo "Step 7: Index already exists. Skipping."
fi

# =============================================================================
# STEP 8: Equilibration (125 ps NVT with position restraints)
# =============================================================================
if [ ! -f "equil.gro" ]; then
    echo ""
    echo "--- Step 8: Equilibration ---"

    $SERIAL $GMX grompp \
        -f step4.1_equilibration.mdp \
        -c min.gro \
        -r min.gro \
        -p topol_single_water.top \
        -n index.ndx \
        -o equil.tpr

    $MDRUN $GMX mdrun -v -deffnm equil

    echo "  Equilibration complete."
else
    echo "Step 8: Equilibration already done. Skipping."
fi

# =============================================================================
# STEP 9: Production MD (150 ns)
# =============================================================================
if [ ! -f "prod.gro" ]; then
    echo ""
    echo "--- Step 9: Production MD (150 ns) ---"

    $SERIAL $GMX grompp \
        -f step5_production.mdp \
        -c equil.gro \
        -p topol_single_water.top \
        -n index.ndx \
        -o prod.tpr

    $MDRUN $GMX mdrun -v -deffnm prod

    echo "  Production complete."
else
    echo "Step 9: Production already done. Skipping."
fi

# =============================================================================
# STEP 10: Fix PBC — Clean Trajectory for Analysis & Visualization
# =============================================================================
echo ""
echo "--- Step 10: Fixing PBC (Periodic Boundary Conditions) ---"

if [ -f "prod.xtc" ] && [ -f "prod.tpr" ] && [ ! -f "prod_noPBC.xtc" ]; then
    # Step A: Remove jumps across box boundaries
    echo "SOLU" | $SERIAL $GMX trjconv \
        -f prod.xtc \
        -s prod.tpr \
        -n index.ndx \
        -o prod_nojump.xtc \
        -pbc nojump

    # Step B: Center protein in box and make molecules whole
    echo -e "SOLU\nSystem" | $SERIAL $GMX trjconv \
        -f prod_nojump.xtc \
        -s prod.tpr \
        -n index.ndx \
        -o prod_noPBC.xtc \
        -pbc mol -center

    # Clean up intermediate file
    rm -f prod_nojump.xtc

    echo "  PBC-corrected trajectory: prod_noPBC.xtc"
else
    if [ -f "prod_noPBC.xtc" ]; then
        echo "  PBC fix already done. Skipping."
    else
        echo "  WARNING: Production files not found — skipping PBC fix."
    fi
fi

# =============================================================================
# STEP 11: Analysis — Radius of Gyration
# =============================================================================
echo ""
echo "--- Step 11: Radius of Gyration Analysis ---"

if [ -f "prod_noPBC.xtc" ] && [ -f "prod.tpr" ]; then
    # Use PBC-corrected trajectory for gyrate
    echo "SOLU" | $SERIAL $GMX gyrate \
        -f prod_noPBC.xtc \
        -s prod.tpr \
        -n index.ndx \
        -o gyrate_single_water.xvg

    echo "  Rg written to gyrate_single_water.xvg"
elif [ -f "prod.xtc" ] && [ -f "prod.tpr" ]; then
    # Fallback: use raw trajectory
    echo "SOLU" | $SERIAL $GMX gyrate \
        -f prod.xtc \
        -s prod.tpr \
        -n index.ndx \
        -o gyrate_single_water.xvg

    echo "  Rg written to gyrate_single_water.xvg (raw trajectory)"
else
    echo "  WARNING: Production files not found — skipping gyrate."
fi

echo ""
echo "================================================================"
echo "  Single Chain Water-Only Simulation COMPLETE"
echo "  Output:        $DIR_MAIN/"
echo "  Clean traj:    $DIR_MAIN/prod_noPBC.xtc"
echo "  Gyrate:        $DIR_MAIN/gyrate_single_water.xvg"
echo "================================================================"

