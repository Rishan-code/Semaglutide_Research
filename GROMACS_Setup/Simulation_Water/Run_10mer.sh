#!/bin/bash

# Stop on err
set -e

echo "----------------------------------------------------------------"
echo "Setting up 10-mer Semaglutide Simulation"
echo "----------------------------------------------------------------"

# Define Step Directories
DIR_SOLV="01_solvation"
DIR_IONS="02_ions"
DIR_MIN="03_minimization"
DIR_EQUI="04_equilibration"
DIR_PROD="05_production"

# Create Directories
mkdir -p $DIR_SOLV $DIR_IONS $DIR_MIN $DIR_EQUI $DIR_PROD

# Helper function to move files if they exist in root but not in target
migrate() {
    local file=$1
    local dest_dir=$2
    if [ -f "$file" ] && [ ! -f "$dest_dir/$file" ]; then
        echo "Moving $file to $dest_dir/"
        mv "$file" "$dest_dir/"
    fi
}

# --- 1. Solvation ---
SOLV_GRO="${DIR_SOLV}/system_10mer_solv.gro"
migrate "system_10mer_solv.gro" "$DIR_SOLV"

if [ ! -f "$SOLV_GRO" ]; then
    echo "[1/5] Solvating system..."
    # Note: Inputs (system_10mer.gro, spc216.gro) are expected in root or standard paths
    # Outputting to directory
    gmx solvate -cp system_10mer.gro -cs spc216.gro -o "$SOLV_GRO" -p topol_10mer.top
else
    echo "[1/5] Solvated system already exists (skipping)"
fi


# --- 2. Add Ions ---
FINAL_GRO="${DIR_IONS}/system_10mer_final.gro"
IONS_TPR="${DIR_IONS}/ions.tpr"

migrate "system_10mer_final.gro" "$DIR_IONS"
migrate "ions.tpr" "$DIR_IONS"

if [ ! -f "$FINAL_GRO" ]; then
    echo "[2/5] Adding Ions (Neutralizing + 0.15M)..."
    
    # FIX: Sanitize topology before running genion/grompp
    python3 fix_topology.py topol_10mer.top

    # Generate TPR for genion
    gmx grompp -f ions.mdp -c "$SOLV_GRO" -p topol_10mer.top -o "$IONS_TPR" -maxwarn 10
    
    echo "SOL" | gmx genion -s "$IONS_TPR" -o "$FINAL_GRO" -p topol_10mer.top -pname POT -nname CLA -neutral -conc 0.15
else
    echo "[2/5] Ionized system already exists (skipping)"
fi


# --- Configuration for Run ---
INIT_GRO="$FINAL_GRO" # The output of Step 2 is the input for minimization

# --- 2b. Generate Index File ---
# Index file usually stays in root or can be moved. Let's keep it in root as it's a shared config.
if [ ! -f "index.ndx" ]; then
    echo "[2.5/5] Generating Index File..."
    echo -e "\"Protein\" | \"Other\"\nname 17 SOLU\n! \"SOLU\"\nname 18 SOLV\nq" | gmx make_ndx -f "$FINAL_GRO" -o index.ndx
fi


# --- 3. Minimization ---
MINI_PREFIX="step4.0_minimization_10mer"
MINI_GRO="${DIR_MIN}/${MINI_PREFIX}.gro"
MINI_TPR="${DIR_MIN}/${MINI_PREFIX}.tpr"

# Migration of all potential minimization files
for ext in gro tpr log edr trr cpt; do
    migrate "${MINI_PREFIX}.${ext}" "$DIR_MIN"
done

if [ ! -f "$MINI_GRO" ]; then
    echo "[3/5] Running Energy Minimization..."
    # cwd prevents clutter? No, we use absolute or relative paths.
    # We will run gmx commands pointing to the files in their dir.
    
    gmx grompp -f step4.0_minimization.mdp -o "$MINI_TPR" -c "$INIT_GRO" -r "$INIT_GRO" -p topol_10mer.top -maxwarn 10
    
    # mdrun output handling: we use -deffnm with the path.
    # Note: GROMACS writes aux files to the same dir as -deffnm base.
    gmx mdrun -v -deffnm "${DIR_MIN}/${MINI_PREFIX}" -ntmpi 1
else
    echo "[3/5] Minimization already complete (skipping)"
fi


# --- 4. Equilibration ---
EQUI_PREFIX="step4.1_equilibration_10mer"
EQUI_GRO="${DIR_EQUI}/${EQUI_PREFIX}.gro"
EQUI_TPR="${DIR_EQUI}/${EQUI_PREFIX}.tpr"

# Migration
for ext in gro tpr log edr trr cpt xtc; do
    migrate "${EQUI_PREFIX}.${ext}" "$DIR_EQUI"
    migrate "${EQUI_PREFIX}_prev.${ext}" "$DIR_EQUI" 2>/dev/null || true
done

if [ ! -f "$EQUI_GRO" ]; then
    echo "[4/5] Running Equilibration..."
    
    gmx grompp -f step4.1_equilibration.mdp -o "$EQUI_TPR" -c "$MINI_GRO" -r "$INIT_GRO" -p topol_10mer.top -n index.ndx
    gmx mdrun -v -deffnm "${DIR_EQUI}/${EQUI_PREFIX}" -ntmpi 1
else
    echo "[4/5] Equilibration already complete (skipping)"
fi


# --- 5. Production ---
PROD_PREFIX="step5_production_10mer"
PROD_GRO="${DIR_PROD}/${PROD_PREFIX}.gro"
PROD_TPR="${DIR_PROD}/${PROD_PREFIX}.tpr"

# Migration
for ext in gro tpr log edr trr cpt xtc; do
    migrate "${PROD_PREFIX}.${ext}" "$DIR_PROD"
    migrate "${PROD_PREFIX}_prev.${ext}" "$DIR_PROD" 2>/dev/null || true
done

echo "----------------------------------------------------------------"
echo "Starting Production Pipeline"
echo "----------------------------------------------------------------"

if [ ! -f "$PROD_GRO" ]; then
    echo "[5/5] Running Production Simulation..."
    
    # Check if we should restart from a checkpoint or start fresh
    # If a checkpoint exists in the dir, mdrun -cpi auto-detects it if we provide it.
    
    # Generate TPR if not exists (or always regenerate? Safe to regenerate if inputs didn't change, but keep it simple)
    if [ ! -f "$PROD_TPR" ]; then
         gmx grompp -f step5_production.mdp -o "$PROD_TPR" -c "$EQUI_GRO" -p topol_10mer.top -n index.ndx
    fi
    
    # Run
    gmx mdrun -v -deffnm "${DIR_PROD}/${PROD_PREFIX}" -ntmpi 1 -cpi "${DIR_PROD}/${PROD_PREFIX}.cpt" -append
else
    echo "[5/5] Production already complete (skipping)"
fi

echo "----------------------------------------------------------------"
echo "Simulation Pipeline Complete!"
echo "----------------------------------------------------------------"
