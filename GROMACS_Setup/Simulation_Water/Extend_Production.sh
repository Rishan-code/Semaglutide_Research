#!/bin/bash
# Usage: ./Extend_Production.sh <time_in_ps>

EXTEND_TIME=$1

if [ -z "$EXTEND_TIME" ]; then
    echo "----------------------------------------------------------------"
    echo "Usage: ./Extend_Production.sh <time_in_ps>"
    echo ""
    echo "Example: ./Extend_Production.sh 1000"
    echo "         (This extends the run by 1000 ps = 1 ns)"
    echo "----------------------------------------------------------------"
    exit 1
fi

# Define paths (matching your organized structure)
DIR_PROD="05_production"
PREFIX="step5_production_10mer"
TPR="${DIR_PROD}/${PREFIX}.tpr"
CPT="${DIR_PROD}/${PREFIX}.cpt"

if [ ! -f "$TPR" ]; then
    echo "Error: TPR file ($TPR) not found. Are you in the root GROMACS_Setup directory?"
    exit 1
fi

if [ ! -f "$CPT" ]; then
    echo "Error: Checkpoint file ($CPT) not found. Cannot restart."
    exit 1
fi

echo "----------------------------------------------------------------"
echo "Extending simulation by $EXTEND_TIME ps..."
echo "----------------------------------------------------------------"

# Backup the original TPR just in case
cp "$TPR" "${TPR}.backup"

# Extend the TPR file
gmx convert-tpr -s "$TPR" -extend "$EXTEND_TIME" -o "$TPR"

echo "TPX file updated. Resuming simulation..."

# Change to the production directory to avoid path mismatches with the checkpoint
cd "$DIR_PROD"

# Resume using the local filenames (since we are now IN the directory)
# We expect the files to be step5_production_10mer.*
gmx mdrun -v -deffnm "$PREFIX" -cpi "${PREFIX}.cpt" -append -ntmpi 1
