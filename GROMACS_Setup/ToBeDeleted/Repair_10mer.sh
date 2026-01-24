#!/bin/bash
set -e

echo "------------------------------------------------"
echo "Repairing 10-mer System"
echo "------------------------------------------------"

# 1. Extract valid single molecule from the original input
# We use Python to robustly grab the first 507 atoms (Protein + KWF1 Tail)
# bypassing GROMACS group selection issues.
if [ -f "Extra_Files/step3_input.gro" ]; then
    echo "Extracting valid monomer (Protein + Tail)..."
    python3 extract_protein.py "Extra_Files/step3_input.gro" "single_correct.gro"
else
    echo "Error: Extra_Files/step3_input.gro not found!"
    exit 1
fi

# 2. Create new 10-mer box
echo "Creating new system_10mer.gro..."
gmx insert-molecules -ci single_correct.gro -nmol 10 -box 10 10 10 -o system_10mer.gro

# 3. CLEANUP (Critical)
# We must delete the old solvated files so Run_10mer.sh regenerates them!
echo "Removing stale files..."
rm -f system_10mer_solv.gro system_10mer_final.gro ions.tpr

echo "------------------------------------------------"
echo "Repair Complete & Stale Files Removed."
echo "Now run: ./Run_10mer.sh"
echo "------------------------------------------------"
