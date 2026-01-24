#!/bin/bash
# measure_monomer.sh

# 1. Create index file with split residues
# We assume 10 chains. 'splitres 1' splits group 1 (Protein) into chains.
# We need to find which group is Protein. Usually group 1.
echo -e "splitres 1\nq" | gmx make_ndx -f 05_production/step5_production_10mer.tpr -o analysis/single_peptide.ndx

# 2. Run gyrate on the FIRST chain (Protein_1)
# The groups created will be named Protein_1, Protein_2 ... Protein_10
# We need to select one. Let's select "Protein_1".
# To find the group number of Protein_1, we can just pass the name if gmx supports it,
# or we have to know the number. make_ndx usually appends them at the end.
# If we have ~18 default groups, Protein_1 might be 19.
# Let's rely on the interactive selection string matching.
# We will echo the GROUP NAME.

echo "Calculating Rg for Protein_1..."
echo "Protein_1" | gmx gyrate -f 05_production/step5_production_10mer.xtc -s 05_production/step5_production_10mer.tpr -n analysis/single_peptide.ndx -o analysis/gyrate_monomer.xvg

echo "Done. Results in analysis/gyrate_monomer.xvg"
