"""Move useless files/dirs from GROMACS_Setup root into Trash/ folder."""
import shutil
import os

BASE = os.path.dirname(os.path.abspath(__file__))
TRASH = os.path.join(BASE, "Trash")

# Ensure Trash directory exists
os.makedirs(TRASH, exist_ok=True)

FILES_TO_MOVE = [
    # One-off Python scripts
    "analyze_itp.py",
    "append_bridge_params.py",
    "cleanup.py",
    "extract_protein.py",
    "extract_semaglutide.py",
    "fix_charge_sum.py",
    "fix_res_rtf.py",
    "fix_topology.py",
    "graft_sidechain.py",
    "gro_to_mol2.py",
    "gro_to_pdb.py",
    "organize_files.py",
    "organize_projects.py",
    "prepare_clean_pdb.py",
    "prepare_solution_input.py",
    "rename_phenol.py",
    "standardize_rtf.py",
    "translate_params.py",
    # Intermediate/unused data files
    "modified_lysine_smiles.txt",
    "semaglutide_backbone_clean.pdb",
    "semaglutide_extracted.pdb",
    "semaglutide_merged.pdb",
    "semaglutide_solution_input.pdb",
    "7K10.pdb",
    "7KI0.pdb",
    # Duplicate / obsolete
    "submit_slurm.sh",
    "SIMULATION_README.md",
]

DIRS_TO_MOVE = [
    "Simulation_Single_Phenol",
    "New_Simulation",
]

moved = 0
for name in FILES_TO_MOVE:
    src = os.path.join(BASE, name)
    dst = os.path.join(TRASH, name)
    if os.path.exists(src):
        shutil.move(src, dst)
        print(f"  Moved: {name}")
        moved += 1
    else:
        print(f"  SKIP (not found): {name}")

for name in DIRS_TO_MOVE:
    src = os.path.join(BASE, name)
    dst = os.path.join(TRASH, name)
    if os.path.exists(src):
        shutil.move(src, dst)
        print(f"  Moved dir: {name}/")
        moved += 1
    else:
        print(f"  SKIP dir (not found): {name}/")

print(f"\nDone! Moved {moved} items to Trash/")
