import os
import shutil
import glob

TRASH_DIR = "Trash"
SIM_DIRS = ["Simulation_Water", "Simulation_Phenol", "Simulation_Single_Phenol", "."]

# Patterns to MOVE to Trash (Glob patterns)
PATTERNS_TO_MOVE = [
    "#*#",              # GROMACS backups
    "step4.0_*",        # Minimization files
    "step4.1_*",        # Equilibration files
    "*.log",            # Logs (except production if whitelisted)
    "*.edr",            # Energies
    "*.trr",            # Trajectories
    "*.cpt",            # Checkpoints
    "ions.tpr",
    "mdout.mdp",
    "*.xvg",            # Graphs
    "box_*.gro",        # Intermediate build steps
    "protein_only.gro",
    "ions.gro",
    "temp_*.top",
]

# Specific files to KEEP (Whitelist)
FILES_TO_KEEP = [
    # Scripts
    "Run_Phenol.sh",
    "Run_Water.sh",
    "Run_Single_Phenol.sh",
    "submit_slurm.sh",
    "cleanup.py",
    "extract_protein.py",
    "fix_topology.py",
    "gro_to_mol2.py",
    "gro_to_pdb.py",
    "organize_files.py",
    "organize_projects.py",
    "rename_phenol.py",
    "test_installation.py",
    
    # Templates & Data
    "phenol_charmm.itp",
    "phenol_charmm.pdb",
    "ions.mdp",
    "step4.0_minimization.mdp",
    "step4.1_equilibration.mdp",
    "step5_production.mdp",
    
    # Docs
    "README.md",
    "SIMULATION_README.md",
    ".gitignore",
    
    # Directories (Keep them!)
    "Simulation_Phenol",
    "Simulation_Water",
    "Simulation_Single_Phenol",
    "toppar",
    "ToBeDeleted",
    "Trash",
    "__pycache__",
    ".git",
    ".vs",
    
    # Production Outputs (Sim_Phenol)
    "step5_production_phenol.xtc",
    "step5_production_phenol.tpr",
    "step5_production_phenol.gro",
    "step5_production_phenol.edr",
    "step5_production_phenol.log",
    "step5_production_phenol.cpt",
    "index_phenol.ndx",
    "topol_phenol.top",
    "step5_production_phenol_clean.xtc",
    
    # Production Outputs (Sim_Water)
    "step5_production_10mer.xtc",
    "step5_production_10mer.tpr",
    "step5_production_10mer.gro",
    "step5_production_10mer.edr",
    "step5_production_10mer.log",
    "step5_production_10mer.cpt",
    "index.ndx",
    "topol_10mer.top",
]

def cleanup_directory(directory):
    trash_path = os.path.join(directory, TRASH_DIR)
    if not os.path.exists(trash_path):
        os.makedirs(trash_path)
        # print(f"Created {trash_path}") 

    print(f"Cleaning {directory}...")
    
    # Get all files in the directory
    all_files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    
    for filename in all_files:
        filepath = os.path.join(directory, filename)
        
        # Check against whitelist FIRST
        if filename in FILES_TO_KEEP:
            continue
            
        move_it = False
        
        # Check against patterns
        for pattern in PATTERNS_TO_MOVE:
            if glob.fnmatch.fnmatch(filename, pattern):
                move_it = True
                break
        
        if move_it:
            dst = os.path.join(trash_path, filename)
            # Handle duplicates
            if os.path.exists(dst):
                base, ext = os.path.splitext(filename)
                count = 1
                while os.path.exists(os.path.join(trash_path, f"{base}_{count}{ext}")):
                    count += 1
                dst = os.path.join(trash_path, f"{base}_{count}{ext}")
            
            try:
                shutil.move(filepath, dst)
                print(f"Moved {filename} -> {TRASH_DIR}/")
            except Exception as e:
                print(f"Error moving {filename}: {e}")

def main():
    for sim_dir in SIM_DIRS:
        if os.path.exists(sim_dir):
            cleanup_directory(sim_dir)
        else:
            print(f"Directory {sim_dir} not found (skipping).")

if __name__ == "__main__":
    main()
