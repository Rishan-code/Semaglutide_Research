import os
import shutil

def organize():
    root = "."
    
    # 1. Define Directories
    dir_water = os.path.join(root, "Simulation_Water")
    dir_phenol = os.path.join(root, "Simulation_Phenol")
    
    os.makedirs(dir_water, exist_ok=True)
    os.makedirs(dir_phenol, exist_ok=True)
    
    # 2. Move Water Simulation Folders
    # Folders 01-05 and analysis
    folders_to_move = [
        "01_solvation", 
        "02_ions", 
        "03_minimization", 
        "04_equilibration", 
        "05_production", 
        "analysis"
    ]
    
    for folder in folders_to_move:
        src = os.path.join(root, folder)
        dst = os.path.join(dir_water, folder)
        if os.path.exists(src):
            print(f"Moving {folder} -> {dir_water}")
            try:
                shutil.move(src, dst)
            except Exception as e:
                print(f"Skipped {folder}: {e}")

    # 3. Move Water Simulation Files
    # specific files usually in root
    files_to_move = [
        "topol_10mer.top",
        "system_10mer.gro",
        "index.ndx",
        "Run_10mer.sh",
        "Extend_Production.sh",
        "submit_slurm.sh", # Assuming this was for the 10-mer
        "measure_monomer.sh",
        "split_index.input"
    ]
    
    # Also move any lingering logs/xvg/gro from root if they look like 10mer files
    # (But be careful not to move shared scripts or phenol files)
    
    for f in files_to_move:
        src = os.path.join(root, f)
        dst = os.path.join(dir_water, f)
        if os.path.exists(src):
            print(f"Moving {f} -> {dir_water}")
            try:
                shutil.move(src, dst)
            except Exception as e:
                print(f"Skipped {f}: {e}")

    # 4. Handle GROMACS_Phenol (if it was created by previous attempts)
    # The user wants "Simulation_Phenol". If "GROMACS_Phenol" exists, we can rename/merge it.
    old_phenol_dir = os.path.join(root, "GROMACS_Phenol")
    if os.path.exists(old_phenol_dir):
        print(f"Merging {old_phenol_dir} -> {dir_phenol}")
        # Move contents
        for item in os.listdir(old_phenol_dir):
            s = os.path.join(old_phenol_dir, item)
            d = os.path.join(dir_phenol, item)
            if not os.path.exists(d):
                shutil.move(s, d)
        # Remove empty dir
        try:
            os.rmdir(old_phenol_dir)
        except:
            pass
            
    print("Organization Complete.")
    print(f"Water Run -> {dir_water}")
    print(f"Phenol Run -> {dir_phenol}")

if __name__ == "__main__":
    organize()
