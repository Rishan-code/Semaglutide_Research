import os
import shutil

# mapping of directory to file patterns or exact names
STRUCTURE = {
    "01_solvation": ["system_10mer_solv.gro"],
    "02_ions": ["system_10mer_final.gro", "ions.tpr"],
    "03_minimization": ["step4.0_minimization_10mer*"],
    "04_equilibration": ["step4.1_equilibration_10mer*"],
    "05_production": ["step5_production_10mer*"]
}

def organize():
    base_dir = os.getcwd()
    print(f"Organizing files in {base_dir}")

    for folder, patterns in STRUCTURE.items():
        # Create folder
        if not os.path.exists(folder):
            try:
                os.makedirs(folder)
                print(f"Created {folder}")
            except Exception as e:
                print(f"Error creating {folder}: {e}")
                continue
        
        # Move files
        for pattern in patterns:
            # Simple glob matching if * is present
            files = []
            if "*" in pattern:
                import glob
                files = glob.glob(pattern)
            else:
                if os.path.exists(pattern):
                    files = [pattern]
            
            for f in files:
                # check if file is already in destination (avoid self-move loop if logic is wrong)
                src = os.path.join(base_dir, f)
                dst = os.path.join(base_dir, folder, f)
                
                if os.path.exists(src) and not os.path.exists(dst):
                    print(f"Moving {f} -> {folder}/")
                    try:
                        shutil.move(src, dst)
                    except Exception as e:
                        print(f"Failed to move {f}: {e}")

if __name__ == "__main__":
    organize()
