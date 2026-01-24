import os
import shutil
import glob

TO_BE_DELETED = "ToBeDeleted"
EXTRA_FILES_DIR = "Extra_Files"

# Files to KEEP in root (besides the main folders)
KEEP_EXACT = {
    "Run_10mer.sh", 
    "step5_production.mdp", 
    "topol_10mer.top", 
    "index.ndx",
    "fix_topology.py",
    "extract_protein.py",
    "organize_files.py",
    "cleanup.py",
    "ions.mdp",
    "step4.0_minimization.mdp",
    "step4.1_equilibration.mdp",
    "system_10mer.gro" # Keep initial input? User said "put all unwanted". Let's assume input gro is wanted.
}

KEEP_DIRS = {
    "01_solvation",
    "02_ions",
    "03_minimization",
    "04_equilibration",
    "05_production",
    "toppar",
    TO_BE_DELETED # obviously
}

def cleanup():
    if not os.path.exists(TO_BE_DELETED):
        os.makedirs(TO_BE_DELETED)
        print(f"Created {TO_BE_DELETED}")

    # 1. Move Extra_Files content
    if os.path.exists(EXTRA_FILES_DIR):
        print(f"Moving content from {EXTRA_FILES_DIR}...")
        for item in os.listdir(EXTRA_FILES_DIR):
            src = os.path.join(EXTRA_FILES_DIR, item)
            dst = os.path.join(TO_BE_DELETED, item)
            
            # If valid usage of cleanup, handle duplicates by renaming
            if os.path.exists(dst):
                base, ext = os.path.splitext(item)
                dst = os.path.join(TO_BE_DELETED, f"{base}_copy{ext}")

            try:
                shutil.move(src, dst)
            except Exception as e:
                print(f"Error moving {item}: {e}")
        
        # Move the empty directory itself? Or leave it?
        # User said "put files from the folder extra files".
        # I'll move the folder itself if empty, or just move the folder INTO ToBeDeleted?
        # "put files from the folder extra files into that new folder" -> implied flat or nested?
        # Flat is messier if collisions. Nested is cleaner.
        # But user said "unwanted files into a folder".
        # I already moved them flat above.
        
        try:
            os.rmdir(EXTRA_FILES_DIR)
            print("Removed empty Extra_Files directory")
        except:
             print("Extra_Files not empty, keeping it")

    # 2. Move files from Root
    print("Cleaning root directory...")
    for item in os.listdir("."):
        if item in KEEP_EXACT or item in KEEP_DIRS:
            continue
        
        # Files like #...#, temp.top..., etc.
        # Also .mdp files? I added mdps to KEEP_EXACT just in case, 
        # but user might consider them "unwanted" if they are copies?
        # "dont touch the solvation ions minimization and etc folders"
        
        # Safest bet: Move everything NOT in the explicit Keep list.
        # Let's review the list.
        
        # Run_10mer_BACKUP.sh -> Move
        # Repair_10mer.sh -> Move? Probably (was from prev attempt)
        # mdout.mdp -> Move (generated)
        # single_correct.gro -> Move?
        # temp.top* -> Move
        # #*# -> Move
        
        src = os.path.join(".", item)
        dst = os.path.join(TO_BE_DELETED, item)
        
        if os.path.exists(dst):
             base, ext = os.path.splitext(item)
             dst = os.path.join(TO_BE_DELETED, f"{base}_root_copy{ext}")

        try:
            shutil.move(src, dst)
            print(f"Moved {item}")
        except Exception as e:
             print(f"Error moving {item}: {e}")

if __name__ == "__main__":
    cleanup()
