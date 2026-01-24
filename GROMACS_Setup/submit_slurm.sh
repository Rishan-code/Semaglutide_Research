#!/bin/bash
#SBATCH --job-name=Phenol_Sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8            # Request 8 CPU cores for GROMACS
#SBATCH --partition=phd_student        # User's cluster partition
#SBATCH --time=24:00:00                # Running for 15 ns requires time!
#SBATCH --output=par_phenol_%j.log
#SBATCH --error=par_phenol_%j.err
#SBATCH --mem=8G

# ==========================================
# GROMACS HPC Submission Script
# ==========================================

# 1. Load GROMACS Module
echo "Loading GROMACS module..."
module load gromacs

# 2. Check Permissions
chmod +x Run_Phenol.sh

# 3. Run the Pipeline
echo "Starting Phenol Simulation Pipeline on Host: $(hostname)"
./Run_Phenol.sh

echo "Job Complete."
