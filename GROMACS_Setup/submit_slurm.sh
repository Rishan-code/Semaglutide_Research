#!/bin/bash
#SBATCH --job-name=Phenol_Sim
#SBATCH --output=par_phenol_%j.log
#SBATCH --error=par_phenol_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --mem=8G

# ==========================================
# GROMACS HPC Submission Script
# ==========================================

# 1. Load GROMACS Module
# ADJUST THIS depending on your cluster!
# Examples:
# module load gromacs/2023.3
# module load gromacs/2023.3-cuda
echo "Loading GROMACS module..."
module load gromacs

# 2. Check Permissions
chmod +x Run_Phenol.sh

# 3. Run the Pipeline
# This script is idempotent:
# - If steps are done, it skips them.
# - If production was interrupted, it resumes (-cpi).
echo "Starting Phenol Simulation Pipeline on Host: $(hostname)"
./Run_Phenol.sh

echo "Job Complete."
