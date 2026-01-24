#!/bin/bash
#SBATCH --job-name=Semaglutide_10mer   # Job name
#SBATCH --partition=gpu                # Partition (queue) name - CHECK YOUR CLUSTER DOCS!
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of MPI tasks (usually 1 for GPU runs)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task
#SBATCH --gres=gpu:1                   # Request 1 GPU
#SBATCH --time=24:00:00                # Time limit (HH:MM:SS)
#SBATCH --output=slurm_%j.out          # Standard output log
#SBATCH --error=slurm_%j.err           # Standard error log

# Load GROMACS module (Check your cluster specific module name!)
# module load gromacs/2023.3-gpu
# OR
# source /path/to/GMXRC

echo "Job started on $(hostname) at $(date)"

# Define paths
DIR_PROD="05_production"
PREFIX="step5_production_10mer"

# Ensure we are in the right place
if [ ! -d "$DIR_PROD" ]; then
    echo "Error: Directory $DIR_PROD not found."
    exit 1
fi

cd "$DIR_PROD" || exit

# Run GROMACS
# Note: -ntmpi 1 for 1 MPI task (common for single GPU runs)
# -ntomp $SLURM_CPUS_PER_TASK uses the requested CPU cores for OpenMP
gmx mdrun -v -deffnm "$PREFIX" -cpi "${PREFIX}.cpt" -append -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK}

echo "Job finished at $(date)"
