# GROMACS Simulation Guide: Semaglutide + Phenol

## Current Status
- **Goal:** Simulate self-assembly of 10-mer Semaglutide with Phenol.
- **Target Time:** 100 ns (needed for convergence).
- **Current Progress:** ~15 ns (job stopped due to 24h walltime limit).

## How to Continue the Simulation
Your scripts are configured to **automatically resume** from the last checkpoint. You do not need to edit any files.

### 1. Submit the Job on the Cluster
Run this command on your HPC cluster:
```bash
sbatch submit_slurm.sh
```
*   This will run for another 24 hours.
*   It will extend the simulation by approximately **15 ns** per run.
*   **Repeat this process** about 5-6 times until you reach 100 ns.

### 2. Check Progress
To see how many picoseconds have been simulated so far:
```bash
gmx check -f Simulation_Phenol/step5_production_phenol.xtc
```

### 3. Analyze Convergence (Radius of Gyration)
After each chunk completes, update the Rg plot to see if the structure has stabilized (flattened out).
```bash
cd Simulation_Phenol
gmx gyrate -s step5_production_phenol.tpr -f step5_production_phenol.xtc -o gyrate_phenol.xvg
# (Select group 1 "Protein" when prompted)
```

## Troubleshooting
**Job Failures:**
- If the job crashes, check the error logs (`par_phenol_*.err`).
- If you need to restart from scratch (reset everything), delete the `Simulation_Phenol` directory and run `./Run_Phenol.sh` manually.

**Convergence:**
- If Rg is still fluctuating after 100 ns, you may need to extend even further (e.g., to 200 ns).
