# Semaglutide Research Simulations

This repository contains GROMACS simulation setups for studying the aggregation behavior of Semaglutide (10-mer) in different environments.

## Project Structure

*   **`Simulation_Water/`**: Contains the Phase 1 simulation of 10-mer Semaglutide in pure water.
*   **`Simulation_Phenol/`**: Contains the Phase 2 simulation of 10-mer Semaglutide in a Phenol preservative formulation.
*   **`toppar/`**: CHARMM36 force field parameters.
*   **`inputs/`**: (If created) Simulation parameter files.
*   **`scripts/`**: Analysis and utility scripts.

## Prerequisites

*   **OS**: Linux (Ubuntu 20.04/22.04 via WSL2 recommended on Windows).
*   **Software**: GROMACS (v2023 or later), Python 3.
*   **Python Libraries**: `matplotlib`, `numpy`, `scipy`.

## How to Run

### 1. Phenol Simulation (Phase 2)
This script sets up the system with ~50 Phenol molecules, minimizes, equilibrates, and runs production (10 ns).

```bash
./Run_Phenol.sh
```

**What it does:**
1.  Extracts the protein structure from the Water simulation.
2.  Inserts 50 Phenol molecules.
3.  Generates a CHARMM36-compatible topology.
4.  Solvates and adds ions (Neutral + 0.15M KCl/NaCl).
5.  Runs Energy Minimization.
6.  Runs Equilibration (125 ps, using custom SOLU/SOLV groups).
7.  Runs Production (10 ns).

### 2. Water Simulation (Phase 1)
To re-run the control simulation:

```bash
./Run_10mer.sh
```

## Analysis

### Radius of Gyration (Rg)
To calculate and plot the Radius of Gyration during the simulation:

1.  **Calculate:**
    ```bash
    cd Simulation_Phenol
    gmx gyrate -f step5_production_phenol.xtc -s step5_production_phenol.tpr -o gyrate_phenol.xvg
    # Select Group 1 (Protein) when prompted
    ```

2.  **Plot:**
    ```bash
    python ../plot_gyrate.py
    ```

## Notes
*   **Resume Capability**: The `Run_Phenol.sh` script is idempotent. If it stops, run it again—it will skip completed steps and resume the simulation from the last checkpoint.
*   **Topology**: Phenol topology (`phenol_charmm.itp`) was manually parameterized for CHARMM36 based on Tyrosine sidechain analogs.
