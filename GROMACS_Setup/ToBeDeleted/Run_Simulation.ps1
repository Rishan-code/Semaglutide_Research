$ErrorActionPreference = "Stop"

# Configuration
$init = "step3_input"
$mini_prefix = "step4.0_minimization"
$equi_prefix = "step4.1_equilibration"
$prod_prefix = "step5_production"

Write-Host "Starting Semaglutide Simulation on Windows..." -ForegroundColor Cyan

# 1. Minimization
Write-Host "1. Energy Minimization..." -ForegroundColor Yellow
gmx grompp -f "${mini_prefix}.mdp" -o "${mini_prefix}.tpr" -c "${init}.gro" -r "${init}.gro" -p topol.top -n index.ndx -maxwarn 10
if ($LASTEXITCODE -ne 0) { Write-Error "Minimization grompp failed"; exit 1 }

gmx mdrun -v -deffnm $mini_prefix -ntmpi 1
if ($LASTEXITCODE -ne 0) { Write-Error "Minimization mdrun failed"; exit 1 }

# 2. Equilibration
Write-Host "2. Equilibration..." -ForegroundColor Yellow
gmx grompp -f "${equi_prefix}.mdp" -o "${equi_prefix}.tpr" -c "${mini_prefix}.gro" -r "${init}.gro" -p topol.top -n index.ndx
if ($LASTEXITCODE -ne 0) { Write-Error "Equilibration grompp failed"; exit 1 }

gmx mdrun -v -deffnm $equi_prefix -ntmpi 1
if ($LASTEXITCODE -ne 0) { Write-Error "Equilibration mdrun failed"; exit 1 }

# 3. Production (Running just step 1 for now)
Write-Host "3. Production Run..." -ForegroundColor Yellow
$prod_step = "${prod_prefix}_1"
gmx grompp -f "${prod_prefix}.mdp" -o "${prod_step}.tpr" -c "${equi_prefix}.gro" -p topol.top -n index.ndx
if ($LASTEXITCODE -ne 0) { Write-Error "Production grompp failed"; exit 1 }

gmx mdrun -v -deffnm $prod_step -ntmpi 1

Write-Host "Simulation Complete!" -ForegroundColor Green
