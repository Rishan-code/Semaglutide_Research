
import re

rtf_path = r"d:\Ozempic\GROMACS_Setup\New_Simulation\sem_standardized.rtf"
output_path = r"d:\Ozempic\GROMACS_Setup\New_Simulation\sem_standardized_fixed.rtf"

with open(rtf_path, 'r') as f:
    lines = f.readlines()

atoms = []
total_charge = 0.0
atom_indices = {}

# Parse Charges
for i, line in enumerate(lines):
    if line.strip().startswith("ATOM"):
        parts = line.split()
        # parts[0]="ATOM", parts[1]=Name, parts[2]=Type, parts[3]=Charge
        name = parts[1]
        charge = float(parts[3])
        atoms.append({'name': name, 'charge': charge, 'line_idx': i})
        atom_indices[name] = i
        total_charge += charge

print(f"Current Total Charge: {total_charge:.6f}")

# Target Integer
target_charge = round(total_charge)
diff = target_charge - total_charge
print(f"Target: {target_charge}, Diff: {diff:.6f}")

# Distribute diff to CB and CG (non-polar sidechain start) avoids messing up backbone or polar head
# CB and CG
targets = ['CB', 'CG']
target_indices = [atom_indices[t] for t in targets if t in atom_indices]

if not target_indices:
    # Fallback to just C11/C12?
    # Try generic carbon sidechain atoms
    targets = ['CA2', 'CB2', 'CG2']
    target_indices = [atom_indices[t] for t in targets if t in atom_indices]

if not target_indices:
    print("Error: Could not find target atoms to adjust.")
    exit()

delta = diff / len(target_indices)

new_lines = list(lines)

for idx in target_indices:
    line = lines[idx]
    parts = line.split()
    old_charge = float(parts[3])
    new_charge = old_charge + delta
    
    # Reconstruct line carefully maintaining spacing
    # ATOM <NAME> <TYPE> <CHARGE> ...
    # Standard format: ATOM %-4s %-4s %6.3f
    # Check original spacing
    # pattern: ATOM\s+(\S+)\s+(\S+)\s+(\S+)(.*)
    m = re.match(r"(ATOM\s+\S+\s+\S+\s+)(\S+)(.*)", line)
    if m:
        prefix = m.group(1)
        suffix = m.group(3)
        # Format charge to 3 decimals
        new_lines[idx] = f"{prefix}{new_charge:6.3f}{suffix}\n"
        print(f"Adjusted {parts[1]}: {old_charge} -> {new_charge}")

# Verify
check_sum = 0.0
for i, line in enumerate(new_lines):
    if line.strip().startswith("ATOM"):
        parts = line.split()
        check_sum += float(parts[3])

print(f"New Total Charge: {check_sum:.6f}")

with open(output_path, 'w') as f:
    f.writelines(new_lines)
    
print(f"Fixed RTF written to {output_path}")
