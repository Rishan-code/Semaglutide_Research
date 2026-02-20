
input_rtf = r"d:\Ozempic\GROMACS_Setup\New_Simulation\charmm-gui-7152362725\res\res.rtf"
output_rtf = r"d:\Ozempic\GROMACS_Setup\New_Simulation\sem_fixed.rtf"

# Mapping from Ligand Reader names to Standard AA names
atom_map = {
    "N1": "N",
    "CA1": "CA",
    "C2": "C",
    "O11": "O",
    "CB1": "CB",
    "CG1": "CG",
    "CD1": "CD",
    "CE": "CE",
    "NE2": "NZ"
}

with open(input_rtf, 'r') as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    # Handle Residue Name
    if line.strip().startswith("RESI res"):
        new_lines.append(line.replace("RESI res", "RESI SEM"))
        continue
        
    # Handle Atom Definitions and Bond Definitions
    # Simple replace might be dangerous if names are substrings (e.g. C2 vs C21)
    # But CGenFF names are usually unique in the context or we can match whole words.
    
    parts = line.split()
    if not parts:
        new_lines.append(line)
        continue
        
    modified_line = line
    
    # Check for atom names in line and replace
    # We sort keys by length descending to avoid substring issues just in case, though precise token matching is better.
    # But for RTF, atoms are space separated.
    
    for old_name, new_name in atom_map.items():
        # Use regex to replace whole words only to avoid replacing C21 with C1
        import re
        modified_line = re.sub(r'\b' + old_name + r'\b', new_name, modified_line)
        
    new_lines.append(modified_line)

with open(output_rtf, 'w') as f:
    f.writelines(new_lines)

print(f"Fixed RTF saved to {output_rtf}")
