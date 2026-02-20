import sys
import os

def fix_topology(filename):
    print(f"Fixing topology file: {filename}")
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        cleaned_lines = []
        molecule_section = False
        molecules = []

        for line in lines:
            stripped = line.strip()
            
            # Detect section
            if stripped.startswith('[') and stripped.endswith(']'):
                if 'molecules' in stripped:
                    molecule_section = True
                    cleaned_lines.append(line.rstrip() + '\n')
                    continue
                else:
                    molecule_section = False
            
            if molecule_section:
                # Parse molecule lines
                if stripped and not stripped.startswith(';'):
                    parts = stripped.split()
                    if len(parts) >= 2:
                        name = parts[0]
                        count = int(parts[1])
                        
                        # Check if we already have this molecule
                        found = False
                        for i, (existing_name, existing_count) in enumerate(molecules):
                            if existing_name == name:
                                # Update count (assuming we want the last one or sum? 
                                # Usually solvate appends, so we might want the last one.
                                # But if we have duplicates, usually gmx solvate adds duplicate lines.
                                # Genion fails if there are duplicate SOL lines sometimes.
                                # Let's keep distinct lines unless identical?
                                # Actually, standard is: PROH 10, SOL 1000. 
                                # If we see SOL 1000 and SOL 200, is it additive? Yes.
                                # But genion expects ONE line for SOL to edit.
                                # We will merge SOL entries.
                                if name == 'SOL':
                                    molecules[i] = (name, existing_count + count) # Merge
                                    found = True
                                    break
                        if not found:
                            molecules.append((name, count))
                    else:
                        # Comments or empty in section
                        pass 
                else:
                    # preserve comments/empty lines within section but we rebuild it later
                    pass
            else:
                cleaned_lines.append(line.rstrip() + '\n')

        # Reconstruct [ molecules ] section
        if molecules:
            for name, count in molecules:
                cleaned_lines.append(f"{name:<15} {count}\n")
        
        # Write back with Unix endings
        with open(filename, 'w', newline='\n') as f:
            f.writelines(cleaned_lines)

        print("Topology fixed and normalized.")

    except Exception as e:
        print(f"Error fixing topology: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fix_topology.py <topol.top>")
    else:
        fix_topology(sys.argv[1])
