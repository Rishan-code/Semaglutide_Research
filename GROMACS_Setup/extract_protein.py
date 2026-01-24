import sys

def extract_protein(input_file, output_file, atom_count=507):
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < atom_count + 3:
            print(f"Error: Input file has fewer lines than expected ({len(lines)})")
            sys.exit(1)

        # 1. Header
        title = lines[0]
        # 2. Atoms (Lines 2 to 2+atom_count-1 -> indices 2 to atom_count+2)
        # GRO input:
        # Line 0: Title
        # Line 1: Num Atoms
        # Line 2: First Atom
        
        atoms = lines[2 : 2 + atom_count]
        
        # 3. Footer (Box vectors - usually the last line)
        box = lines[-1]

        # Verify we got the right number
        if len(atoms) != atom_count:
            print(f"Error during extraction: Captured {len(atoms)} atoms, expected {atom_count}")
            sys.exit(1)

        with open(output_file, 'w') as f:
            f.write("Semaglutide Monomer Corrected\n")
            f.write(f"{atom_count}\n")
            for line in atoms:
                f.write(line)
            f.write(box)
        
        print(f"Successfully extracted {atom_count} atoms to {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_protein.py <input.gro> <output.gro>")
        sys.exit(1)
    
    extract_protein(sys.argv[1], sys.argv[2])
