import sys

def clean_pdb(input_gro_path, output_pdb_path):
    # Read GRO
    lines = []
    with open(input_gro_path, 'r') as f:
        lines = f.readlines()
        
    print(f"Reading from {input_gro_path} ({len(lines)} lines)")
    
    # Extract only the first chain (assume residues 1-31 or molecule 1)
    # GRO Format: ResidueNum ResidueName AtomName AtomNum X Y Z
    # We will convert KWF1 -> LYS and only keep backbone + CB sidechain atoms that match standard LYS
    
    # Standard LYS atoms (PDB names): N, CA, C, O, CB, CG, CD, CE, NZ
    # KWF1 has N, CA, C, O... and many others.
    # We will map KWF1 atoms to LYS names if they match standard pattern (N, CA, C, O are same)
    # For sidechain beyond CB, we might just truncate if names don't match, or map if possible.
    # Better: keep only N, CA, C, O, CB. Let CHARMM-GUI add the rest.
    
    valid_atoms = {
        'N': 'N', 'CA': 'CA', 'C': 'C', 'O': 'O', 
        'CB': 'CB', 'CG': 'CG', 'CD': 'CD', 'CE': 'CE', 'NZ': 'NZ',
        'H': 'H', 'HA': 'HA' 
        # Add more H if needed
    }
    
    pdb_lines = []
    atom_serial = 1
    
    residue_map = {
        'KWF1': 'LYS',
        'AIB': 'AIB', # Keep AIB if standard, might need ALA if CHARMM doesn't know AIB (it does ideally)
        # Actually AIB is standard-ish.
    }
    
    current_res_id = -1
    res_counter = 0
    
    print("Processing atoms...")
    
    for line in lines[2:-1]: # Skip title/natoms and box
        try:
            res_num = int(line[0:5].strip())
            res_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            x = float(line[20:28]) * 10.0 # nm to Angstrom
            y = float(line[28:36]) * 10.0
            z = float(line[36:44]) * 10.0
        except ValueError:
            continue
            
        # Check if new molecule (residue number resets or jumps)
        if res_num < current_res_id:
            # We only want the first chain (residues 1-31 approx)
             break
        
        current_res_id = res_num
        
        # Map residue name
        final_res_name = residue_map.get(res_name, res_name)
        
        # Filter KWF1 atoms
        if res_name == 'KWF1':
            # Map atom names if possible
            if atom_name in valid_atoms:
                final_atom_name = valid_atoms[atom_name]
            elif atom_name == 'N2': # Might be mapped to NZ? No, NZ is the end of Lysine.
                # In KWF1, N2 might be the epsilon nitrogen if defined that way.
                # Let's check our previous analysis...
                # KWF1 atoms: N, CA, C, O, CB, ... N2(271)??
                # Wait, typically standard Lysine N is epsilon.
                pass
            else:
                 # Skip non-standard atoms for now to let CHARMM rebuild
                 if atom_name not in ['N', 'CA', 'C', 'O', 'CB']:
                     # We keep backbone and beta carbon to anchor
                     continue

        # Format PDB Line
        # ATOM      1  N   HIS A   1      49.462  64.281  30.436  1.00  0.00           N
        # 0-6: Record name "ATOM  "
        # 6-11: Serial
        # 12-16: Atom Name
        # 17-20: Residue Name
        # 21: Chain ID (A)
        # 22-26: Residue Seq
        # 30-38: X
        # 38-46: Y
        # 46-54: Z
        
        # Pad atom name properly (e.g. " CA " vs "N   ")
        if len(atom_name) < 4:
            clean_atom_name = f" {atom_name:<3}" 
        else:
            clean_atom_name = atom_name[:4]
            
        pdb_line = f"ATOM  {atom_serial:>5} {clean_atom_name:<4} {final_res_name:>3} A{res_num:>4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {atom_name[0]}"
        pdb_lines.append(pdb_line)
        atom_serial += 1
        
    # Write output
    with open(output_pdb_path, 'w') as f:
        f.write("\n".join(pdb_lines))
        f.write("\nTER\nEND\n")
        
    print(f"Wrote {len(pdb_lines)} atoms to {output_pdb_path}")

clean_pdb("Simulation_Water/05_production/step5_production_10mer.gro", "semaglutide_backbone_clean.pdb")
