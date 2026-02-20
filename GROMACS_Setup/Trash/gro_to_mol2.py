import sys

def write_mol2(gro_file, mol2_file):
    with open(gro_file, 'r') as f:
        lines = f.readlines()

    atom_count = int(lines[1])
    atoms = []
    
    # Parse atoms
    # GRO: resid(5) resname(5) atomname(5) atomid(5) x(8.3) y(8.3) z(8.3)
    # Mapping logic for Phenol (C6H5OH)
    # Based on standard ordering in the input GRO (likely from LigParGen):
    # 1: O00 (Hydroxyl Oxygen) -> O.3
    # 2: C01 (Attached to O)   -> C.ar
    # 3-7: C02-C06 (Ring)      -> C.ar
    # 8-13: Hydrogens          -> H
    
    for i in range(2, 2 + atom_count):
        line = lines[i]
        name = line[10:15].strip()
        x = float(line[20:28]) * 10.0 # nm to Angstrom
        y = float(line[28:36]) * 10.0
        z = float(line[36:44]) * 10.0
        
        # Determine Atom Type (Tripos)
        atype = "H"
        if name.startswith("O"):
            atype = "O.3"
        elif name.startswith("C"):
            atype = "C.ar"
            
        atoms.append({'id': i-1, 'name': name, 'x': x, 'y': y, 'z': z, 'type': atype, 'charge': 0.0})

    # Bond Connectivity (Explicit)
    # 1(O00) - 2(C01)
    # 1(O00) - 8(H07) (Assuming H07 is the hydroxyl H)
    # Ring: 2-3-4-5-6-7-2
    # Hydrogens on Ring:
    # 3-9, 4-10, 5-11, 6-12, 7-13
    
    bonds = [
        (1, 2, "1"),  # O - C (Single)
        (1, 8, "1"),  # O - H (Single)
        (2, 3, "ar"), (3, 4, "ar"), (4, 5, "ar"), (5, 6, "ar"), (6, 7, "ar"), (7, 2, "ar"), # Ring (Aromatic)
        (3, 9, "1"), (4, 10, "1"), (5, 11, "1"), (6, 12, "1"), (7, 13, "1") # C - H (Single)
    ]
    
    # Write MOL2
    with open(mol2_file, 'w') as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write("PHENOL\n")
        f.write(f"{len(atoms)} {len(bonds)} 0 0 0\n")
        f.write("SMALL\n")
        f.write("USER_CHARGES\n")
        f.write("\n")
        
        f.write("@<TRIPOS>ATOM\n")
        for a in atoms:
            # atom_id atom_name x y z atom_type subst_id subst_name charge
            f.write(f"{a['id']:4d} {a['name']:<4s} {a['x']:10.4f} {a['y']:10.4f} {a['z']:10.4f} {a['type']:<5s} 1 UNK {a['charge']:.4f}\n")
            
        f.write("@<TRIPOS>BOND\n")
        bid = 1
        for b1, b2, btype in bonds:
            f.write(f"{bid:4d} {b1:4d} {b2:4d} {btype}\n")
            bid += 1

    print(f"Written {mol2_file}")

if __name__ == "__main__":
    write_mol2("phenol.gro", "phenol_SwissParam.mol2")
