import sys

def rename_phenol_pdb(input_pdb, output_pdb):
    """
    Renames atoms in LigParGen generated Phenol PDB to match 
    CHARMM36 Tyrosine-like naming scheme used in phenol_charmm.itp.
    """
    
    # Mapping Dictionary
    # LigParGen (OPLS) -> CHARMM36 (Manual ITP)
    # Based on standard ring ordering.
    # O00 -> OH1
    # C01 -> CZ  (Attached to Oxygen)
    # C02 -> CE1
    # C03 -> CD1
    # C04 -> CG  (Para to CZ)
    # C05 -> CD2
    # C06 -> CE2
    # Hydrogens
    # H07 -> H   (Hydroxyl H)
    # H08 -> HE1
    # H09 -> HD1
    # H0A -> HG1 (Attached to CG)
    # H0B -> HD2
    # H0C -> HE2
    
    mapping = {
        'O00': 'OH1',
        'C01': 'CZ',
        'C02': 'CE1',
        'C03': 'CD1',
        'C04': 'CG',
        'C05': 'CD2',
        'C06': 'CE2',
        'H07': 'H',
        'H08': 'HE1',
        'H09': 'HD1',
        'H0A': 'HG1',
        'H0B': 'HD2',
        'H0C': 'HE2'
    }

    with open(input_pdb, 'r') as fin, open(output_pdb, 'w') as fout:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # PDB Atom Name columns: 12-16
                old_name = line[12:16].strip()
                
                # Check if we have a mapping
                new_name = mapping.get(old_name, old_name)
                
                # Residue Name Rename (UNK -> PHEN)
                # Columns 17-20
                
                # Reconstruct line
                # 0-12: Start
                # 12-16: Atom Name (4 chars, usually space padded. " OH1" or " CZ ")
                # 16-17: AltLoc
                # 17-20: ResName
                
                # Formatting atom name: Center or left align? PDB standard is specific.
                # 4-character names usually start at col 13. 
                # " OH1"
                fmt_name = f"{new_name:<4s}" if len(new_name) == 4 else f" {new_name:<3s}"
                
                # Replace Residue Name "UNK" with "PHEN"
                line_pre = line[:12]
                line_post_name = line[16:].replace("UNK", "PHEN")
                
                new_line = f"{line_pre}{fmt_name}{line_post_name}"
                fout.write(new_line)
            else:
                fout.write(line)

    print(f"Mapped {input_pdb} to {output_pdb}")

if __name__ == "__main__":
    # Use the PDB we generated earlier
    rename_phenol_pdb("phenol_SwissParam.pdb", "phenol_charmm.pdb")
