
input_pdb = r"d:\Ozempic\GROMACS_Setup\semaglutide_backbone_clean.pdb"
output_pdb = r"d:\Ozempic\GROMACS_Setup\semaglutide_solution_input.pdb"

with open(input_pdb, 'r') as f:
    lines = f.readlines()

new_lines = []
found_lys26 = False

for line in lines:
    if line.startswith("ATOM"):
        res_seq = int(line[22:26])
        res_name = line[17:20]
        
        if res_seq == 26:
            # Change Residue Name to SEM
            # ATOM      1  N   LYS A  26
            # 17-20 is ResName
            new_line = line[:17] + "SEM" + line[20:]
            
            # Keep standard atom names (N, CA, C, O, CB...) because we mapped RTF to match these!
            # But standard Lysine has NZ relative to CE.
            # Our "SEM" RTF has NZ relative to CE.
            # Perfect.
            
            new_lines.append(new_line)
            found_lys26 = True
        else:
            new_lines.append(line)
    else:
        new_lines.append(line)

with open(output_pdb, 'w') as f:
    f.writelines(new_lines)

print(f"Created {output_pdb} with LYS26 -> SEM mutation")
