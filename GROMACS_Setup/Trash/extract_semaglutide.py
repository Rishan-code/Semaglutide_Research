
import sys

def extract_chain_p(input_pdb, output_pdb):
    lines = []
    with open(input_pdb, 'r') as f:
        lines = f.readlines()

    chain_p_lines = []
    conect_lines = []
    
    # New serial number tracker
    new_serial = 1
    # Mapping old serial -> new serial
    atom_serial_map = {}
    
    # Residue offset: 7 -> 1
    offset = -6
    
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain_id = line[21]
            if chain_id == 'P':
                try:
                    res_seq_str = line[22:26].strip()
                    res_seq = int(res_seq_str)
                except ValueError:
                    continue
                
                # Filter for Residues 7-36 and WF1 101
                if (7 <= res_seq <= 36) or (res_seq == 101):
                    # Keep chain ID as P? Or change to X? Usually fine as P.
                    
                    # Renumber residue
                    if res_seq == 101:
                        new_res_seq = 101
                    else:
                        new_res_seq = res_seq + offset
                    
                    old_serial = int(line[6:11])
                    atom_serial_map[old_serial] = new_serial
                    
                    # Reconstruct new line using list conversion to handle strict positioning
                    # This is safer than f-string for fixed width
                    line_list = list(line)
                    
                    # Update Serial (columns 7-11, indices 6-11)
                    # Serial is right justified
                    serial_str = f"{new_serial:5d}"
                    line_list[6:11] = list(serial_str)
                    
                    # Update Residue Sequence (columns 23-26, indices 22-26)
                    # ResSeq is right justified
                    res_seq_str = f"{new_res_seq:4d}"
                    line_list[22:26] = list(res_seq_str)
                    
                    new_line = "".join(line_list)

                    
                    chain_p_lines.append(new_line)
                    new_serial += 1
        
        elif line.startswith("CONECT"):
            parts = line.split()
            if len(parts) < 2: continue
            
            try:
                center_serial = int(parts[1])
            except ValueError:
                continue
            
            if center_serial in atom_serial_map:
                new_center = atom_serial_map[center_serial]
                new_others = []
                for other in parts[2:]:
                    try:
                        other_serial = int(other)
                        if other_serial in atom_serial_map:
                            new_others.append(atom_serial_map[other_serial])
                    except ValueError:
                        pass
                
                if new_others:
                    # Construct CONECT line
                    # Usually "CONECT" then 5 chars per serial
                    conect_str = f"CONECT{new_center:5d}"
                    for neighbor in new_others:
                        conect_str += f"{neighbor:5d}"
                    conect_lines.append(conect_str + "\n")

    with open(output_pdb, 'w') as out:
        out.writelines(chain_p_lines)
        out.writelines(conect_lines)
        out.write("END\n")

if __name__ == "__main__":
    extract_chain_p("d:/Ozempic/GROMACS_Setup/7KI0.pdb", "d:/Ozempic/GROMACS_Setup/semaglutide_extracted.pdb")
