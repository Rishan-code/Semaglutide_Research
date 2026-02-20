
import math

def parse_pdb(file_path):
    atoms = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = {
                    "line": line,
                    "serial": int(line[6:11]),
                    "name": line[12:16].strip(),
                    "resName": line[17:20].strip(),
                    "chain": line[21],
                    "resSeq": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": line[76:78].strip()
                }
                atoms.append(atom)
    return atoms

def distance(a, b):
    return math.sqrt((a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2)

def subtract(a, b):
    return [a['x']-b['x'], a['y']-b['y'], a['z']-b['z']]

def add(v1, v2):
    return [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]

def dot(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def cross(v1, v2):
    return [
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    ]

def normalize(v):
    l = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if l == 0: return [0,0,0]
    return [v[0]/l, v[1]/l, v[2]/l]

def mat_vec_mult(m, v):
    return [
        m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
        m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
        m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]
    ]

# Kabsch Algorithm implementation specifically for 3 points
def get_rotation_matrix(p_coords, q_coords):
    # p_coords: target (protein)
    # q_coords: source (ligand)
    # Centroids
    pc = [sum(c[0] for c in p_coords)/3, sum(c[1] for c in p_coords)/3, sum(c[2] for c in p_coords)/3]
    qc = [sum(c[0] for c in q_coords)/3, sum(c[1] for c in q_coords)/3, sum(c[2] for c in q_coords)/3]
    
    # Center points
    p_centered = [[c[0]-pc[0], c[1]-pc[1], c[2]-pc[2]] for c in p_coords]
    q_centered = [[c[0]-qc[0], c[1]-qc[1], c[2]-qc[2]] for c in q_coords]
    
    # Covariance matrix H
    H = [[0,0,0], [0,0,0], [0,0,0]]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                H[j][k] += q_centered[i][j] * p_centered[i][k]
                
    # SVD of H (Approximation or use numpy if available, but doing explicit for safety)
    # Since it's just 3 points, we can construct basis vectors.
    
    # Alternative: Construct basis vectors for P and Q
    # v1 = p2 - p1, v2 = p3 - p1
    # n = v1 x v2
    # Construct orthonormal basis U = [u1, u2, u3]
    
    def get_basis(coords):
        v1 = [coords[1][0]-coords[0][0], coords[1][1]-coords[0][1], coords[1][2]-coords[0][2]] # CA - N
        v2 = [coords[2][0]-coords[1][0], coords[2][1]-coords[1][1], coords[2][2]-coords[1][2]] # C - CA
        u1 = normalize(v1)
        
        # u2 = normalize(v2 - proj_u1(v2))
        d = dot(v2, u1)
        proj = [d*u1[0], d*u1[1], d*u1[2]]
        u2_raw = [v2[0]-proj[0], v2[1]-proj[1], v2[2]-proj[2]]
        u2 = normalize(u2_raw)
        
        u3 = cross(u1, u2)
        return [u1, u2, u3] # Basis matrix (rows)

    B_p = get_basis(p_coords)
    B_q = get_basis(q_coords)
    
    # Rotation R takes B_q to B_p
    # R * B_q^T = B_p^T  => R = B_p^T * B_q
    
    # Transpose bases
    BpT = [[B_p[j][i] for j in range(3)] for i in range(3)] 
    
    # R = BpT * Bq
    R = [[0,0,0], [0,0,0], [0,0,0]]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                R[i][j] += BpT[i][k] * B_q[k][j]
                
    return R, qc, pc

# --- MAIN ---
protein_pdb = r"d:\Ozempic\GROMACS_Setup\semaglutide_solution_input.pdb"
ligand_pdb = r"d:\Ozempic\GROMACS_Setup\New_Simulation\charmm-gui-7152362725\ligandrm.pdb"
output_pdb = r"d:\Ozempic\GROMACS_Setup\semaglutide_merged.pdb"
rtf_path = r"d:\Ozempic\GROMACS_Setup\New_Simulation\sem_fixed.rtf"

# 1. Load Atoms
p_atoms = parse_pdb(protein_pdb)
l_atoms = parse_pdb(ligand_pdb)

# 2. Identify Alignment Anchors
# Protein Res 26: N, CA, C
# Ligand: N1, CA1, C2 (Based on RTF connectivity, C2 matches protein C)
# Let's verify Ligand Anchor Names
p_anchors = {}
p_res26_atoms = [a for a in p_atoms if a['resSeq'] == 26]
for a in p_res26_atoms:
    if a['name'] in ['N', 'CA', 'C']:
        p_anchors[a['name']] = [a['x'], a['y'], a['z']]

l_anchors = {}
for a in l_atoms:
    if a['name'] == 'N1': l_anchors['N'] = [a['x'], a['y'], a['z']]
    if a['name'] == 'CA1': l_anchors['CA'] = [a['x'], a['y'], a['z']]
    if a['name'] == 'C2': l_anchors['C'] = [a['x'], a['y'], a['z']]

if len(p_anchors) != 3 or len(l_anchors) != 3:
    print(f"Error: Missing anchors. P: {p_anchors.keys()}, L: {l_anchors.keys()}")
    exit()

# 3. Calculate Rotation/Translation
p_coords = [p_anchors['N'], p_anchors['CA'], p_anchors['C']]
l_coords = [l_anchors['N'], l_anchors['CA'], l_anchors['C']]

R, centroid_L, centroid_P = get_rotation_matrix(p_coords, l_coords)

# 4. Transform Ligand Atoms
l_transformed = []
for a in l_atoms:
    # Shift to origin (minus centroid L)
    v = [a['x'] - centroid_L[0], a['y'] - centroid_L[1], a['z'] - centroid_L[2]]
    # Rotate
    v_rot = mat_vec_mult(R, v)
    # Shift to target (plus centroid P)
    v_final = [v_rot[0] + centroid_P[0], v_rot[1] + centroid_P[1], v_rot[2] + centroid_P[2]]
    
    a_new = a.copy()
    a_new['x'], a_new['y'], a_new['z'] = v_final
    l_transformed.append(a_new)

# 5. Merge PDBs
# Read RTF to get map of Ligand Atom Name -> Fixed Atom Name
# Because ligandrm.pdb uses CGenFF names (N1, CA1...), protein uses standard (N, CA...)
# We want output to use STANDARD names (matching sem_fixed.rtf)
atom_map_rev = {
    "N1": "N", "CA1": "CA", "C2": "C", "O11": "O", "CB1": "CB",
    "CG1": "CG", "CD1": "CD", "CE": "CE", "NE2": "NZ"
}
# Add H maps from fixed RTF logic if needed, but H usually rebuilt.
# Key is heavy atoms.

final_lines = []
atom_serial = 0

# Write Protein Atoms (Except Res 26 Sidechain, keeping Res 26 Backbone from Protein or Ligand? 
# Use Protein Backbone to be safe with bonds to Res 25/27.
# Use Ligand Sidechain.

for a in p_atoms:
    if a['resSeq'] == 26:
        # Only keep Backbone + O (if standard) or we replace ALL of Res 26 with Transformed Ligand?
        # Replacing ALL is safer for internal geometry consistency, but we must ensure N and C connect to prev/next.
        # Since we aligned N/CA/C, they should match closely.
        continue
    
    atom_serial += 1
    # Format line
    # ATOM  xxxxx  xxxx RRR C xxxx    xxxxxxx xxxxxxx xxxxxxx
    line = f"ATOM  {atom_serial:5d} {a['name']:^4s} {a['resName']:3s} {a['chain']}{a['resSeq']:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00           {a['element']:>2s}\n"
    final_lines.append(line)
    
    # If we just finished Res 25, insert Res 26 here
    if a['resSeq'] == 25 and a['name'] == 'C': # Roughly after Res 25
         pass # Actually we loop purely by order. Res 26 comes after.

# Insert Modified Res 26
# We need to preserve order.
# Let's iterate again properly.

final_lines = []
atom_serial = 0

# Split protein into Pre-26 and Post-26
pre_26 = [a for a in p_atoms if a['resSeq'] < 26]
post_26 = [a for a in p_atoms if a['resSeq'] > 26]

# Write Pre-26
for a in pre_26:
    atom_serial += 1
    line = f"ATOM  {atom_serial:5d} {a['name']:^4s} {a['resName']:3s} {a['chain']}{a['resSeq']:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00           {a['element']:>2s}\n"
    final_lines.append(line)

# Write Res 26 (From Transformed Ligand)
# Filter Ligand atoms: Ligand usually has H. We can keep them or drop them. CHARMM can rebuild H.
# But keeping them is good.
# Also need to Rename atoms according to map.

for a in l_transformed:
    name = a['name']
    mapped_name = atom_map_rev.get(name, name) # specific overrides
    
    # Update res info
    resName = "SEM"
    chain = "A"
    resSeq = 26
    
    atom_serial += 1
    # Note: Column spacing for 4-char atom names like 'HG11'
    if len(mapped_name) == 4:
        name_str = f"{mapped_name:4s}" # Shift left 
    else:
        name_str = f" {mapped_name:^3s}" # Center/Rightish
        
    # Standard PDB format: Name at 12-15. 
    # If 4 chars, 12-15. If <4, usually 13-?.
    
    if len(mapped_name) == 4:
        name_field = f"{mapped_name:<4s}" # e.g. "HG11"
    else:
        name_field = f" {mapped_name:<3s}" # e.g. " N  "
        
    # Careful formatting
    line = f"ATOM  {atom_serial:5d} {name_field} {resName:3s} {chain}{resSeq:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00           {a['element']:>2s}\n"
    final_lines.append(line)

# Write Post-26
for a in post_26:
    atom_serial += 1
    line = f"ATOM  {atom_serial:5d} {a['name']:^4s} {a['resName']:3s} {a['chain']}{a['resSeq']:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00           {a['element']:>2s}\n"
    final_lines.append(line)
    
last_ter = f"TER   {atom_serial+1:5d}      {post_26[-1]['resName']:3s} {post_26[-1]['chain']}{post_26[-1]['resSeq']:4d}\n"
final_lines.append(last_ter)
final_lines.append("END\n")

with open(output_pdb, 'w') as f:
    f.writelines(final_lines)

print(f"Done. Merged PDB saved to {output_pdb}")
