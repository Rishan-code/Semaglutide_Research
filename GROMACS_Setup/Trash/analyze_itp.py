import re

def parse_itp(file_path, residue_name="KWF1"):
    atoms = {} # id -> name
    bonds = []
    
    in_atoms = False
    in_bonds = False
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
                
            if line.startswith('[ atoms ]'):
                in_atoms = True
                in_bonds = False
                continue
            elif line.startswith('[ bonds ]'):
                in_atoms = False
                in_bonds = True
                continue
            elif line.startswith('['):
                in_atoms = False
                in_bonds = False
                continue
                
            if in_atoms:
                parts = line.split()
                # nr type resnr residu atom cgnr charge mass
                # 263 NH1 26 KWF1 N 263 -0.47 14.007
                if len(parts) >= 5:
                    atom_nr = int(parts[0])
                    res_name = parts[3]
                    atom_name = parts[4]
                    if res_name == residue_name:
                        atoms[atom_nr] = atom_name
            
            if in_bonds:
                parts = line.split()
                if len(parts) >= 2:
                    a1 = int(parts[0])
                    a2 = int(parts[1])
                    if a1 in atoms and a2 in atoms:
                        bonds.append((a1, a2))
    
    # Analyze Graph
    print(f"Residue: {residue_name}")
    print(f"Atom Count: {len(atoms)}")
    print(f"Bond Count: {len(bonds)}")
    
    # Reconstruct connectivity from Backbone CA (Carbon Alpha)
    # Typically Lys CA is the branching point.
    # Find CA
    ca_id = None
    for aid, name in atoms.items():
        if name == 'CA':
            ca_id = aid
            break
            
    print(f"CA Atom ID: {ca_id}")
    
    if ca_id:
        # Build adjacency
        adj = {a: [] for a in atoms}
        for a1, a2 in bonds:
            adj[a1].append(a2)
            adj[a2].append(a1)
            
        # Recursive DFS to find the longest path from CA
        # This will show the main chain length typically
        
        def get_heavy_atom_path(start_node, visited):
            # Only traverse heavy atoms (exclude H)
            max_path = []
            
            is_heavy = lambda nid: not atoms[nid].startswith('H')
            
            neighbors = [n for n in adj[start_node] if n not in visited and is_heavy(n)]
            
            if not neighbors:
                return [atoms[start_node]]
                
            for n in neighbors:
                new_visited = visited.copy()
                new_visited.add(n)
                sub_path = get_heavy_atom_path(n, new_visited)
                if len(sub_path) > len(max_path):
                    max_path = sub_path
                    
            return [atoms[start_node]] + max_path

        # Start from CB (first sidechain atom)
        cb_id = None
        for aid, name in atoms.items():
            if name == 'CB':
                cb_id = aid
                break
                
        if cb_id:
            print("\nLongest Sidechain Path (Heavy Atoms):")
            path = get_heavy_atom_path(cb_id, {ca_id, cb_id}) # Exclude CA to go down sidechain
            print(" -> ".join(path))
            print(f"Path Length: {len(path)}")
            
        # Also print all heavy atoms in sidechain to verify branching
        print("\nAll Heavy Atoms in Sidechain:")
        queue = [cb_id]
        visited_sc = {ca_id, cb_id}
        heavy_atoms = []
        while queue:
            curr = queue.pop(0)
            if not atoms[curr].startswith('H'):
                heavy_atoms.append(f"{atoms[curr]}({curr})")
            
            for n in adj[curr]:
                if n not in visited_sc:
                    visited_sc.add(n)
                    queue.append(n)
        print(", ".join(heavy_atoms))


    print("\nAtoms List:")
    for aid in sorted(atoms.keys()):
        print(f"{aid}: {atoms[aid]}")

parse_itp("d:/Ozempic/GROMACS_Setup/toppar/PROH.itp")
