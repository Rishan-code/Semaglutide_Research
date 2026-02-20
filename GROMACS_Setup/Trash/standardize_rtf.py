
import re

input_rtf = r"d:\Ozempic\GROMACS_Setup\New_Simulation\sem_fixed.rtf"
output_rtf = r"d:\Ozempic\GROMACS_Setup\New_Simulation\sem_standardized.rtf"

# Standard Backbone Charges (from CHARMM36 LYS)
# N: -0.47, HN: 0.31, CA: 0.07, HA: 0.09
# C: 0.51, O: -0.51
# Total Backbone Charge: -0.47 + 0.31 + 0.07 + 0.09 + 0.51 - 0.51 = 0.00
# Perfect.

# Atoms to KEEP (and remap types/charges)
backbone_map = {
    "N":  {"type": "NH1", "charge": -0.47},
    "H":  {"type": "H",   "charge":  0.31}, # Was H1/H2
    "CA": {"type": "CT1", "charge":  0.07},
    "HA1": {"type": "HB1", "charge": 0.09, "rename": "HA"}, # Standard naming
    "C":  {"type": "C",   "charge":  0.51},
    "O":  {"type": "O",   "charge": -0.51}
}

# Atoms to DELETE (Terminal atoms from Ligand definition)
# H2 (Standard N has only H), OXT1, OXT2, HXT1, HXT2
delete_atoms = ["H2", "OXT1", "OXT2", "HXT1", "HXT2", "O1", "O11"] 
# Note: Ligand might have C1=O1, C2=O11, C1-OXT1 etc.
# Need to check structure.
# In sem_fixed.rtf:
# ATOM C1 CG2O2 ... BOND C1 O1 ... BOND C1 OXT1 ..
# ATOM C2 CG2O2 ... BOND C2 O11 ... BOND C2 OXT2 ..
# Wait, C1 and C2 in Ligand Reader are likely the carboxylic tails if it's a zwitterion or di-acid?
# Let's check view_file of RTF again.
# ATOM CA2 connects to C1. ATOM CA1 connects to C2.
# Backbone path?
# In sem_fixed.rtf, we mapped names.
# N -> N1 -> CA -> CA1 -> C -> C2 (Wait, check RTF mapping)
# BOND CA1 N1
# BOND CA1 C2
# So C2 is the Backbone C.
# BOND CA2 C1
# CA2 is sidechain? No.
# Let's re-read the mapping in sem_fixed.rtf from previous turn.
# Line 209: BOND CA1 C2. 
# So C2 was the backbone C. We remapped C2 -> C.
# And N1 -> N.
# So "C" is our backbone carbonyl.
# "C" is bonded to O11 (remapped to O) and OXT2.
# So we need to KEEP "O" (was O11).
# We need to DELETE "OXT2" and "HXT2".
# What about C1/O1/OXT1?
# BOND CA2 C1. CA2 is sidechain?
# BOND CB2 CA2.
# N1-CA1-CB1-CG1...
# Sidechain branches?
# Let's check the structure image or PDB.
# Ligand: N1-CA1-tail...
# But also N1-CA1-C2 (acid).
# And CA1-CB1...
# Where is the other end?
# Ah, Peptide backbone is N-CA-C.
# Semaglutide is a modified LYSINE.
# So it has N, CA, C (backbone).
# And CB, CG, CD, CE, NZ (Sidechain).
# And NZ is attached to the Linker (Glu-gamma...).
# In RTF:
# ATOM N1 (N)
# ATOM CA1 (CA)
# ATOM CB1 (CB) ...
# ATOM NE2 (NZ) -> Matches Lysine NZ.
# BOND CA1 CB1.
# BOND CB1 CG1.
# ...
# BOND CE NE2 (NZ).
# BOND NE2 CD2. -> This is the start of the linker!
# So everything "downstream" from NZ is the sidechain modification.
# WHAT IS CA2/C1/N2 etc?
# Line 23-26: CB2, CA2, N2, C10.
# BOND NZ CD2. BOND CD2 OE1/CG2. BOND CG2 CB2. BOND CB2 CA2.
# It seems the linker goes NZ -> CD2 -> ... -> CA2.
# And CA2 is an alpha carbon of the *linker* (Gamma-Glu)?
# Gamma-Glu backbone: N-CA-C(alpha_COOH)-CB-CG(linked to Lys)-...
# Yes.
# So C1/O1/OXT1 are the Alpha-Carboxyl of the Gamma-GLU linker.
# THIS SHOULD BE KEPT (protonated or deprotonated).
# Since it's free in solution (pH 7), it should be COO-.
# So Clean Up:
# Backbone: N, CA, C, O (Keep, Convert to Protein Types).
# Termini of Backbone: Remove H2 (N-term), Remove OXT2 (C-term). "C" connects to next residue.
# Sidechain Termini: The gamma-Glu carboxyl (C1, O1, OXT1). Keep as COO- or COOH.
# Ligand Reader likely made it COOH or COO-.
# C1 is CG2O2. O1 is OG2D1 (=O). OXT1 is OG311 (-OH or -O-).
# If HXT1 exists, it's COOH.
# We should probably keep it as is from CGenFF sidechain params.

lines = []
with open(input_rtf, 'r') as f:
    lines = f.readlines()

new_lines = []
group_open = False
atom_section = True 

# Standardize Charge Logic
# We need to calculate how much charge we removed/added to the backbone.
# Original N: -1.016, H1:0.368, H2:0.368 -> N group ~ -0.28
# New N: -0.47, H: 0.31 -> -0.16. 
# Original CA: 0.224, HA: 0.09.
# New CA: 0.07, HA: 0.09.
# Original C: 0.744, O: -0.544.
# New C: 0.51, O: -0.51.
# We need to ensure the TOTAL residual charge remains integer (usually -1 or -2 for Semaglutide sidechain).
# Strategy: Update backbone charges to standard. sum the difference. distribute difference onto CB/CG.

original_charge_sum = 0.0
new_charge_sum = 0.0
modified_atoms = []

# First pass: Read original charges for backbone atoms being modified
# Valid backbone atoms: N, CA, C, O, HA1(HA), H1(H).
# Atoms being removed: H2, OXT2, HXT2. (Backbone termini)
# Standard Residue does NOT have H2, OXT2.
# So we compare:
# CHARMM LYS Backbone (Net 0) vs CGenFF Ligand Backbone (Net ?).
# We effectively "cut out" the CGenFF backbone and replace with CHARMM backbone.
# The interface is CA-CB bond.
# We need to preserve electrostatics at CB.
# Best bet: Keep CB charge as is. Just normalize the total residue charge to nearest integer.

# Regex to parse ATOM lines
# ATOM <NAME> <TYPE> <CHARGE>
atom_re = re.compile(r"ATOM\s+(\w+)\s+(\w+)\s+([-\d\.]+)")

final_atom_lines = []
bonds = []
impropers = []

skip_atoms = ["H2", "OXT2", "HXT2"] # Backbone termini removal
# C1/O1/OXT1/HXT1 are sidechain (gamma-Glu). Keep them.

for line in lines:
    line = line.strip()
    if line.startswith("RESI"):
        # RESI SEM 0.000 or similar.
        # We'll set net charge later or keep 0.000 placeholder?
        # CHARMM-GUI often recalculates.
        new_lines.append("RESI SEM           -1.000\n") # Assuming net -1? (Glu-gamma coo-, Lys+... wait amide link. Lys is neutral amide. Glu is -1. AEEA neutral. End Acid?). PDB pka check usually handles this. keeping 0.0 or let grompp complain. 
        # Actually safest is to sum charges.
        continue
        
    if line.startswith("GROUP"):
        new_lines.append(line + "\n")
        continue
        
    m = atom_re.match(line)
    if m:
        name = m.group(1)
        atype = m.group(2)
        charge = float(m.group(3))
        
        if name in skip_atoms:
            continue
            
        if name in backbone_map:
            # Updating Backbone Atom
            info = backbone_map[name]
            new_name = info.get("rename", name)
            new_type = info["type"]
            new_charge = info["charge"]
            
            # Formatted ATOM line.
            # ATOM NAME  TYPE   CHARGE
            # ATOM N     NH1    -0.470
            new_lines.append(f"ATOM {new_name:<4s} {new_type:<6s} {new_charge:6.3f} ! Standard Backbone\n")
        else:
            # Sidechain atom. Keep as is.
            # Wait, if we removed atoms, checking bond section later.
            new_lines.append(line + "\n")
            
    elif line.startswith("BOND") or line.startswith("DOUBLE"):
        # Filter bonds involving deleted atoms
        parts = line.split()
        # BOND A B
        # Could be multiple pairs? "BOND A B C D" - CHARMM allows multiple.
        # Our RTF likely one pair per line or check.
        # CGenFF usually: BOND A B
        
        # Simple parser for "BOND A B" pairs
        # Tokenize whole line
        tokens = line.split()
        keyword = tokens[0] 
        pairs = tokens[1:]
        
        valid_pairs = []
        for i in range(0, len(pairs), 2):
            if i+1 >= len(pairs): break
            a1 = pairs[i]
            a2 = pairs[i+1]
            if a1 in skip_atoms or a2 in skip_atoms:
                continue
            # Handle Renames
            if a1 == "HA1": a1 = "HA"
            if a1 == "H1": a1 = "H"
            if a2 == "HA1": a2 = "HA"
            if a2 == "H1": a2 = "H"
            
            valid_pairs.append((a1, a2))
            
        if valid_pairs:
            s = f"{keyword} "
            for p in valid_pairs:
                s += f"{p[0]:<4s} {p[1]:<4s} "
            new_lines.append(s.strip() + "\n")
            
    elif line.startswith("IMPR"):
        # Filter impropers
        tokens = line.split()
        keyword = tokens[0]
        atoms = tokens[1:] # 4 atoms
        
        if any(a in skip_atoms for a in atoms):
            continue
            
        # Renames
        atoms = ["HA" if a=="HA1" else "H" if a=="H1" else a for a in atoms]
        
        new_lines.append(f"IMPR {' '.join(atoms)}\n")
        
    elif line.startswith("DELETE"):
        continue # Don't want DELETE in RTF usually
        
    else:
        new_lines.append(line + "\n")

# Write output
with open(output_rtf, 'w') as f:
    f.writelines(new_lines)

print(f"Standardized RTF written to {output_rtf}")
