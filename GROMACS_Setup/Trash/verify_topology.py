"""Verify Semaglutide (PROH) topology against known sequence."""

itp_file = r"d:\Ozempic\GROMACS_Setup\toppar\PROH.itp"

lines = open(itp_file).readlines()

# Parse atoms section
in_atoms = False
atoms = []
for line in lines:
    if "[ atoms ]" in line:
        in_atoms = True
        continue
    if in_atoms and line.strip().startswith("["):
        break
    if in_atoms and not line.startswith(";") and line.strip():
        parts = line.split()
        if len(parts) >= 8 and parts[0].isdigit():
            atoms.append({
                "nr": int(parts[0]),
                "type": parts[1],
                "resnr": int(parts[2]),
                "resname": parts[3],
                "atom": parts[4],
                "charge": float(parts[6]),
                "mass": float(parts[7]),
            })

# Residue sequence
residues = []
seen = set()
for a in atoms:
    key = a["resnr"]
    if key not in seen:
        seen.add(key)
        residues.append((a["resnr"], a["resname"]))

# Reference Semaglutide sequence (GLP-1 positions 7-36)
ref = [
    "His", "Aib", "Glu", "Gly", "Thr", "Phe", "Thr", "Ser", "Asp", "Val",
    "Ser", "Ser", "Tyr", "Leu", "Glu", "Gly", "Gln", "Ala", "Ala", "Lys*",
    "Glu", "Phe", "Ile", "Ala", "Trp", "Leu", "Val", "Arg", "Gly", "Arg",
]

mapping = {
    "His": "HSD", "Aib": "AIB", "Glu": "GLU", "Gly": "GLY", "Thr": "THR",
    "Phe": "PHE", "Ser": "SER", "Asp": "ASP", "Val": "VAL", "Tyr": "TYR",
    "Leu": "LEU", "Gln": "GLN", "Ala": "ALA", "Lys*": "KWF1",
    "Trp": "TRP", "Arg": "ARG", "Ile": "ILE",
}

print("=" * 60)
print("SEMAGLUTIDE TOPOLOGY VERIFICATION")
print("=" * 60)

# Sequence comparison
print(f"\n{'Pos':>4}  {'GLP1#':>5}  {'Sema':>6}  {'ITP':>6}  {'Status':>8}")
print("-" * 40)
all_match = True
for i, ((resnr, resname), ref_name) in enumerate(zip(residues, ref), 1):
    expected = mapping.get(ref_name, "???")
    ok = expected == resname
    if not ok:
        all_match = False
    status = "OK" if ok else "MISMATCH!"
    print(f"{i:4d}  {resnr:5d}  {ref_name:>6}  {resname:>6}  {status:>8}")

if len(residues) != len(ref):
    print(f"\nWARNING: Length mismatch! ITP has {len(residues)} residues, expected {len(ref)}")
    all_match = False

# Charge analysis
total_charge = sum(a["charge"] for a in atoms)
total_mass = sum(a["mass"] for a in atoms)

# Per-residue charge
print(f"\n{'Res':>4}  {'Name':>5}  {'Charge':>8}  {'Atoms':>5}")
print("-" * 30)
for resnr, resname in residues:
    res_atoms = [a for a in atoms if a["resnr"] == resnr]
    res_charge = sum(a["charge"] for a in res_atoms)
    print(f"{resnr:4d}  {resname:>5}  {res_charge:>8.3f}  {len(res_atoms):5d}")

# C-terminus check
print("\n--- C-terminus (last 5 atoms) ---")
for a in atoms[-5:]:
    print(f"  atom {a['nr']:4d}  {a['type']:>6}  {a['resname']:>4}  {a['atom']:>4}  q={a['charge']:+.3f}")

# N-terminus check
print("\n--- N-terminus (first 5 atoms) ---")
for a in atoms[:5]:
    print(f"  atom {a['nr']:4d}  {a['type']:>6}  {a['resname']:>4}  {a['atom']:>4}  q={a['charge']:+.3f}")

print(f"\n{'='*60}")
print(f"Total atoms:    {len(atoms)}")
print(f"Total residues: {len(residues)}")
print(f"Total charge:   {total_charge:+.4f}")
print(f"Total mass:     {total_mass:.2f} Da")
print(f"Sequence match: {'ALL OK' if all_match else 'ISSUES FOUND'}")

# Check KWF1 (modified Lys26) details
print(f"\n--- KWF1 (Modified Lysine @ pos 26) ---")
kwf_atoms = [a for a in atoms if a["resname"] == "KWF1"]
kwf_charge = sum(a["charge"] for a in kwf_atoms)
kwf_mass = sum(a["mass"] for a in kwf_atoms)
print(f"Atoms: {len(kwf_atoms)}")
print(f"Charge: {kwf_charge:+.4f}")
print(f"Mass: {kwf_mass:.2f} Da")
# Print atom names in KWF1
print("Atom names:", ", ".join(a["atom"] for a in kwf_atoms))
