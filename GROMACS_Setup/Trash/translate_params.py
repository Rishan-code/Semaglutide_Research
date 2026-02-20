
input_prm = r"d:\Ozempic\GROMACS_Setup\New_Simulation\charmm-gui-7152362725\res\res.prm"
output_prm = r"d:\Ozempic\GROMACS_Setup\New_Simulation\sem_standardized.prm"

# Mapping from CGenFF Types (Old) to Protein Types (New)
# Based on sem_standardized.rtf
type_map = {
    "NG321": "NH1",   # Backbone N
    "HGPAM2": "H",    # Backbone H (Amide)
    "CG311": "CT1",   # Backbone CA
    "HGA1": "HB1",    # Backbone HA
    "CG2O2": "C",     # Backbone C (Carbonyl - check if mapped correctly)
    "OG2D1": "O",     # Backbone O (Carbonyl)
    # Note: CGenFF might use OG2D1 for both Backbone O and Sidechain O?
    # In RTF: Backbone O is O8, O2, etc? No.
    # Backbone O is "O" (was OG2D1).
    # Sidechain Carbonyls (O2, O5, O8) are also OG2D1.
    # Sidechain Carboxyls (O1, O9) are OG2D1.
    # Sidechain Hydroxyls (O10, OXT1) are OG311.
    # Implication: If we map OG2D1 -> O globally, we change sidechain carbonyls to type "O" (Peptide Carbonyl).
    # This might be Okay (chemically similar amide/ester carbonyls).
    # BUT standard "O" atom type is very specific to Peptide Backbone (H-bonding etc).
    # A sidechain ester C=O might use "OB" or "OS"?
    # BETTER STRATEGY: 
    # Only map parameters if the OTHER atoms in the bond/angle are BACKBONE atoms.
    # If a parameter is `CG321 - OG2D1` (Sidechain C - Sidechain O), we should keep it as `CG321 - OG2D1`?
    # Or should we replicate it as `CG321 - O` just in case?
    # Actually, we didn't change the atom types of the Sidechain atoms in the RTF.
    # Sidechain atoms KEPT their original types (NG2S1, CG2O1, OG2D1 etc).
    # We ONLY changed Backbone Atoms (N, CA, C, O, H, HA).
    # So we only need to map parameters involving `NG321`, `CG311`, `CG2O2` (Backbone C only?), `HGPAM2`.
    # Wait, is `CG2O2` unique to Backbone C?
    # In RTF:
    # ATOM C40 CG2O2 (Sidechain Carboxyl C) -> We kept this type.
    # ATOM C1 CG2O2 (Linker Carboxyl C) -> We kept this type.
    # ATOM C CG2O2 (Backbone C) -> We CHANGED to "C".
    # So `CG2O2` is used for BOTH Backbone and Sidechain.
    # If we map `CG2O2` -> `C` globally in PRM, we create parameters for `C` (Backbone type) which is correct for Backbone C.
    # But we define `C40` as type `CG2O2`. It needs parameters for `CG2O2`.
    # So `res.prm` MUST keep the original `CG2O2` parameters.
    # And we MUST ADD duplicated parameters where `CG2O2` is replaced by `C`.
    # BUT only for the Backbone C context?
    # Actually, chemically `C` (Peptide) and `CG2O2` (Carboxyl) are similar.
    # If we blindly duplicate `CG2O2` lines to `C` lines, we provide params for Backbone C.
    # Sidechain C40 will still use original `CG2O2` params (which correspond to the remaining lines).
    # So duplication is safe.
    # The only risk is if `C` type behaves differently than `CG2O2` for the Sidechain interactions? no, sidechain uses CG2O2.
}

lines = []
with open(input_prm, 'r') as f:
    lines = f.readlines()

new_lines = []
new_lines.append("! Additional parameters translated from CGenFF for Protein Backbone compatibility\n")

seen_params = set()

for line in lines:
    line = line.strip()
    if line.startswith("!") or not line:
        continue
        
    # Check section (BONDS, ANGLES, DIHEDRALS, IMPROPERS)
    # usually CGenFF prm has "BONDS", "ANGLES" headers.
    # We can just process every line that looks like a parameter.
    # Parameter line format: ATOM1 ATOM2 [ATOM3] [ATOM4] K ... 
    
    parts = line.split()
    if not parts: continue
    
    # Heuristic: Uppercase atom names/types usually.
    # We look for our keys in the first 4 columns.
    
    # Bonds (2 atoms), Angles (3 atoms), Dihedrals/Impropers (4 atoms).
    # Nonbonded (atom eps rmin) - usually at end.
    
    # We iterate and check if any atom type matches our map keys.
    # If so, we create a copy of the line with substitutions.
    
    # Note: Valid lines start with atom types.
    if len(parts) < 2: 
        # Section header "BONDS" etc?
        new_lines.append(line + "\n")
        continue

    # Try to map
    mapped_types = []
    has_mapping = False
    
    # Determine number of atom columns?
    # Bonds: 2. Angles: 3. Dihedrals: 4. Improper: 4.
    # But followed by numbers.
    # We can try to map the first 4 tokens if they are strings.
    # (Actually simpler: just replace text if it matches robustly?)
    # "CG311" -> "CT1".
    # But be careful not to replace partial strings. Use split/join.
    
    current_line_tokens = list(parts)
    modified_tokens = list(parts)
    
    # Identify how many atoms? 
    # We can just check key matches in the tokens.
    
    for i, token in enumerate(current_line_tokens):
        # Stop processing if we hit numbers (params)
        # But atom types can resemble? No, CGenFF types are alphanumeric.
        # Params are floats.
        try:
            float(token)
            # It's a number, stop checking for types here (usually types come first)
            # Exception: "10.0" k value.
            # But we only map the TYPES.
            continue 
        except ValueError:
            pass
            
        if token in type_map:
            modified_tokens[i] = type_map[token]
            has_mapping = True
            
    # Include original line
    new_lines.append(line + "\n")
    
    # Include translated line if mapping occurred
    if has_mapping:
        new_line = " ".join(modified_tokens)
        # Avoid duplication if possible?
        if new_line != line:
            new_lines.append(new_line + " ! Translated\n")

# Write output
with open(output_prm, 'w') as f:
    f.writelines(new_lines)

print(f"Translated PRM written to {output_prm}")
