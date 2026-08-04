import parmed as pmd
import csv

def apply_protonation_state(amber, protonation_dict):
    # Define atom names by residue and itauto
    atom_map = {
        "GLU": {3: "HE1", 4: "HE2"},
	"ASP": {3: "HD1", 4: "HD2"},
        "HSP": {1: "HD1", 2: "HE2"},
    }

    for res in amber.residues:
        resid = res.number + 1 # Because Amber index starts at 1
        if resid not in protonation_dict:
            continue

        state = protonation_dict[resid]['State']
        site = protonation_dict[resid]['Protonation Site']
        resname = res.name
        print(resname, resid)

        if state == 'deprotonated' and resname in ["ASP", "GLU"]:
            # Remove both protons if deprotonated
            for name in atom_map.get(resname, {}).values():
                amber.strip(f':{resid}@{name}')
        elif state == 'protonated' and resname in ["ASP", "GLU"]:
            # Remove *other* proton
            for it, name in atom_map.get(resname, {}).items():
                if it != site:
                    amber.strip(f':{resid}@{name}') 
                #amber.strip(f':{resid}@HE1') if resname == "GLU" else amber.strip(f':{resid}@HD1') #Because Gromacs gmx2pdb does not distinguished HD2 and HD1 in CHARMM36m
        elif state == 'deprotonated' and resname in ["HSP"]:
            # Remove *other* proton
            for it, name in atom_map.get(resname, {}).items():
                if it != site:
                    amber.strip(f':{resid}@{name}')
        elif state == 'intermediate':
            print(f"Skipping intermediate residue {resname}{resid}")

def load_protonation_csv(csv_filename):
    protonation_dict = {}
    with open(csv_filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            resid = int(row["Residue"])
            state = row["State"]
            site = row["Protonation Site"]
            print(site)
            protonation_dict[resid] = {
                "State": state,
                "Protonation Site": int(site) if site.isdigit() else None
            }
    return protonation_dict

def change_prot_name(amber, protonation_dict): #CHARMM36m
    for res in amber.residues:
        resid = res.number + 1
        if resid not in protonation_dict:
            continue

        state = protonation_dict[resid]["State"]
        site = protonation_dict[resid]["Protonation Site"]
        print(resid,res.name,state)	

        if res.name in ["HSP"]: #Amber, HIP
            if state == "deprotonated" and site == 1:
                res.name = "HSD"  # proton at ND1, Amber: HID
            elif state == "deprotonated" and site == 2:
                res.name = "HSE"  # proton at NE2, Amber: HIE
        elif res.name in ["GLU"]:
            if state in ["protonated"]:
                res.name = "GLUP"  # protonated GLU, Amber: GLH
            else:
                res.name = "GLU"
        elif res.name in ["ASP"]:
            if state in ["protonated"]:
                res.name = "ASPP"  # protonated ASH
            else:
                res.name = "ASP"
                
def fix_charges(amber, protonation_dict, tolerance=1e-3):
    for res in amber.residues:
        resid = res.number + 1
        if resid not in protonation_dict:
            continue  # Only check titratable residues

        rname = res.name
        if rname not in CHARMM36M_CHARGES:
            print(f" ⚠️ No CHARMM charges for residue {resid} ({rname}) — skipping.")
            continue

        charmm_ref = CHARMM36M_CHARGES[rname]
        total_charge = 0.0

        for atom in res.atoms:
            atom_name = atom.name
            if atom_name not in charmm_ref:
                print(f"    ⚠️ Atom {atom_name} not found in CHARMM reference for {rname}")
                total_charge += atom.charge
                continue

            ref_charge = charmm_ref[atom_name]
            if abs(atom.charge - ref_charge) > tolerance:
                print(f"    ↪️ Fixing charge of {atom_name} in residue {resid} ({rname}): {atom.charge:.3f} → {ref_charge:.3f}")
                atom.charge = ref_charge
            total_charge += atom.charge

        expected = EXPECTED_CHARGES.get(rname)
        print(f"[{resid:4d}] {rname} total charge after fix: {total_charge:.3f}")
        if expected is not None and abs(total_charge - expected) > 0.1:
            print(f"    ⚠️ Charge mismatch: got {total_charge:.3f}, expected {expected:.2f}")  

CHARMM36M_CHARGES = {
    "GLU": {"CD": 0.6200, "OE1": -0.7600, "OE2": -0.7600, "CG": -0.2800}, 
    "GLUP": {"CD": 0.7500, "OE1": -0.5500, "OE2": -0.6100, "HE2": 0.4400, "HE1": 0.4400}, # Expand with HE1
    "ASP": {"CG": 0.6200, "OD1": -0.7600, "OD2": -0.7600, "CG": -0.2800},
    "ASPP": {"CG": 0.7500, "OD1": -0.5500, "OD2": -0.6100, "HD2": 0.4400, "HD1": 0.4400}, # Expand with HD1
    "HSP": {"ND1": -0.5100, "HD1": 0.4400, "NE2": -0.5100, "HE2": 0.4400, "CE1": 0.3200, "HE1": 0.1800, "CD2": 0.1900, "HD2": 0.1300, "CG": 0.1900},
    "HSE": {"ND1": -0.7000, "CG": 0.2200, "CE1": 0.2500, "HE1": 0.1300, "NE2": -0.3600, "HE2": 0.3200, "CD2": -0.0500, "HD2": 0.0900},
    "HSD": {"ND1": -0.3600, "HD1": 0.3200, "CG": -0.0500, "CE1": 0.2500, "HE1": 0.1300, "NE2": -0.7000, "CD2": 0.2200, "HD2": 0.1000}
}
EXPECTED_CHARGES = {
    "GLU": -1.00, "GLUP": 0.00,
    "ASP": -1.00, "ASPP": 0.00,
    "HSD": 0.00, "HSE": 0.00, "HSP": +1.00
}              
# Load the structure
amber = pmd.load_file("cphmd.parm7", "frame4690.pdb")

# Load protonation data from CSV
protonation_dict = load_protonation_csv("protonation_state.csv")

# Apply protonation edits
apply_protonation_state(amber, protonation_dict)
change_prot_name(amber, protonation_dict)
fix_charges(amber, protonation_dict)

# Save modified structure
amber.save("frame4690_fixed.pdb", overwrite=True)
amber.save("frame4690_fixed.top", format="gromacs", overwrite=True)
amber.save("frame4690_fixed.gro", format="gro",overwrite=True)
