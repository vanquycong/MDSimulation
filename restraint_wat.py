## The script creates new molecules type for restrained waters
## Input: file with water residue number

import parmed as pmd

# Load topology and coordinates
s = pmd.load_file("hydronium.top", xyz="hydronium.gro")

# ---- Step 1: Read selected residue IDs ----
selected_resids = set()
with open("restraint_wat_resid.txt") as f:
    for line in f:
        for item in line.split():
            selected_resids.add(int(item)-1)

print("Residues to rename:", selected_resids)


# ---- Step 3: Rename selected water residues ----
for res in s.residues:
    if res.number in selected_resids and res.name == "WAT":
        res.name = "rWAT"

# ---- Step 4: Write updated files ----
s.save("hydronium_res.gro", format="gro", overwrite=True)
s.save("hydronium_res.pdb", overwrite=True)
s.save("hydronium_res.top", format="gromacs", overwrite=True)
