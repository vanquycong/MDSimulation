#Use to cut the supercell back into the original cell. 

import parmed as pmd
import numpy as np

struct = pmd.load_file("methylImizodal_super.top",xyz= "100ns.gro")

box = struct.box[:3]

keep_residues = []

for res in struct.residues:
    keep = False

    for atom in res.atoms:
        x, y, z = atom.xx, atom.xy, atom.xz

        fx = x / box[0]
        fy = y / box[1]
        fz = z / box[2]

        if (1/3 < fx < 2/3) and (1/3 < fy < 2/3) and (1/3 < fz < 2/3):
            keep = True
            break

    if keep:
        keep_residues.append(res)

# 🔴 IMPORTANT FIX: convert atoms → indices
atom_indices = []
for res in keep_residues:
    for atom in res.atoms:
        atom_indices.append(atom.idx)

central = struct[atom_indices]

central.save("central_cell.gro", format="gro", overwrite=True)
central.save("central_cell.top", format="gromacs",  overwrite=True)
