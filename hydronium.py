##########https://parmed.github.io/ParmEd/html/api/parmed/parmed.html?highlight=angle#parmed.Angle

import parmed as pmd
from parmed.tools import change
from parmed import Atom, Bond, Angle
import numpy as np
from parmed.topologyobjects import BondType, AngleType

def create_hydronium(structure, wat_id=116568):
    # Constants
    bond_length = 0.98  # typical O–H bond length in Å
    angle_deg = 113.4   # ideal H–O–H angle in degrees
    angle_rad = np.deg2rad(angle_deg)

    wat_res = structure.atoms[wat_id].residue
    print(wat_res.name, wat_res.number)

    # Select water molecule to convert
    res = structure.residues[wat_res.number]  # zero-based index
    print("Original atoms:", res.atoms)
    
    # Rename residue
    res.name = "HO3"
    
    # Rename atoms
    res.atoms[0].name = "O"
    res.atoms[1].name = "H1"
    res.atoms[2].name = "H2"
    
    # Get coordinates
    O = np.array([res.atoms[0].xx, res.atoms[0].xy, res.atoms[0].xz])
    H1 = np.array([res.atoms[1].xx, res.atoms[1].xy, res.atoms[1].xz])
    H2 = np.array([res.atoms[2].xx, res.atoms[2].xy, res.atoms[2].xz])
    print(O,H1,H2)
    
    # Compute local frame
    v1 = H1 - O
    v2 = H2 - O
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    
    # Bisector direction (between H1 and H2)
    bisector = (v1 + v2)
    bisector /= np.linalg.norm(bisector)
    
    # Normal to the H1–O–H2 plane
    normal = np.cross(v1, v2)
    normal /= np.linalg.norm(normal)
    
    # Rotate bisector by angle to move out of plane
    # H3 lies at angle θ from the bisector in the plane perpendicular to it
    cos_theta = np.cos(angle_rad)
    sin_theta = np.sin(angle_rad)
    
    # Use Rodrigues' rotation formula to get final direction
    H3_dir = cos_theta * bisector + sin_theta * normal
    print(H3_dir)
    H3_coord = O + bond_length * H3_dir
    
    # Create the third H atom
    H3 = Atom(name="H3", type="HT", atomic_number=1, charge=1.0, mass=1) #,  Charge = 1.0 instead of 0.417, do not use :type="HT",
    H3.xx, H3.xy, H3.xz = H3_coord
    
    print(H3_coord)
    
    print(H3.number)  
    print(res.number)  # Should be True
    
    #structure.strip(f':{wat_res.number} & :WAT')
    
    # Add H3 to the same residue
    #structure.atoms.append(H3)
    #res.atoms.append(H3)
    #structure.res.append(H3)
    structure.add_atom(H3, res.name, res)
    
    H3.residue = res
    #res.add_atom(H3)
    
    # Add bonds
    bond_length = 0.09572   # in nm (0.9572 Å)
    bond_k = 462750.0         # in kJ/mol·nm² (Amber unit conversion)
    
    #structure.bonds.append(Bond(res.atoms[0], res.atoms[1]))  # O–H1
    #structure.bonds.append(Bond(res.atoms[0], res.atoms[2]))  # O–H2
    structure.bonds.append(Bond(res.atoms[0], H3, type=BondType(bond_k, bond_length)))# O–H3 (new)
    
    #Add angle
    angle_value = 113.40     # degrees for H–O–H in H3O+
    angle_k = 460.24        # kJ/mol·rad² (or similar, depending on FF)
    structure.angles.append(Angle(res.atoms[1], res.atoms[0], H3, type=AngleType(angle_k, angle_value)))
    structure.angles.append(Angle(res.atoms[2], res.atoms[0], H3, type=AngleType(angle_k, angle_value)))

structure = pmd.load_file("fix_wat.top", xyz="mini.gro")
create_hydronium(structure, wat_id=116568) #Need index of the water oxygen as input e.g. see in VMD. 
create_hydronium(structure, wat_id=117065)

#structure.save("hydronium.gro", format="gro", overwrite=True)
#structure.save("hydronium.pdb", overwrite=True)
structure.save("hydronium.psf", format="psf", overwrite=True)
#structure.save("hydronium.top", format="gromacs", overwrite=True)
