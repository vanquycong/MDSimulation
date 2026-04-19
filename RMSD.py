#Sometimes, CAN NOT SOLVE PERIODIC CONDITION USING NoJump() if dimer jumps at the end
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import pickle
from MDAnalysis.transformations import NoJump, unwrap, center_in_box
import glob

def rmsd_for_atomgroups(universe, reference, selection1, selection2=None):
    # Align trajectory to the reference
    alignment = align.AlignTraj(mobile=universe, reference=reference, select=selection1, in_memory=True)
    alignment.run()

    # Calculate RMSD
    rmsd_analysis = rms.RMSD(universe, reference=reference, select=selection1, groupselections=selection2)
    rmsd_analysis.run()
    
    # Create a DataFrame to store RMSD results
    columns = [selection1, *selection2] if selection2 else [selection1]
    rmsd_df = pd.DataFrame(np.round(rmsd_analysis.results.rmsd[:, 2:], 2), columns=columns)
    rmsd_df.index.name = "frame"
    return rmsd_df

dimer_backbone = {} # of backbone

# Get all *.dcd files in the current directory
traj_files = glob.glob("system_*.dcd")

# Sort the trajectory files (alphabetical order by default)
traj_files = sorted(traj_files)
top_file = 'system.psf'

# Load the universe with all trajectory files
# Load the topology and trajectory
u = mda.Universe(top_file, *traj_files)
u.trajectory.add_transformations(NoJump())

protein_selection = "protein and backbone"  # Selection for protein backbone atoms

# Set reference as the first frame of the forward trajectory
ref = mda.Universe(top_file, traj_files[0]) 

# Calculate RMSD for the protein
rmsd_df = rmsd_for_atomgroups(u, ref, protein_selection)
dimer_backbone = rmsd_df[protein_selection].values  # Extract the RMSD values
time = np.arange(len(dimer_backbone)) * 1e-1
