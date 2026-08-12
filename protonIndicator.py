import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.transformations import NoJump, unwrap, wrap, center_in_box
import pickle
from MDAnalysis.analysis import align

###########################
base_dir = "/p/scratch/vglut2free/QMMM/03_gating/fromEq/pH5-7prot/HSD133/setup4/f2525_1Cl/hydronium7"
traj_files = [f"{base_dir}/02_prod/{i}ps/mimic.trr" for i in [1,2]]

u = mda.Universe(f"{base_dir}/build/hydronium.gro", traj_files)


ref = mda.Universe(f"{base_dir}/build/hydronium.gro")  # coordinates of first frame
align.AlignTraj(u, ref, select="protein and backbone", in_memory=True).run()
u.trajectory.add_transformations(NoJump())

e191 = u.select_atoms("protein and resid 133 and name NE2")
oxygens = u.select_atoms("resname WAT HO3 HSD HSE and name O NE2 ND1 and around 12 (protein and resid 133 368 256)")  # or "name O*" for water oxygens
hydrogens = u.select_atoms("resname WAT HO3 and name H1 H2 H3 and around 12 (protein and resid 133 368 256)")

############################
# Parameters
cutoff = 1.30  # Å, O–H covalent-like distance
excess_coords = []
rel_vectors = []
rel_distances = []

prev_pos = None  # to track continuity

for ts in u.trajectory:
    
    e191_center = e191.center_of_mass()

    # ---- 1. Check if E191 is protonated ----
    dist_E_H = mda.lib.distances.distance_array(
        e191.positions, hydrogens.positions
    )

    min_dist = dist_E_H.min()

    if min_dist < cutoff:
        # Proton is on E191
        proton_idx = np.unravel_index(np.argmin(dist_E_H), dist_E_H.shape)[1]
        proton_pos = hydrogens.positions[proton_idx]

        excess_coords.append(proton_pos.copy())
        rel_vec = proton_pos - e191_center

    else:
        # ---- 2. Proton is still on water ----
        dists = mda.lib.distances.distance_array(
            oxygens.positions, hydrogens.positions
        )

        coord_numbers = np.sum(dists < cutoff, axis=1)
        hydronium_idx = np.argmax(coord_numbers)

        proton_pos = oxygens.positions[hydronium_idx]

        excess_coords.append(proton_pos.copy())
        rel_vec = proton_pos - e191_center

    rel_vectors.append(rel_vec.copy())
    rel_distances.append(np.linalg.norm(rel_vec))

excess_coords = np.array(excess_coords)
rel_vectors = np.array(rel_vectors)
rel_distances = np.array(rel_distances)
times = np.arange(len(rel_vectors))*0.24*0.1

data = {"times": times, "coord": excess_coords, "coord_rel": rel_vectors, "dis_rel": rel_distances}
with open("coord_HO3.pkl", "wb") as f:
    pickle.dump(data, f)

##################
fig, ax = plt.subplots(3,1, figsize=(3, 3), sharex=True)
labels = ["x", "y", "z"]
for i, label in enumerate(labels):
    ax[i].plot(times, excess_coords[:, i], label=label, color=f"C{i}")
    ax[i].set_ylabel(f"{label} (Å)")
    ax[i].legend()
    
ax[-1].set_xlabel("Time (ps)")
plt.tight_layout()
plt.show()

#####################
fig, ax = plt.subplots(3,1, figsize=(6, 6), sharex=True, sharey=True, gridspec_kw={"hspace": 0})
labels = ["x", "y", "z"]

for i, label in enumerate(labels):
    ax[i].plot(times, rel_vectors[:, i], label=label, color=f"C{i}")
    #ax[i].set_ylabel(f"${label}_{\mathrm{H+}}-{label}_{$\mathrm{E191}$}$ (Å)")
    ax[i].set_ylabel(f"${label}_{{\\mathrm{{H^+}}}}-{label}_{{\\mathrm{{E191}}}}$ (Å)")
    ax[i].legend()

ax[-1].set_xlabel("Time (ps)")
plt.tight_layout()
plt.savefig("coord_HO3.png", dpi=300)
plt.show()

####################
fig, ax = plt.subplots(figsize=(6,3))
labels = ["x", "y", "z"]

for i, label in enumerate(labels):
    ax.plot(times, rel_vectors[:, i], label=label, color=f"C{i}")
    ax.set_ylabel(f"${label}_{{\\mathrm{{H^+}}}}-{label}_{{\\mathrm{{E191}}}}$ (Å)")
    ax.legend(frameon=False)
    
ax.plot(times, rel_distances, color='black')

ax.set_xlabel("Time (ns)")
ax.set_ylabel(r"$d (\mathrm{H_3O^+}-\mathrm{E191})$ (Å)")
plt.tight_layout()
plt.show()
