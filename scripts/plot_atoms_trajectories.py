import numpy as np
import os
from ase import Atoms
from ase.io import read, write
import matplotlib.pyplot as plt
from aimdprobe.init_data import init_data, get_raw_traj
from aimdprobe.other_probes.functions import track_atoms

"""
plot the movement trajectories of studied atoms in simulations, e.g., Li+ ion diffusion in water
"""

fp = os.getcwd()
fn = 'vasprun.xml'

# get raw data
raw_data = init_data(fp, fn)
raw_traj = get_raw_traj(raw_data)

"""
choose atoms you are interested in, e.g., here are Au in the slab, and O and H atoms in H2O, 
with atomic numbers as 1, 80 and 120
"""
atoms_list = [1, 80,120]

# get trajectories for atoms in atoms_list
trajectories = track_atoms(raw_traj, atoms_list)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for atom in atoms_list:
    xs, ys, zs = np.array(trajectories[atom],dtype=object)
    ax.plot3D(xs,ys,zs, marker = 'o')

ax.set_xlabel('X (r$\AA$)')
ax.set_ylabel('Y (r$\AA$)')
ax.set_zlabel('Z (r$\AA$)')

fig.tight_layout()
fig.savefig('track_atoms.png',dpi=100)