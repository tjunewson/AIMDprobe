###plot radial distribution of solvents
import os
import csv
import json
import numpy as np
from ase import Atoms
from ase.io import  read, write
import matplotlib.pyplot as plt 

ads_list = []
slab_list = np.arange(64)

fp = os.getcwd()
fn = 'au_clean1_vasprun.xml'
#fn = 'OUTCAR'
raw_data = init_data(fp, fn)
raw_traj = get_raw_traj(raw_data)

nbins = 100
O_bins = np.zeros(nbins)
H_bins = np.zeros(nbins)
for traj in raw_traj:
    walls, o_bins, h_bins = get_rdf_3d(raw_data, traj, ads_list, slab_list, nbins)
    O_bins += o_bins
    H_bins += h_bins

O_rdf = O_bins/len(raw_traj)
H_rdf = H_bins/len(raw_traj)
walls = walls

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(walls, O_rdf, lw = 2, color = 'red', label = 'O')
ax.plot(walls, H_rdf, lw = 2, color = 'blue', label = 'H')

ax.set_xlabel('Z (A)', fontsize = 20)
ax.set_ylabel('g(z)', fontsize = 20)
ax.tick_params(labelsize=20)
ax.legend(fontsize = 20)
fig.tight_layout()
fig.savefig('rdf.png', dpi = 150)