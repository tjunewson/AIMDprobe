"""
plot numbers of hydrogen bonds in solvents (H2O) using AIMDprobe
"""

import os
import numpy as np
from ase import Atoms
from ase.io import  read, write
import matplotlib.pyplot as plt 
from aimdprobe.init_data import init_data, get_raw_traj
from aimdprobe.structure_probe.probe_h_bonds import get_hbonds_sol

fp = os.getcwd()
fn = 'vasprun.xml'
# fn = 'OUTCAR'

# get raw data
raw_data = init_data(fp, fn)
raw_traj = get_raw_traj(raw_data)

# parameters
ads_list = [] # no adsorbate in the system
slab_list = np.arange(64) # metal slab has 64 Au atoms
hdist = 2.55 # Angstrom, a general hydrogen bond length H-O---H

H_bonds = []

for traj in raw_traj:
    h_bonds = get_hbonds_sol(raw_data, traj, ads_list, slab_list, hdist)
    H_bonds.append(h_bonds)

H_bonds_avg = get_cumulative_avg(h_bonds)
time = np.arange(len(H_bonds_avg))*0.001 # ps

fig, ax = plt.subplots(figsize=(8,6))

ax.plot(np.arange(0,runtime)*0.001, H_bonds_avg, lw = 1, color = 'grey', alpha = 0.5)
ax.plot(np.arange(0,runtime)*0.001, H_bonds_avg, lw = 3, color = 'black')

ax.annotate('Average H bonds: '+str(round(H_bonds_avg[-1],2)), (1, 2.5), fontsize = 15)

ax.set_xlabel('Time (ps)', fontsize = 20)
ax.set_ylabel('N', fontsize = 20)
ax.tick_params(labelsize=20)
ax.set_xlim(0,10)
ax.set_ylim(-0.2,4.2)
fig.tight_layout()
fig.savefig('n_hbonds.png', dpi = 150)

