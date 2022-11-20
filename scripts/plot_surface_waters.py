"""
plot surface adsorbed solvents (H2O)
"""

import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from aimdprobe.init_data import init_data, get_raw_traj
from aimdprobe.structure_probe.probe_surface_waters import get_adsorbed_h2o
from aimdprobe.useful_functions import get_cumulative_avg

fp = os.getcwd()
fn = 'vasprun.xml'
# fn = 'OUTCAR'

# get raw data
raw_data = init_data(fp, fn)
raw_traj = get_raw_traj(raw_data)

# parameters
nbins = 100
ads_list = [] # no adsorbate in the system
slab_list = np.arange(64) # metal slab has 64 Au atoms
dist = 3 # Angstrom, a general bond length for water adsorption on transition metal surfaces

N_w_ads = []

for traj in raw_traj:
    n_w_ads = get_adsorbed_h2o(raw_data, traj, ads_list, slab_list, dist)
    N_w_ads.append(n_w_ads)

N_w_ads_avg = get_cumulative_avg(N_w_ads)
time = np.arange(len(N_w_ads_avg))*0.001 # ps

fig, ax = plt.subplots(figsize=(8,6))

ax.plot(np.arange(0,runtime)*0.001, N_w_ads, lw = 1, color = 'grey', alpha = 0.5)
ax.plot(np.arange(0,runtime)*0.001, N_w_ads_avg, lw = 3, color = 'black')

ax.annotate('Average H bonds: '+str(round(N_w_ads_avg[-1],2)), (1, 2.5), fontsize = 15)

ax.set_xlabel('Time (ps)', fontsize = 20)
ax.set_ylabel('N / per adsorbate', fontsize = 20)
ax.tick_params(labelsize=20)
ax.set_xlim(0,10)
ax.set_ylim(0,0.3)
fig.tight_layout()
fig.savefig('n_w_ads.png', dpi = 150)

