import numpy as np
import os
from ase import Atoms
from ase.io import read, write
import matplotlib.pyplot as plt
from aimdprobe.init_data import init_data, get_raw_traj
from aimdprobe.other_probes.functions import diffusion_rate

"""
plot the movement trajectories of studied atoms in simulations, e.g., Li+ ion diffusion in water
"""

fp = os.getcwd()
fn = 'vasprun.xml'

# get raw data
raw_data = init_data(fp, fn)
raw_traj = get_raw_traj(raw_data)

"""
choose atoms you are interested in, e.g., O and H atoms in H2O, 
with atomic numbers as 80 and 120
"""
atoms_list = [80,120]

time_step = 1 #fs
diffusion_rates, diffusion_rates_avg, steps = diffusion_rate(raw_data, raw_traj, atoms_list, time_step)


fig, axes = plt.subplots(2,1,figsize = (6,8))

for i, atom in enumerate(atoms_list):
    time = np.arange(1,steps)
    axes[i].plot(time, diffusion_rates[i], marker = 'o') # diffusion rates
    axes[i].plot(time, diffusion_rates_avg[i], color='k') # time-averaged rates
    axes[i].set_xlabel('Time (fs)')
    axes[i].set_ylabel('Rate (r$\AA$/fs)')

fig.tight_layout()
fig.savefig('diffusion_rate.png',dpi=100)