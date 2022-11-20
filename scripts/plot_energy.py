"""
plot energies and temperatures
"""
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import  read, write

from aimdprobe.energy_probe.functions import get_energy, get_temperatures

#get raw_data
fn = 'OUTCAR' ## change to your raw data file
raw_data = read(fn,':')

pot_energy, pot_energy_avg = get_energy(raw_data)
temperatures, temperatures_avg = get_temperatures(fn)
time = np.arange(len(raw_data))*0.001 ## ps

###for plot
fig, (ax, ax1) = plt.subplots(1,2,figsize = (12,4))

ax.plot(time,pot_energy,lw=1)
ax.plot(time,pot_energy_avg,lw=3)

ax1.plot(time,temperatures,lw=1)
ax1.plot(time,temperatures_avg,lw=3)

for ax in [ax, ax1]:
    ax.set_xlabel('Time (ps)')
    if ax == ax1:
        ax.set_ylabel('T (K)')
    else:
        ax.set_ylabel('Energy (eV)')

fig.tight_layout()
fig.savefig('energy.png',dpi=100)