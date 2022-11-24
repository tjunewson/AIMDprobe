"""
plot energies and temperatures
"""
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import  read, write
from aimdprobe.energy_probe.functions import get_energy, get_temperatures
from aimdprobe.useful_functions import get_cumulative_avg

#in this case, we have 6 vasprun files, with name as 'vasprun_n.xml'
N = range(1,6)

Pot_energy = []
for n in N:
   fn = 'vasprun_'+str(n)+'.xml'
   raw_data = read(fn,':')
   pot_energy, pot_energy_avg = get_energy(raw_data)
   Pot_energy.extend(pot_energy)
print('Total steps: ' + str(len(Pot_energy)))

Pot_energy_avg = get_cumulative_avg(np.array(Pot_energy))
time = np.arange(len(Pot_energy))*0.001

#plotting
fig, ax = plt.subplots(figsize = (8,4))

ax.plot(time,Pot_energy,lw=1)
ax.plot(time,Pot_energy_avg,lw=3)

#ax1.plot(len(temperatures),temperatures,lw=1)
#ax1.plot(len(temperatures),temperatures_avg,lw=3)

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Energy (eV)')

fig.tight_layout()
fig.savefig('energy.png',dpi=100)
