## Before you start using the scripts...

<br>

- The main functions of AIMDprobe help analyze results from VASP-based AIMD simulations with input file: OUTCAR, vasprun.xml, etc.

- In structure_probe, we split the atomic system into: **slab (slab_list), adsorbate (ads_list) and solvent (solv_list)**. *_list is a list of atomic numbers of respestive component. e.g., slab_list = numpy.arange(64) in a 64Au-40H2O system, which is based on your own ASE-based structures, e.g., POSCAR.

- If you would like to analyze a series of input files, e.g., you have 10 chronological OUTCARs, then you could create a loop while using the functions.

<br>

Let's take energy_probe as an example:

```
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import  read, write
from aimdprobe.energy_probe.functions import get_energy 
from aimdprobe.useful_functions import get_cumulative_avg

outcars = ['outcar1', 'outcar2', ... , 'outcarN']

Pot_energy = []
for outcar in outcars:
    raw_data = read(outcar,':')
    pot_energy, pot_energy_avg = get_energy(raw_data)
    Pot_energy.extend(pot_energy)

Pot_energy_avg = get_cumulative_avg(Pot_energy)
Time = np.arange(len(Pot_energy))

fig, ax = plt.subplots()
ax.plot(Time, Pot_energy)
ax.plot(Time, Pot_energy_avg)

```

An example script to get energies from multiple vasprun.xml files is given as 'plot_energy_multi_input.py'