#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import  read, write
from aimdprobe.energy_probe.functions import get_energy, get_temperatures
from aimdprobe.useful_functions import get_cumulative_avg
import argparse

if __name__ == '__main__':
   """
   plot energies and temperatures
   """
   # In this case, we have n vasprun files, with name as 'vasprun_1.xml', 'vasprun_2.xml', ..., 'vasprun_n.xml'
   # We apply an argument to define the file name and number.

   # 'python plot_energy_multi_input.py n'

   # configure all considered calculations
   parser = argparse.ArgumentParser(description='Give folder name and repetition time')
   #parser.add_argument('filename', default='vasprun', help='configure file name')
   parser.add_argument('n',type=int, default=0, help='configure number of runs n')
   args = parser.parse_args()

   #foldername = args.filename
   n_repeat = args.n
   print(n_repeat)

   N = range(1,n_repeat+1)

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
