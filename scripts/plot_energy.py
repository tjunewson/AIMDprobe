#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import  read, write
from aimdprobe.energy_probe.functions import get_energy, get_temperatures


if __name__ == '__main__':
    """
    plot energies and temperatures
    """

    fn = 'OUTCAR' ## change to your raw data file
    raw_data = read(fn,':')

    pot_energy, pot_energy_avg = get_energy(raw_data)
    temperatures, temperatures_avg = get_temperatures(fn)
    time = np.arange(len(pot_energy)) # fs *0.001 -> ps

    # save the output!
    titles = ['runtime','temperature','pot_energy','pot_energy_avg']
    summary = np.array([time, temperatures, pot_energy, pot_energy_avg]).T
    str_arr = summary.astype(str)

    table = '\t'.join(titles)+'\n'
    for smr in summary:
        table += '\t'.join(smr)+'\n'
    f = open('pot_energy_table.txt','w')
    f.write(table)
    f.close()

    ### check convergences
    conv = if_converge(pot_energy_avg)
    print(conv)
    ### for plot
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

