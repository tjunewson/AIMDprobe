#!/usr/bin/env python

import os
import numpy as np
from ase import Atoms
from ase.io import  read, write
import matplotlib.pyplot as plt 
from aimdprobe.init_data import init_data, get_raw_traj
from aimdprobe.structure_probe.probe_water_distribution import get_rdf 

if __name__ == "__main__":
    """
    plot radial distribution of solvents (H2O) using AIMDprobe
    """
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

    O_rdf = np.zeros(nbins)
    H_rdf = np.zeros(nbins)

    for traj in raw_traj:
        walls, o_rdf, h_rdf = get_rdf(raw_data, traj, ads_list, slab_list, nbins)
        O_rdf += o_rdf
        H_rdf += h_rdf

    O_rdf = O_rdf/len(raw_traj)
    H_rdf = H_rdf/len(raw_traj)
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