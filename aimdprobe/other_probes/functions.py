"""
functions for other probes, in development
"""

import os
import json
import numpy as np
from ase import Atoms
from ase.io import read, write, animation
from aimdprobe.useful_functions import get_cumulative_avg, get_real_distance

"""
track the trajectory of certain atoms during simulations, e.g., Li+ ion diffusion in water
"""

def track_atoms(raw_traj, atoms_list):
    """
    retrieve the 3D trajectory of atoms in atoms_list ('[atomic numbers]')
    as a dictionary with keys as numbers in atoms_list and values as [X,Y,Z] for each atom
    """
    bins = []
    
    for atom in atoms_list:
        ibinx = []
        ibiny = []
        ibinz = []
        for traj in raw_traj:
            ibinx.append(traj[atom][0])
            ibiny.append(traj[atom][1])
            ibinz.append(traj[atom][2])
        bins.append([ibinx, ibiny, ibinz])
    
    trajectories = dict(zip(atoms_list, bins))
    return trajectories

def diffusion_rate(raw_data, raw_traj, atoms_list, time_step):
    """
    calculate the diffusion rate of certain atoms with a time step (e.g., fs)
    """
    steps = len(raw_traj)
    diffusion_rates = []
    diffusion_rates_avg = []
    for atom in atoms_list:
        Rate = []
        for step in np.arange(1, steps):
            coordinate = raw_traj[step][atom]
            coordinate_pre = raw_traj[step-1][atom]
            distance = get_real_distance(raw_data, coordinate, coordinate_pre)
            rate = distance/time_step #the distance an atom covers during 1-step simulation
            Rate.append(rate)
        diffusion_rates.append(Rate)
        diffusion_rates_avg.append(get_cumulative_avg(Rate))

    return diffusion_rates, diffusion_rates_avg, steps

def animate(filename, raw_data, itv, size, timeframe, rotation):
    """
    make an animation based on AIMD trajectories;
    images = raw_data
    itv = slices of images
    size = (a,b,c), enlargement of the unit cell by (a, b, c)
    """
    data = [raw_data[::itv][i]*size for i in np.arange(len(raw_data[::itv]))]
    animation.write_gif(filename+'.gif', data, interval=timeframe, rotation=rotation)