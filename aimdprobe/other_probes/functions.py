"""
functions for other probes, in development
"""

import os
import json
import numpy as np
from ase import Atoms
from ase.io import read, write
from aimdprobe.useful_functions import get_cumulative_avg

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
