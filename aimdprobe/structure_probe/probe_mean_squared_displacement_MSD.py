import numpy as np
from ase import Atoms
from ase.io import read, write
from aimdprobe.useful_functions import get_real_distance
from aimdprobe.structure_probe.functions import get_positions, get_elements_atoms
from aimdprobe.structure_probe.probe_mass_centers import mass_centers

###### main functions ######
"""
probe the mean sqaured displacement of an atom or atom ensembles using mass centers.

both MSD at xyz and xy dimentions are calculated to indicate desorption tendency
"""

def get_msd(raw_data, atoms_list):

    mcenters = mass_centers(raw_data, atoms_list)
    frames = len(mcenters)
    msd_xyz = []
    msd_xy = []
    for lag in range(1, frames):
        displacement_xyz = []
        displacement_xy = []
        for i in range(frames - lag):
            dx = mcenters[i+lag][0] - mcenters[i][0]
            dy = mcenters[i+lag][1] - mcenters[i][1]
            dz = mcenters[i+lag][2] - mcenters[i][2]
            displacement_xyz.append(dx*dx + dy*dy + dz*dz)
            displacement_xy.append(dx*dx + dy*dy)
        msd_xyz.append(np.mean(displacement_xyz))
        msd_xy.append(np.mean(displacement_xy))
    
    return msd_xyz, msd_xy
        








