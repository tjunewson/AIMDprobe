import numpy as np
import os
from ase import Atoms
from ase.io import read, write
import matplotlib.pyplot as plt
from aimdprobe.init_data import init_data, get_raw_traj


def mass_centers(raw_data, atoms_list):
    """
    retrieve the mass centers (x,y,z) of targeted atomic ensembles
    """
    mass_centers = []

    for data in raw_data:
        ensemble = Atoms([data[i] for i in atoms_list])
        mass_center = ensemble.get_center_of_mass(scaled=False) #If scaled=True the center of mass in scaled coordinates is returned (from ASE).
        mass_centers.append(mass_center)

    return mass_centers