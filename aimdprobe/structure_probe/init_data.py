import os
import numpy as np
from ase import Atoms
from ase.io import read, write

## Retrieve AIMD results from vasprun.xml files

def init_data(filepath, filename, fmat='vasp'):
    fpath = filepath
    fname = filename ##e.g., vasprun.xml
    raw_data = read(fpath+'/'+fname,':')
    return raw_data

def get_raw_traj(raw_data):
    """
    function to retrieve positions of all structures
    """
    raw_traj=[]
    for rd in raw_data:
        raw_traj.append(rd.get_positions())
    return raw_traj