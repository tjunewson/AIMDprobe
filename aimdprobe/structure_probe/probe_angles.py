import numpy as np
from ase import Atoms
from ase.io import read, write
from aimdprobe.useful_functions import get_cumulative_avg


###### main functions ######
"""
Probe the angles of certain atoms (>=3) and facets (>=2) in the simulated system
"""

def get_atomic_angles(raw_data, atoms_list):
    """
    Calculate angle in degrees between vectors between atoms a2->a1 and a2->a3, 
    where a1, a2, and a3 are in each row of indices.
    """

    assert len(atoms_list) == 3
    angles = []
    n = 3
    indices = []
    
    for i in range(n):
        for j in range(n):
            if j != i:
                for k in range(n):
                    if k != i and k!= j:
                        indices.append((atoms_list[i],atoms_list[j],atoms_list[k]))
    
    for traj in raw_data:
        angles.append(traj.get_angles(indices, mic=True))
    angles_ = np.transpose(angles)    
    angles_avg = [get_cumulative_avg(i) for i in angles_]
    return angles_avg


def get_dihetral_angles(raw_data, atoms_list):
    """
    calculate dihedral angle (in degrees) between the vectors a0->a1 and a2->a3
    """

    assert atoms_list == 4
    angles = []
    n = 4
    indices = []
    
    for i in range(n):
        for j in range(n):
            if j != i:
                for k in range(n):
                    if k != i and k != j:
                        for m in range(n):
                            if m != i and m != j and m != k:
                                indices.append((atoms_list[i],atoms_list[j],atoms_list[k],atoms_list[m]))
    
    for traj in raw_data:
        angles.append(traj.get_angles(indices, mic=True))
    angles_ = np.transpose(angles)    
    dihetral_angles_avg = [get_cumulative_avg(i) for i in angles_]

    return dihetral_angles_avg