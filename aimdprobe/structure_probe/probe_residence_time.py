import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from aimdprobe.useful_functions import get_real_distance
from aimdprobe.structure_probe.functions import get_positions, get_elements_atoms

###### main functions ######
"""
probe the residence time of the possible adsorption
"""

def residence_time(raw_data, raw_traj, slab_list, adsorbate_list, cutoff):
    ## cutoff is in the form of {'atom_A':a, atom_B':b,...}
    ads_elements = get_elements_atoms(raw_data, adsorbate_list)
    rt = np.zeros(len(ads_elements))
    ads_elements_order = [a+str(i) for i, a in enumerate(ads_elements)]

    print(ads_elements_order)

    resi_time = {a:b for a, b in zip(ads_elements_order, rt)} 

    for traj in raw_traj:
        top_slabs = [traj[i] for i in slab_list]
        adsorbate = [traj[i] for i in adsorbate_list]
        for i, ads_atom in enumerate(adsorbate):
            ads_element = ads_elements[i]
            ads_element_order = ads_elements_order[i]
            for top_atom in top_slabs:
                if get_real_distance(raw_data, ads_atom, top_atom) <= cutoff[ads_element]:
                    resi_time[ads_element_order] += 1
    return resi_time







