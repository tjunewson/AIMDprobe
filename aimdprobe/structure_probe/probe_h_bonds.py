import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from aimdprobe.useful_functions import get_elements, get_solvent_traj, get_cumulative_avg
from aimdprobe.init_data import get_raw_traj

"""
probe h-bonds formation between water - adsorbates/ water - water
""" 

###### main functions ######
def get_hbonds_ads(raw_data, traj, slab_list, ads_list, h_dist): #e.g., h_dist = 2.55+-0.05 for HO--H
    """
    function to retrieve numbers of hydrogen bonds of adsorbate with solvents;
    criteria for H-bonds: O-O distances / cutoff 2.55 AA / O-H-O angle > 140*;
    e.g., FCHO* - H2O
    """
    print('Start counting h-bonds between adsorbate and water')

    #e.g., slab elements: [Au 64, O 42, C 5, H 84]
        
    # get adsorbate atom and element
    all_elements = get_elements(raw_data)
    ads_symb = [all_elements[i] for i in ads_list]
    # get water atom and element
    solvent_traj, solvent_symb, o_h2o_traj, h_h2o_traj = get_solvent_traj(raw_data, traj, ads_list, slab_list)
    count = 0
    ads_atom_to_form_hbond_o = []
    ads_atom_to_form_hbond_h = []

    for i, symb in enumerate(ads_symb):
        if symb == 'O':
            for h_pos in h_h2o_traj:
                if abs(np.linalg.norm(ads_traj[i] - h_pos) - h_dist) <= 0.05:
                    count += 1
                    ads_atom_to_form_hbond_o.append(i)
        elif symb == 'H':
            for o_pos in o_h2o_traj:
                if abs(np.linalg.norm(ads_traj[i] - o_pos) - h_dist) <= 0.05:
                    count += 1
                    ads_atom_to_form_hbond_h.append(i)
        else:
            count += 0
    hbonds_ads = count
    return hbonds_ads, ads_atom_to_form_hbond_o, ads_atom_to_form_hbond_h

def get_hbonds_ads_all(raw_data, raw_traj, traj, slab_list, ads_list, h_dist):

    raw_traj = get_raw_traj(raw_data)

    hbonds_ads_all = 0
    ads_atom_to_form_hbond_o_all = 0
    ads_atom_to_form_hbond_h_all = 0
    for traj in raw_traj:
        hbonds_ads, ads_atom_to_form_hbond_o, ads_atom_to_form_hbond_h = get_hbonds_ads(raw_data, traj, slab_list, ads_list, h_dist)
        hbonds_ads_all += hbonds_ads
        ads_atom_to_form_hbond_o_all += ads_atom_to_form_hbond_o
        ads_atom_to_form_hbond_h_all += ads_atom_to_form_hbond_h
    hbonds_ads_avg = get_cumulative_avg(hbonds_ads_all)
    ads_atom_to_form_hbond_o_all_avg = get_cumulative_avg(ads_atom_to_form_hbond_o_all)
    ads_atom_to_form_hbond_h_all_avg = get_cumulative_avg(ads_atom_to_form_hbond_h_all)
    return hbonds_ads_avg, ads_atom_to_form_hbond_o_all_avg, ads_atom_to_form_hbond_h_all_avg

# probe h-bonds formation between water molecules
def get_hbonds_sol(raw_data, traj, slab_list, ads_list, h_dist): #e.g., h_dist = 2.55+-0.05 for HO--H
    """
    function to retrieve numbers of hydrogen bonds in solvents;
    criteria for H-bonds: O-O distances / cutoff 2.55 AA / O-H-O angle > 140*;
    O-H bond length shoud be longer than a covalent O-H bond, cutoff 1.5 AA
    """
    print('Start counting h-bonds between waters')

    #e.g., slab elements: [Au 64, O 42, C 5, H 84]
        
    # get adsorbate atom and element
    all_elements = get_elements(raw_data)
    ads_symb = [all_elements[i] for i in ads_list]

    # get water atom and element
    solvent_traj, solvent_symb, o_h2o_traj, h_h2o_traj = get_solvent_traj(raw_data, traj, ads_list, slab_list)
    X,Y,Z = raw_data[0].get_cell()
    X = X[0], Y = Y[1], Z = Z[2]
    count = 0
    for i, o_pos in enumerate(o_h2o_traj):
        for h_pos in h_h2o_traj:
            dist_o_h = get_real_distance(o_pos, h_pos)
            if abs(dist_o_h - h_dist) <= 0.05 and dist_o_h - 1.5 >= 0.05:
                    count += 1
    hbonds_sol = count
    return hbonds_sol