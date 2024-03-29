"""
functions to analyze structures of AIMD results

note that we assume the interface takes an order of surface-solvent-adsorbate
e.g., CO* on 4x4x4 Cu(111) surface Cu_64_O_20_H_40_C1_O1
"""

import os
import numpy as np
import csv, json
from ase import Atoms
from ase.io import read, write
from ase.geometry.analysis import Analysis


def get_elements(raw_data):
    all_elements = raw_data[0].symbols
    return all_elements

def get_elements_atoms(raw_data, atoms_list):
    all_elements = raw_data[0].symbols
    atoms_elements = [all_elements[i] for i in atoms_list]
    return atoms_elements

def store_raw_traj(raw_data, raw_traj): #store raw_traj coordinates to a csv file, with 1st row of elements
    all_elements = get_elements(raw_data)
    with open('raw_traj.csv','w',encoding = 'UTF8') as f:
        fwriter = csv.writer(f)
        fwriter.writerow(all_elements)
        fwriter.writerow(raw_traj)
    f.close()
    print('A file containing raw_traj has been saved!')
    return

def get_real_distance(raw_data, a, b):
    """
    avoid boundary issue for distance calculation between two atoms a and b in the cell
    lattice constant = X, Y, Z
    """
    X,Y,Z = raw_data[0].get_cell()
    X = X[0]
    Y = Y[1]
    Z = Z[2]
    c = [aa - bb for aa, bb in zip(a,b)]
    # avoid periodic distances
    for i, cc in enumerate(c):
        if i == 0:
            if np.abs(cc) > X:
                c[i] = np.abs(cc) - X
        elif i == 1:
            if np.abs(cc) > Y:
                c[i] = np.abs(cc) - Y
        else:
            if np.abs(cc) > Z:
                c[i] = np.abs(cc) - Z
    dist_a_b = np.linalg.norm(c)
    return dist_a_b


def get_separate_traj(traj,slab_list,ads_list):
    slab_traj = [traj[i] for i in slab_list]
    ads_traj = [traj[i] for i in ads_list]
    solvent_traj = [traj[i] for i in np.arange(len(traj)) if i not in slab_list if i not in ads_list]   
    return slab_traj, ads_traj, solvent_traj

def get_top_slab_mean_z_positions(traj, slab_list, n_s):
    slab_pos_z = [traj[i][2] for i in slab_list] #get z coordinates for all slab atoms
    surf_pos_top_z = np.sort(slab_pos_z)[-1*n_s:] # get the top-layer slab z coordinates
    mean_top_z = np.mean(surf_pos_top_z)
    return mean_top_z

def get_positions(traj, atoms_list):
    atoms_traj = [traj[i] for i in atoms_list]
    return atoms_traj

def get_solvent_traj(raw_data, traj, ads_list, slab_list):
    all_elements = get_elements(raw_data)
    solvent_traj = [traj[i] for i in np.arange(len(traj)) if i not in slab_list if i not in ads_list]   
    solvent_symb = [all_elements[i] for i in np.arange(len(all_elements)) if i not in slab_list if i not in ads_list]
    o_h2o_traj = []
    h_h2o_traj = []
    for i, s in enumerate(solvent_symb):
        if s == 'O':
            o_h2o_traj.append(solvent_traj[i])
        else:
            h_h2o_traj.append(solvent_traj[i]) 
    return solvent_traj, solvent_symb, o_h2o_traj, h_h2o_traj


def get_solvent_z_positions(traj, ads_list, slab_list): ###use water here
    sol_pos_z = []
    solvent_traj, solvent_symb = get_solvent_traj(traj, ads_list, slab_list)
    sol_pos_z.append(solvent_traj[:,2])
    return sol_pos_z



