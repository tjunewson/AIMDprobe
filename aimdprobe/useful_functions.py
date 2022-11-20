"""
useful functions for analyze energy and structure of AIMD simulations
"""

import os
import numpy as np
import csv, json
from ase import Atoms
from ase.io import read, write

def count_time(raw_data, time_step):
    """
    function to retrieve the AIMD simulation time
    """
    steps = len(raw_data)
    runtime = steps*time_step

    if time_step == 1:
        print('Runtime: ' + runtime + 'fs')
    elif time_step == 0.001:
        print('Runtime: ' + runtime + 'ps')
    else:
        print('Runtime: ' + runtime + 'ns')
    return steps, runtime

def get_elements(raw_data):
    all_elements = raw_data[0].symbols
    return all_elements

def make_lattice(a,b,c):
    """
    function to make a 3D lattice mesh with array a, b, c
    """
    lattice = []
    for aa in a:
        for bb in b:
            for cc in c:
                coord = [aa,bb,cc]
                lattice.append(coord)
    return lattice

def remove_zero(list, xx, yy, zz):
    """
    function to remove all zeros in list
    """
    zeros = [ind for ind in range(0,len(list)) if list[ind] == 0]
    list_nonzero = [l for l in list if l != 0]
    xx_nonzero = [x for i, x in enumerate(xx) if i not in zeros]
    yy_nonzero = [y for i, y in enumerate(yy) if i not in zeros]
    zz_nonzero = [z for i, z in enumerate(zz) if i not in zeros]
    return list_nonzero, xx_nonzero, yy_nonzero, zz_nonzero

def get_cumulative_avg(values):
    return np.cumsum(values) / np.arange(1, len(values) + 1)

def isfloat(value): 
    """
    remove bugged values in calculations
    """
    if value is None:
        return False
    else:
        try:
            float(value)
            return True
        except ValueError:
            return False

def get_real_distance(raw_data, a, b):
    """
    avoid boundary issue for distance calculation between two atoms a and b in the cell
    lattice constant = X, Y, Z
    """
    X,Y,Z = raw_data[0].get_cell()
    X = X[0], Y = Y[1], Z = Z[2]
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

def get_solvent_traj(raw_data, traj, ads_list, slab_list):
    """
    extract solvent structures for further analysis, e.g., H-bond counting or rdf
    """
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

def get_ads_mean_z_positions(traj, ads_list):
    ads_pos_z = [traj[i][2] for i in ads_list]
    mean_z_ads = np.mean(ads_pos_z)
    return mean_z_ads

def get_top_slab_mean_z_positions(traj, slab_list):
    slab_traj = [traj[i] for i in slab_list]
    surf_pos_z = np.array(slab_traj)[:,2] #get z coordinates for all slab atoms
    surf_pos_top_z = np.sort(surf_pos_z)[-1*len(slab_list):] # get the top-layer slab z coordinates
    mean_top_z = np.mean(surf_pos_top_z)
    return mean_top_z