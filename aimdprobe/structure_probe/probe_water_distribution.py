import os
import numpy as np
from ase import Atoms
from ase.io import read, write

#get the radial distrubution function for water solvent
def get_rdf(raw_data, traj, n_s, N_w, ads_list, slab_list, nbins): ##counting solvents along z-axis with nbins density
    solvent_traj, o_h2o_traj, h_h2o_traj = get_h2o_separate_traj(raw_data, traj, N_w, ads_list, slab_list)
    sol_pos_z = [st[2] for st in solvent_traj]
    o_pos_z = [st[2] for st in o_h2o_traj]
    h_pos_z = [st[2] for st in h_h2o_traj]
    X, Y, Z = raw_data[0].get_cell()
    Z = Z[2]
    slab_top_z = get_top_slab_mean_z_positions(traj, n_s, slab_list)
    assert slab_top_z < Z
    sol_pos_z_norm = [s - slab_top_z for s in sol_pos_z]
    o_pos_z_norm = [s - slab_top_z for s in o_pos_z]
    h_pos_z_norm = [s - slab_top_z for s in h_pos_z]
    Z_norm = Z - slab_top_z

    #make empty bins to put in atoms (e.g., O and H in H2O)
    walls = np.linspace(0, Z_norm, nbins) #set bottom and top values for nbins
    o_bins = np.zeros(nbins)
    h_bins = np.zeros(nbins)

    for i, bot in enumerate(walls[:nbins-1]): #to avoid overflow
        top = walls[i+1]
        for oz in o_pos_z_norm:
            if bot < oz <= top:
                o_bins[i] += 1
        for hz in h_pos_z_norm:
            if bot < hz <= top:
                h_bins[i] += 1

    o_rdf = [np.mean(i) for i in o_bins]
    h_rdf = [np.mean(i) for i in h_bins]
    return walls, o_rdf, h_rdf

