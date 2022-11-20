import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from aimdprobe.useful_functions import get_solvent_traj, get_top_slab_mean_z_positions


###### main functions ######
"""
get the radial distrubution function for water solvent
"""
def get_rdf(raw_data, traj, ads_list, slab_list, nbins): ##counting solvents along z-axis with nbins density
    solvent_traj, solvent_symb, o_h2o_traj, h_h2o_traj = get_solvent_traj(raw_data, traj, ads_list, slab_list)
    sol_pos_z = [st[2] for st in solvent_traj]
    o_pos_z = [st[2] for st in o_h2o_traj]
    h_pos_z = [st[2] for st in h_h2o_traj]
    X, Y, Z = raw_data[0].get_cell()
    Z = Z[2]
    slab_top_z = get_top_slab_mean_z_positions(traj, slab_list)
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

    o_rdf = o_bins
    h_rdf = h_bins
    return walls, o_rdf, h_rdf

def get_rdf_3d(raw_data, traj, ads_list, slab_list, nbins): 
    
    """
    counting solvents in the space with nbins density
    """
    solvent_traj, solvent_symb, o_h2o_traj, h_h2o_traj = get_solvent_traj(raw_data, traj, ads_list, slab_list)
   
    X, Y, Z = raw_data[0].get_cell()
    X = X[0]
    Y = Y[1]
    Z = Z[2]
    slab_top_z = get_top_slab_mean_z_positions(traj, slab_list)
    assert slab_top_z < Z
    
    ## remove slab Z
    for o in o_h2o_traj:
        o[2] -= slab_top_z
    
    for h in h_h2o_traj:
        h[2] -= slab_top_z
    
    Z_norm = Z - slab_top_z
 
    #make nbins empty 3D boxes to put in atoms (e.g., O and H in H2O)
    #set bottom and top values for nbins
    walls_z = np.linspace(0, Z_norm, nbins+1) 
    walls_x = np.linspace(0, X, nbins+1)
    walls_y = np.linspace(0, Y, nbins+1)

    bins = (walls_x, walls_y, walls_z)
    lattice = make_lattice(walls_x[:-1], walls_y[:-1], walls_z[:-1])

    counts_o, bins_o = np.histogramdd(np.array(o_h2o_traj), bins = bins)
    counts_h, bins_h = np.histogramdd(np.array(h_h2o_traj), bins = bins)
    counts_o_f = counts_o.flatten()
    counts_h_f = counts_h.flatten()
    
    return counts_o_f, counts_h_f, lattice



