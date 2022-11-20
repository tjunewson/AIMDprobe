import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from aimdprobe.useful_functions import get_solvent_traj, get_top_slab_mean_z_positions


###### main functions ######
def get_adsorbed_h2o(raw_data,traj, ads_list, slab_list, dist): ##N_w_ads/N_s = N_w_ads site-1
    """
    retrieve the cumulative number of adsorbed solvents by calculating the
    z-direction distance between solvent molecule and top-slab plane;
    here, we use O in H2O to estimate the position of water molecules;
    """
    solvent_traj, solvent_symb, o_h2o_traj, h_h2o_traj = get_solvent_traj(raw_data, traj, ads_list, slab_list)
    sol_pos_z = [t[2] for t in o_h2o_traj] #use O to represent H2O
    mean_top_z = get_top_slab_mean_z_positions(traj, slab_list)
    n_s = len(slab_list)
    count = 0
    for spz in sol_pos_z:
        if abs(abs(spz - mean_z_top) - dist) <= 0.05:
            count += 1
        else:
            count += 0
        N_w_ads.append(count)
    return np.mean(N_w_ads)/n_s ##normalized to per site

def get_adsorbed_solvent_new(raw_data,traj, ads_list, slab_list, dist): ##N_w_ads/N_s = N_w_ads site-1
    """
    retrieve the cumulative number of adsorbed solvents by calculating the
    z-direction distance between solvent molecule and top-slab plane;
    here, we use O and H in H2O to roughly estimate the position of water molecules;
    and in a water molecule, if either O or H has a distance within dist to the top 
    surface, we count it as adsorbed water;
    Note: the atomic number (AN) of O in H2O should be in line with the Hs in the same H2O,
    e.g., O is the 15th in list of O_AN, Hs should be 30th and 31th in H_AN.
    """

    solvent_traj, solvent_symb, o_h2o_traj, h_h2o_traj = get_solvent_traj(raw_data, traj, ads_list, slab_list)
    sol_pos_z = [t[2] for t in o_h2o_traj] #use O to represent H2O
    mean_top_z = get_top_slab_mean_z_positions(traj, slab_list)
    n_s = len(slab_list)
    N_w = len(sol_pos_z)
    assert len(sol_pos_z) == len(mean_top_z)

    N_w_ads = []
    for i,spz in enumerate(sol_pos_z):
        count = 0
        o_an=[]
        for s in spz[:N_w]: ###O z coordinate
            if s - mean_z_top[i] <= dist:
                count += 1
                o_an.append(i)
        
        h_list = list(range(N_w,N_w*3))
        h_an_out = [i+N_w for i in o_an]
        h_an_out_ = [i+N_w+1 for i in o_an]
        h_an_out.extend(h_an_out_)
        h_list = [i for i in h_list if i not in h_an_out]
        ##remove counted H2O based on O-surface distance
        for hi in h_list:
            if hi % 2 == 0:
                if spz[hi]-mean_z_top[i] <= dist*0.8:
                    count += 1
            else:
                if spz[hi]-mean_z_top[i] <= dist*0.8 and spz[hi-1]-mean_z_top[i] > dist*0.8:
                    count += 1
        N_w_ads.append(count)

    return np.mean(N_w_ads)/n_s


### still in test ###
def get_adsorbed_solvent_mass(raw_data,traj, ads_list, slab_list, dist): ##N_w_ads/N_s = N_w_ads site-1
    """
    retrieve the cumulative number of adsorbed solvents by calculating the
    z-direction distance between mass center of solvent molecule and top-slab plane;
    here, we use O in H2O to roughly estimate the position of water molecules;
    1. return individual H2O traj;
    2. get_center_of_mass(); 
    3. if z_mass_center < dist, then count += 1
    """
    solvent_traj, solvent_symb, o_h2o_traj, h_h2o_traj = get_solvent_traj(raw_data, traj, ads_list, slab_list)
    sol_pos_z = [t[2] for t in o_h2o_traj] #use O to represent H2O
    mean_top_z = get_top_slab_mean_z_positions(traj, slab_list)
    N_s = len(slab_list)
    N_w_ads = []
    
    assert len(sol_pos_z) == len(mean_z_top)
    for id, data in enumerate(raw_data):
        count = 0
        traj = data.get_positions()[N_s:] ##get the xyz positions of solvents
        X, Y, Z = data.get_cell() ##get the x, y of the unit cell size

        #symb = data.get_chemical_symbols()[N_s:]
        """
        slice positions for O and H respectively
        """
        pos_o = traj[:N_w]
        pos_h = traj[N_w:N_w*3]
        #print(mean_z_top[id])
        for i, pos_o in enumerate(pos_o): ##O positions
            for j, pos_h_1 in enumerate(pos_h): ##H positions
                list(pos_h).pop(j)
                h2o_ = []
                for k, pos_h_2 in enumerate(pos_h):
                    """
                    consider the boundary effect of H-O-H, 
                    when H and O locate at different unit cells,
                    thus if d(O-H)x or y > len(unit cell), d(O-H)x or y += -len(unit cell)
                    """
                    h2o_.append([i,j,k,get_real_dist(pos_o,pos_h_1,X,Y,Z)+get_real_dist(pos_o,pos_h_2,X,Y,Z)])
                h2o_ = sorted(h2o_, key=lambda x: x[3])
                o,h1,h2 = h2o_[0][:3]
                h2o = Atoms('H2O', positions = [tuple(traj[a]) for a in [h1,h2,o]])
                if h2o.get_center_of_mass()[2] - mean_z_top[id] <= dist:
                    #print(h2o.get_center_of_mass()[2])
                    count += 1
        N_w_ads.append(count/2) ##in counting process, each H2O was counted twice, so /2 is needed here.

    return np.mean(N_w_ads)/n_s








                    


               








