###functions to analyze structures of AIMD results###

"""
note that we assume the interface takes an order of surface-solvent-adsorbate
e.g., CO* on 4x4x4 Cu(111) surface Cu_64_O_20_H_40_C1_O1
"""

"""
N_s, n_s, N_w, N_ads -> 
number of slab atoms, number of top slab atoms, number of solvents, number of adsorbate atoms
"""
import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.geometry.analysis import Analysis

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


def get_elements(raw_data):
    all_elements = raw_data[0].symbols
    return all_elements

def store_raw_traj(raw_data, raw_traj): #store raw_traj coordinates to a csv file
    all_elements = get_elements(raw_data)
    with open('raw_traj.csv','w',encoding = 'UTF8') as f:
        fwriter = csv.writer(f)
        fwriter.writerow(all_elements)
        fwriter.writerow(raw_traj)
    f.close()
    return 'Notice: raw_traj saved!'



## Structural analysis



def get_separate_traj(traj,slab_list,ads_list):
    slab_traj = [traj[i] for i in slab_list]
    ads_traj = [traj[i] for i in ads_list]
    solvent_traj = [traj[i] for i in np.arange(len(traj)) if i not in slab_list if i not in ads_list]   
    return slab_traj, ads_traj, solvent_traj

def get_top_slab_mean_z_positions(traj, n_s, slab_list):###
    slab_traj = [traj[i] for i in slab_list]
    surf_pos_z = np.array(slab_traj)[:,2] #get z coordinates for all slab atoms
    surf_pos_top_z = np.sort(surf_pos_z)[-1*n_s:] # get the top-layer slab z coordinates
    mean_top_z = np.mean(surf_pos_top_z)
    return mean_top_z

def get_solvent_traj(raw_data, traj, ads_list, slab_list):
    all_elements = get_elements(raw_data)
    solvent_traj = [traj[i] for i in np.arange(len(traj)) if i not in slab_list if i not in ads_list]   
    solvent_symb = [all_elements[i] for i in np.arange(len(all_elements)) if i not in slab_list if i not in ads_list]   
    return solvent_traj, solvent_symb

def get_h2o_separate_traj(raw_data, traj, N_w, ads_list, slab_list):
    solvent_traj, solvent_symb = get_solvent_traj(raw_data, traj, ads_list, slab_list)
    
    o_h2o_traj = []
    h_h2o_traj = []
    for i, s in enumerate(solvent_symb):
        if s == 'O':
            o_h2o_traj.append(solvent_traj[i])
        else:
            h_h2o_traj.append(solvent_traj[i])

    return solvent_traj, o_h2o_traj, h_h2o_traj

def get_solvent_z_positions(traj, ads_list, slab_list): ###use water here
    sol_pos_z = []
    solvent_traj, solvent_symb = get_solvent_traj(raw_data, traj, ads_list, slab_list)
    sol_pos_z.append(solvent_traj[:,2])
    return sol_pos_z



def get_rdf_solvent(raw_data, elements, rmax, nbins, N_s): 
    """
    We monitor the metal-water interface via a radial distribution function (RDF) 
    between the top-layer metal atoms and 
    the hydrogen (H) and oxygen atoms (O) of the water molecules.
    This function returns a dictionary of RDFs of each element with the key of element.
    """
    ana = Analysis(raw_data)
    rdf = {}
    for e in elements:
        raw_rdf = ana.get_rdf(rmax, nbins, imageIdx=None, elements = e, return_dists=False)
        mean_rdf = sum(raw_rdf)/(len(raw_rdf)*N_s)
        rdf[e] = mean_rdf
    
    return rdf


def get_adsorbed_h2o(raw_data,traj, N_w, n_s, ads_list, slab_list, dist): ##N_w_ads/N_s = N_w_ads site-1
    """
    retrieve the cumulative number of adsorbed solvents by calculating the
    z-direction distance between solvent molecule and top-slab plane;
    here, we use O in H2O to estimate the position of water molecules;
    """
    solvent_traj, o_h2o_traj, h_h2o_traj = get_h2o_separate_traj(raw_data,traj, N_w, ads_list, slab_list)
    sol_pos_z = [t[2] for t in o_h2o_traj]
    mean_top_z = get_top_slab_mean_z_positions(traj, n_s, slab_list)
    count = 0
    for spz in sol_pos_z:
        if abs(abs(spz - mean_z_top) - dist) <= 0.05:
            count += 1
        else:
            count += 0
        N_w_ads.append(count)
    return np.mean(N_w_ads)/n_s ##normalized to per site

def get_adsorbed_solvent_new(raw_traj, n_s, N_w, dist): ##N_w_ads/N_s = N_w_ads site-1
    """
    retrieve the cumulative number of adsorbed solvents by calculating the
    z-direction distance between solvent molecule and top-slab plane;
    here, we use O and H in H2O to roughly estimate the position of water molecules;
    and in a water molecule, if either O or H has a distance within dist to the top 
    surface, we count it as adsorbed water;
    Note: the atomic number (AN) of O in H2O should be in line with the Hs in the same H2O,
    e.g., O is the 15th in list of O_AN, Hs should be 30th and 31th in H_AN.
    """

    sol_pos_z = get_solvent_z_positions_new(raw_traj,N_w,N_s)
    mean_z_top = get_top_slab_mean_z_positions(raw_traj,N_s,n_s)
    assert len(sol_pos_z) == len(mean_z_top)

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



def get_adsorbed_solvent_mass(raw_data, N_w, N_s, n_s, dist): ##N_w_ads/N_s = N_w_ads site-1
    """
    retrieve the cumulative number of adsorbed solvents by calculating the
    z-direction distance between mass center of solvent molecule and top-slab plane;
    here, we use O in H2O to roughly estimate the position of water molecules;
    1. return individual H2O traj;
    2. get_center_of_mass(); 
    3. if z_mass_center < dist, then count += 1
    """
    
    N_w_ads = []
    #sol_pos_z = get_solvent_z_positions_new(raw_traj,N_w,N_s)
    mean_z_top = get_top_slab_mean_z_positions(raw_traj,N_s,n_s)
    #assert len(sol_pos_z) == len(mean_z_top)
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

# def check_honds_angles(a, b, a_elem, N_s, N_w, N_ads, angle): #a is the atom in adsorbate, b is the atom in water
#     if a_elem == 'O':

#     else:
def get_ads_mean_z_positions(raw_traj, N_ads):
    mean_z_ads = []
    for traj in raw_traj:
        ads_pos_z = traj[-1*N_ads:,2]
        mean_z_ads.append(np.mean(ads_pos_z))
    return mean_z_ads

def get_hbonds_ads(raw_traj, raw_data, N_s, N_w, N_ads, h_dist): #e.g., h_dist = 2.55+-0.05 for HO--H
    """
    function to retrieve numbers of hydrogen bonds of adsorbate with solvents;
    criteria for H-bonds: O-O distances / cutoff 3.5 AA / O-H-O angle > 140*;
    e.g., FCHO* - H2O
    """
    print('Start counting h-bonds between adsorbate and water')
    hbonds_ads = []
    #slab elements: [Au 64, O 42, C 5, H 84]
    for traj in raw_traj:
        # get adsorbate atom and element
        ads_traj = []
        ads_traj.append(traj[N_s:N_s+2]) # O in adsorbate
        ads_traj.append(traj[N_s+N_w+2:N_s+N_w+2+5+4]) # C and H in adsorbate
        ads_traj = [item for items in ads_traj for item in items]
        ads_symb = raw_data[0].get_chemical_symbols()[N_s:N_s+2] +  raw_data[0].get_chemical_symbols()[N_s+N_w+2:N_s+N_w+2+5+4]# list(ads_traj.symbols)
        ads_symb = [item for items in ads_symb for item in items]
        # get water atom and element
        O_traj = traj[N_s+2:N_s+2+N_w]
        H_traj = traj[-N_w*2:]
        count = 0
        for i, symb in enumerate(ads_symb):
            if symb == 'O':
                for h_pos in H_traj:
                    if abs(np.linalg.norm(ads_traj[i] - h_pos) - h_dist) <= 0.05:
                        count += 1
            elif symb == 'H':
                for o_pos in O_traj:
                    if abs(np.linalg.norm(ads_traj[i] - o_pos) - h_dist) <= 0.05:
                        count += 1
            else:
                count += 0
        hbonds_ads.append(count)
    #print(hbonds_ads)
    return hbonds_ads, np.mean(hbonds_ads)



                    


               

















