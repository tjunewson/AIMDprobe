import os
import numpy as np
from ase import Atoms
from ase.io import read, write

# probe h-bonds formation between water and adsorbates
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

# probe h-bonds formation between water molecules