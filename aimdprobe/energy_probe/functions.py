###functions to analyze energies and other AIMD results###
import os
import json
import numpy as np
from ase import Atoms
from ase.io import read, write
from aimdprobe.useful_functions import get_cumulative_avg

###### main functions ######
def get_energy(raw_data):
    """
    1. Get potential energies from vasprun.xml or outcar
    2. Note that here is only potential energy at 0K
    """
    pot_energy = []
    for rd in raw_data:
        pe = rd.get_potential_energy()
        if pe <= 0: # to remove 'INFINITY' energies in the traj
            pot_energy.append(pe)
    pot_energy_avg = get_cumulative_avg(pot_energy)
    return pot_energy, pot_energy_avg

def get_temperatures(outcar):
    """
    1. Get temperature in outcar in lines of 'kin. lattice  EKIN_LAT=         0.000000  (temperature  293.40 K)'
    2. Note that vasprun only shows results at 0K, thus no T and kinetic_energy are stored
    """
    with open(outcar, 'r') as file:
        temperatures = []
        for line in file:
            if 'EKIN_LAT=' in line.split():
                #print(line.split())
                numbs = [float(w) for w in line.split() if w[0].isdigit()] # get floats in the line
                t = numbs[1]
                temperatures.append(t)
 
    temperatures_avg = get_cumulative_avg(temperatures)
    return temperatures, temperatures_avg  

            
def get_kinetic_energy(outcar):
    """
    Get kinetic energies in outcar in lines of 'kinetic energy EKIN   =         6.068099'
    """
    with open(outcar, 'r') as file:
        kinetic_energy = []
        for line in file:
            if 'kinetic' and 'EKIN' in line.split():
                #print(line.split())
                ekin = [float(w) for w in line.split() if w[0].isdigit()] # get the float in the line
                kinetic_energy.append(ekin[0])
 
    kinetic_energy_avg = get_cumulative_avg(kinetic_energy)
    return kinetic_energy, kinetic_energy_avg  
    
def if_converge(energy, cutoff = 10000, stdev = 0.05):
    """
    check the convergence using the energies 
    in the cutoff region (e.g., last 10000 fs) with stdev (e.g., 0.02 eV) as the criterion
    """
    if len(energy) < cutoff:
        return 'NotLongEnough'
    
    elif np.std(energy[cutoff:]) < stdev:
        converged_energy = energy[-1]
        print('The converged energy is ' + str(converged_energy) + ' eV')
        return 'Converged' 
    
    else:
        print('The energy is not converged!')
        return 'NotConverged'



    
    


    


