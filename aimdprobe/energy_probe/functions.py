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
        pot_energy.append(rd.get_potential_energy())
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
    
    
    


    


