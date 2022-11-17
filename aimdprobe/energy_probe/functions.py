###functions to analyze energies and other AIMD results###
import os
import json
import numpy as np
from ase import Atoms
from ase.io import read, write

def get_cumulative_avg(values):
    return np.cumsum(values) / np.arange(1, len(values) + 1)

def isfloat(value): ##remove bugged energies in calculations
    if value is None:
        return False
    try:
        float(value)
        return True
    except ValueError:
        return False


class EnergyAnalyzer:

    raw_energy_file: str
    ca_energy_file: str
    temperature_file: str

    def __init__(self):
        self.analyze_energy()
        self.store_energy()

    def get_energy(self, raw_data): #get potential energies from vasprun.xml
        potential_energy = []
        temperature = []
        for rd in raw_data:
            potential_energy.append(rd.get_potential_energy())
            temperature.append(rd.get_temperature())
        yield potential_energy, temperature

    def analyze_energy(self):
        tidy_energy = []
        for i, energy in enumerate(self.get_energy()):
            if isfloat(energy):
                tidy_energy.append(Energy)
            else:
                tidy_energy.append(self.get_energy()[0])
        
        self.tidy_energy = tidy_energy
        self.tidy_energy_ca = get_cumulative_avg(tidy_energy)
    
    def store_energy(self):
    # Save the file as a json
        with open(self.raw_energy_file, 'w') as handle:
            json.dump(self.tidy_energy, handle)
        with open(self.ca_energy_file,'w') as handle:
            json.dump(self.tidy_energy_ca, handle)
        with open(self.temperature_file,'w') as handle:
            json.dump(self.temperature, handle)

            

    
    
    


    


