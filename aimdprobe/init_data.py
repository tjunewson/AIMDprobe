import os
import numpy as np
import json
from dataclasses import dataclass
from ase import Atoms
from ase.io import read, write
from aimdprobe.energy_probe.functions import get_energy, get_temperatures, get_kinetic_energy
from aimdprobe.structure_probe.functions import get_elements

"""
Retrieve AIMD results from OUTCAR/vasprun.xml files
"""

def init_data(filepath, filename, fmat='vasp'):
    fpath = filepath
    fname = filename ##e.g., vasprun.xml
    raw_data = read(fpath+'/'+fname,':') #Atoms object
    return raw_data

def get_raw_traj(raw_data):
    """
    function to retrieve positions of all structures
    """
    raw_traj=[]
    for rd in raw_data:
        raw_traj.append(rd.get_positions())
    return raw_traj

class MakeDatafile:
    """
    create a json file with a dictionary containing:
    runtime, temperature, potential_energy, element, traj
    key: runtime
    source: OUTCAR (if not temperatures, you could use vasprun.xml as well)
    """

    filepath : str
    filename : str 
    output_file : str

    def __init__(self, filepath, filename, output_file):
        self.data = init_data(filepath, filename, fmat='vasp')
        self.output_file = output_file
        self.get_data(filepath, filename)
        self.store_data( filepath, filename)
        
    def get_data(self, filepath, filename):
        """ summarize all values for the datafile """
        runtime = np.arange(len(self.data))
        temperatures, temperatures_avg = get_temperatures(filepath+'/'+filename) 
        pot_energy, pot_energy_avg = get_energy(self.data)
        kin_energy, kin_energy_avg = get_kinetic_energy(filepath+'/'+filename)
        elements = get_elements(self.data)
        raw_traj = get_raw_traj(self.data)

        #attention to json format, but you could alsp use 'MyEncoder' to avoid errors afterwards
        values = [[str(elements), temperatures[i], pot_energy[i], kin_energy[i], raw_traj[i].tolist()] for i in runtime]
        yield runtime, values

    def store_data(self, filepath, filename):
        """ store data dictionaries for the datafile.json """

        self.data_dict = []
        for runtime, values in self.get_data(filepath, filename):
            self.data_dict.append(dict(zip(runtime.tolist(), values)))

        with open(self.output_file, 'w') as handle:
            json.dump(self.data_dict[0], handle)

        print('Congratulations, dear AIMDers! \n \n The json file of your AIMD results has been created!')
