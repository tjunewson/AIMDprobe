import os
import numpy as np
import json
from ase import Atoms
from ase.io import read, write

## Retrieve AIMD results from OUTCAR/vasprun.xml files

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
    creat a json file with a dictionary containing: 
    runtime, temperature, potential_energy, element, traj
    key: runtime
    source: OUTCAR
    """
    output_file = str

    def __init__(self,filepath, filename):
        self.data = init_data(filepath, filename, fmat='vasp')
        self.get_values()
        self.store_values()
    
    def get_data(self):
        """ summarize all values for the datafile """
        runtime = range(len(self.data))
        temperatures, temperatures_avg = get_temperatures() 
        pot_energy, pot_energy_avg = get_energy(self.data)
        elements = get_elements(self.data)
        raw_traj = get_raw_traj(self.data)

        #attention to json format
        values = [[str(elements), pot_energy[i], pot_energy_avg[i],raw_traj[i].tolist()] for i in runtime]
        yield runtime, values

    def store_data(self):
        
        data_dict = dict(zip(self.data.get_data()))
        json_object = json.dumps(data_dict)
        with open(self.output_file, "w") as handle:
            handle.write(json_object)
        print('A json file containing major raw_data has been created!')
