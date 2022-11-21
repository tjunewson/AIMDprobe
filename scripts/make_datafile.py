"""
store major raw data from AIMD simulations to a json file for postprocessing
"""
from aimdprobe import init_data

if __name__ == '__main__':

    method = init_data.MakeDatafile(
        filepath= 'your_file_location',
        input_file='OUTCAR',
        output_file='data_name.json',
    )