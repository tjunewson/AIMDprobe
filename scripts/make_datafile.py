"""
store data to a json file for postprocessing
"""
from aimdprobe import init_data

if __name__ == '__main__':

    method = init_data.MakeDatafile(
        input_file='OUTCAR', ## change to your file name
        output_file='data.json',
    )