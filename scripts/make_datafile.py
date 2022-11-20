from aimdprobe import init_data

if __name__ == '__main__':

    method = init_data.MakeDatafile(
        input_file='OUTCAR',
        output_file='data.json',
    )