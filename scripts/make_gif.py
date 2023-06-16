#!/usr/bin/env python
import numpy as np
import os
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import  read, write
from aimdprobe.init_data import init_data
from aimdprobe.other_probes.functions import animate


if __name__ == '__main__':
    """
    make gif animations of AIMD trajectory based on vasprun.xml

    """
    fp = os.getcwd()
    fn = 'vasprun.xml'
    
    # get raw data
    raw_data = init_data(fp, fn)

    # filename, images, intervals, size, timeframe, rotation
    filename = 'traj'
    step = 100
    timeframe = 200
    size = (2,2,1)
    rotation = '10z,-80x'
    animate(filename,raw_data, step, size, timeframe, rotation)

    # output 'traj.gif' file with 200 ms/frame