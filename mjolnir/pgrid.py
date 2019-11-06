#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import hamarr as ham
import sys
import argparse
import h5py
from imp import reload
reload(ham)
import time
import subprocess as spr

first = time.time()
###########################################################################
# Script function purely for setting up the pressure levels used in
# regridding.
# This makes it so that the pressure grid is uniform in time in the
# regridded files.
# Reads in files and computes the average pressure in each level, then
# averages them in time. Writes the resulting mean pressures
# to a plain text file.
###########################################################################


parser = argparse.ArgumentParser()
parser.add_argument('resultsf',metavar='nview',nargs='*',help='Results directory')
parser.add_argument("-s","--simulation_ID",nargs=1,default=['auto'],help='Name of simulation (e.g., planet name)')
parser.add_argument("-i","--initial_file",nargs=1,default=[10],type=int,help='Initial file id number (integer)')
parser.add_argument("-l","--last_file",nargs=1,default=['init'],type=int,help='Last file id number (integer)')
parser.add_argument("-stride","--stride",nargs=1,default=['init'],type=int,help='Stride (cadence) of files to use to set grid')
parser.add_argument("-w","--overwrite",action='store_true',help='force overwrite existing pgrid file')
args = parser.parse_args()
resultsf = args.resultsf[0]
ntsi     = args.initial_file[0]  # initial file id number

if args.last_file[0] == 'init':
    nts = ntsi
else:
    nts = args.last_file[0]     # last file id number

if ntsi > nts:
    nts = ntsi

if args.simulation_ID[0] == 'auto':
    outname = spr.check_output('ls '+resultsf+'/esp_output_*_0.h5',shell=True)
    file0 = outname.decode().split(sep='/')[-1]
    simulation_ID = file0.split(sep='_')[2]
else:
    simulation_ID = args.simulation_ID[0]

stride = args.stride[0]

ham.define_Pgrid(resultsf,simulation_ID,ntsi,nts,stride,overwrite=args.overwrite)

last = time.time()
print(last-first)
