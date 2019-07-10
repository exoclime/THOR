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
import scipy.interpolate as interp


first = time.time()
###########################################################################
# Script function purely for regridding the icosahedral grid.
# Provides more control/options than mjolnir's regrid function
# -t, --type
#    type of horizontal interpolation to use:
#        gd or GD = scipy.interpolate.griddata (fastest, but noisy)
#        sh or SH = spherical harmonics from pySHTools (slower, but very smooth)
#        spl or SPL = scipy.interpolate.SmoothSphereBivariateSpline (currently not working:
#                     very fussy about smoothing factor and often produces artifacts near poles.
#                     if you've some idea how to fix, please let me know!)
###########################################################################


#things to add:
# lmax
# rotation of grid (for tidally locked planets, for example)
# force overwrite
parser = argparse.ArgumentParser()
parser.add_argument('resultsf',metavar='nview',nargs='*',help='Results directory')
parser.add_argument("-s","--simulation_ID",nargs=1,default=['auto'],help='Name of simulation (e.g., planet name)')
parser.add_argument("-t","--type",nargs=1,default=['gd'],choices=['gd','GD','sh','SH'],help='Horizontal interpolation type')
parser.add_argument("-vc","--vcoord",nargs=1,default=['pressure'],help='Vertical coordinate to use (pressure or height)')
parser.add_argument("-pmin","--pressure_min",nargs=1,default=['default'],help='Lowest pressure value in vertical interpolation (vc=pressure only)')
parser.add_argument("-i","--initial_file",nargs=1,default=[10],type=int,help='Initial file id number (integer)')
parser.add_argument("-l","--last_file",nargs=1,default=['init'],type=int,help='Last file id number (integer)')
args = parser.parse_args()
resultsf = args.resultsf[0]

ntsi     = args.initial_file[0]  # initial file id number

if args.last_file[0] == 'init':
    nts = ntsi
else:
    nts = args.last_file[0]     # last file id number

if ntsi > nts:
    nts = ntsi

# resultsf = args.file[0]
if args.simulation_ID[0] == 'auto':
    outname = spr.check_output('ls '+resultsf+'/esp_output_*_0.h5',shell=True)
    file0 = outname.decode().split(sep='/')[-1]
    simulation_ID = file0.split(sep='_')[2]
else:
    simulation_ID = args.simulation_ID[0]

if args.vcoord[0] == 'pressure':
    use_p = True
elif args.vcoord[0] == 'height':
    use_p = False
else:
    raise ValueError('%s not a valid vcoord. Valid options are "pressure" or "height"'%args.vcoord[0])

ham.regrid(resultsf,simulation_ID,ntsi,nts,pressure_vert=use_p,type=args.type[0],pressure_min=args.pressure_min[0])

last = time.time()
print(last-first)
