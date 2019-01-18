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
#
# Options
#
# pview   plot
# uver    Averaged (time and longitude) zonal winds.
#         2D Map Latitude Vs Pressure.
#
# Tver    Averaged (time and longitude) temperatures.
#         2D Map Latitude Vs Pressure.
#
# Tulev   Averaged (time) temperature and wind field.
#         2D Map Longitude Vs Latitude.
#
# PTver   Averaged (time and longitude) potential temperatures.
#         2D Map Latitude Vs Pressure.
#
# ulev    Averaged (time) wind fields in long and lat.
#         2D Map Longitude Vs Latitude.
#
# PVver   Averaged (time and longitude) potential vorticity.
#         2D Map Latitude Vs Pressure.
#
# PVlev   Averaged (time) potential vorticity.
#         2D Map Longitude Vs Latitude.
###########################################################################

parser = argparse.ArgumentParser()
parser.add_argument('pview',metavar='nview',nargs='*',help='Type of plot to make (integer)')
parser.add_argument("-f","--file",nargs=1,default=['results'],help='Results folder to use for plotting')
parser.add_argument("-s","--simulation_ID",nargs=1,default=['auto'],help='Name of simulation (e.g., planet name)')
parser.add_argument("-i","--initial_file",nargs=1,default=[10],type=int,help='Initial file id number (integer)')
parser.add_argument("-l","--last_file",nargs=1,default=['init'],type=int,help='Last file id number (integer)')
parser.add_argument("-p","--pressure_lev",nargs=1,default=[2.5e2],help='Pressure level to plot in temperature/velocity/vorticity field (mbar)')
parser.add_argument("-pmin","--pressure_min",nargs=1,default=['default'],help='Lowest pressure value to plot in vertical plots')
parser.add_argument("-slay","--split_layer",nargs=1,default=['no_split'],help='Split conserved quantities into weather and deep layers at this pressure')
parser.add_argument("-coord","--coordinate_sys",nargs=1,default=['icoh'],help='For KE spectrum, use either icoh grid or llp grid')
parser.add_argument("-ladj","--lmax_adjust",nargs=1,default=[0],help='For KE spectrum, icoh grid, adjust number of wave numbers to fit')
parser.add_argument("-slice","--slice",nargs=1,default=['avg'],help='Plot a long/lat slice or average over all values')
args = parser.parse_args()
pview = args.pview

valid = ['uver','wver','wprof','Tver','Tulev','PTver','ulev','PVver','PVlev',
            'TP','RVlev','cons','stream','pause','tracer','PTP','regrid','KE','SR','uprof','cfl','dPW']

rg_needed = ['Tver','uver','wver','Tulev','PTver','ulev','PVver','PVlev',
            'RVlev','stream','tracer','PTP','KE']  #these types need regrid

openrg = 0
if 'all' in pview:
    pview = valid
    openrg = 1
else:
    for p in pview:
        if openrg == 0 and p in rg_needed:
            openrg = 1
        if p not in valid:
            raise ValueError('%s not a valid plot option. Valid options are '%p+', '.join(valid))

ntsi     = args.initial_file[0]  # initial file id number

if args.last_file[0] == 'init':
    nts = ntsi
else:
    nts = args.last_file[0]     # last file id number

if ntsi > nts:
    nts = ntsi

resultsf = args.file[0]
if args.simulation_ID[0] == 'auto':
    outname = spr.check_output('ls '+resultsf+'/esp_output_*_0.h5',shell=True)
    file0 = outname.decode().split(sep='/')[-1]
    simulation_ID = file0.split(sep='_')[2]
else:
    simulation_ID = args.simulation_ID[0]

# regrid function is a special case and interrupts the other stuff
if 'regrid' in pview:
    ham.regrid(resultsf,simulation_ID,ntsi,nts)
    exit()

outall = ham.GetOutput(resultsf,simulation_ID,ntsi,nts,openrg=openrg)

##########
# Planet #
##########

input = outall.input

########
# Grid #
########

grid = outall.grid

###############
# Diagnostics #
###############

output = outall.output

##################
# Regridded data #
##################
if openrg == 1:
    rg = outall.rg

# quick and dirty break point, if you just want to look at the output data
# does not load regridded data!
if 'pause' in pview:
    import pdb; pdb.set_trace()

#########
# Plots #
#########

# Sigma (normalized pressure) values for the plotting
if (args.pressure_min[0]=='default'):
    args.pressure_min[0] = np.max(output.Pressure[:,grid.nv-1,:])/100
if np.max(input.P_Ref)/np.float(args.pressure_min[0]) > 1000:
    sigmaref = np.logspace(np.log10(input.P_Ref),np.log10(np.float(args.pressure_min[0])*100),20)/input.P_Ref
else:
    sigmaref = np.linspace(input.P_Ref,np.float(args.pressure_min[0])*100,20)/input.P_Ref

#--- Vertical plot types-------------------------------
if 'uver' in pview:
    z = {'value':rg.U, 'label':r'Velocity (m s$^{-1}$)', 'name':'u', 'cmap':'viridis'}
    # Averaged zonal winds (latitude vs pressure)
    #ham.u(input,grid,output,rg,sigmaref,slice=args.slice[0])
    ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice[0])
if 'wver' in pview:
    z = {'value':rg.W, 'label':r'Velocity (m s$^{-1}$)', 'name':'w', 'cmap':'viridis'}
    # Averaged vertical winds (latitude vs pressure)
    #ham.w_ver(input,grid,output,rg,sigmaref)
    ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice[0])
if 'Tver' in pview:
    z = {'value':rg.Temperature, 'label':r'Temperature (K)', 'name':'temperature', 'cmap':'magma'}
    # Averaged temperature (latitude vs pressure)
    #ham.temperature(input,grid,output,rg,sigmaref)
    ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice[0])
if 'PTver' in pview:
    # Averaged potential temperature (latitude vs pressure)
    ham.potential_temp(input,grid,output,rg,sigmaref)
if 'PVver' in pview: # RD: needs some work!
    #sigmaref = np.arange(1,0,-0.05)
    ham.potential_vort_vert(input,grid,output,sigmaref)
if 'stream' in pview: # RD: needs some work!
    if np.max(input.P_Ref)/np.float(args.pressure_min[0]) > 1000:
        sigmaref = np.logspace(np.log10(input.P_Ref),np.log10(np.float(args.pressure_min[0])*100),grid.nv)/input.P_Ref
    else:
        sigmaref = np.linspace(input.P_Ref,np.float(args.pressure_min[0])*100,grid.nv)/input.P_Ref
    ham.streamf(input,grid,output,sigmaref)

#--- Horizontal plot types-------------------------------
if 'Tulev' in pview:
    # Averaged temperature and wind field (longitude vs latitude)
    # PR_LV - Pressure level (Pa)
    PR_LV = np.float(args.pressure_lev[0])*100
    ham.temperature_u_lev(input,grid,output,rg,PR_LV)
if 'ulev' in pview:
    PR_LV = np.float(args.pressure_lev[0])*100
    ham.uv_lev(input,grid,output,rg,PR_LV)
if 'PVlev' in pview:  # RD: needs some work!
    PR_LV = np.float(args.pressure_lev[0])*100
    ham.potential_vort_lev(input,grid,output,PR_LV)
if 'RVlev' in pview:  # RD: needs some work!
    PR_LV = np.float(args.pressure_lev[0])*100
    ham.rela_vort_lev(input,grid,output,PR_LV)
if 'tracer' in pview:  # RD: needs some work!
    PR_LV = np.float(args.pressure_lev[0])*100
    ham.tracer_u_lev(input,grid,output,PR_LV,'ch4')
    ham.tracer_u_lev(input,grid,output,PR_LV,'co')
    ham.tracer_u_lev(input,grid,output,PR_LV,'h2o')
    ham.tracer_u_lev(input,grid,output,PR_LV,'co2')
    ham.tracer_u_lev(input,grid,output,PR_LV,'nh3')

#--- Pressure profile types-------------------------------
if 'TP' in pview:
    z = {'value': output.Pressure/input.Rd/output.Rho, 'label':'Temperature (K)', 'name':'T' }
    #ham.TPprof(input,grid,output,sigmaref,1902)
    ham.profile(input,grid,output,z)
if 'PTP' in pview:
    ham.PTPprof(input,grid,output,sigmaref,1902)
if 'wprof' in pview:  # RD: needs some work!
    z = {'value': output.Wh[:,1:,:]/output.Rho, 'label':r'Vertical velocity (m s$^{-1}$)', 'name':'W' }
    ham.profile(input,grid,output,z,stride=20)
    # Averaged vertical winds (latitude vs pressure)
    # ham.w_prof(input,grid,output)
if 'uprof' in pview:  # RD: needs some work!
    u = (-output.Mh[0]*np.sin(grid.lon[:,None,None])+output.Mh[1]*np.cos(grid.lon[:,None,None]))/output.Rho
    z = {'value': u, 'label':r'Zonal velocity (m s$^{-1}$)', 'name':'U' }
    ham.profile(input,grid,output,z,stride=20)
if 'cfl' in pview:
    dt = output.time[0]/output.nstep[0]*86400
    dx = np.sqrt(np.min(grid.areasT))
    cs = np.sqrt(input.Cp/(input.Cp-input.Rd)*output.Pressure/output.Rho)
    z = {'value': cs*dt/dx, 'label':'CFL number for (horizontal) acoustic waves', 'name':'CFL' }
    ham.profile(input,grid,output,z,stride=20)
if 'dPW' in pview:
    x = grid.Altitude
    x0 = grid.Altitudeh[:-1]
    x1 = grid.Altitudeh[1:]
    y0 = output.Wh[:,:-1,:]
    y1 = output.Wh[:,1:,:]
    W = (y0*(x1-x)[None,:,None]+y1*(x-x0)[None,:,None])/(x1-x0)[None,:,None]
    dp = np.gradient(output.Pressure,axis=1)
    z = {'value':  dp*W/output.Rho, 'label':r'$\Delta P W$', 'name':'dPW' }
    ham.profile(input,grid,output,z,stride=20)


#--- Global diagnostics -----------------------------------
if 'cons' in pview:  # RD: needs some work!
    if args.split_layer[0] == 'no_split':
        split = False
    else:
        split = np.float(args.split_layer[0])*100
    ham.conservation(input,grid,output,split)
if 'KE' in pview:  # RD: needs some work!
    PR_LV = np.float(args.pressure_lev[0])*100
    ham.KE_spect(input,grid,output,rg,PR_LV,coord=args.coordinate_sys[0],lmax_adjust=args.lmax_adjust[0])
if 'SR' in pview:
    ham.SRindex(input,grid,output)

last = time.time()
print(last-first)
