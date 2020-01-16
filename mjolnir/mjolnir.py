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
parser.add_argument('pview',metavar='nview',nargs='*',help='Type of plot to make')
parser.add_argument("-f","--file",nargs=1,default=['results'],help='Results folder to use for plotting')
parser.add_argument("-s","--simulation_ID",nargs=1,default=['auto'],help='Name of simulation (e.g., planet name)')
parser.add_argument("-i","--initial_file",nargs=1,default=[10],type=int,help='Initial file id number (integer)')
parser.add_argument("-l","--last_file",nargs=1,default=['init'],type=int,help='Last file id number (integer)')
parser.add_argument("-lev","--horizontal_lev",nargs=1,default=[2.5e2],help='Horizonal level to plot in temperature/velocity/vorticity field (mbar or km)')
parser.add_argument("-vtop","--vertical_top",nargs=1,default=['default'],help='Location of top of plot (vertical type) in mbar (pressure) or fractional height (height)')
parser.add_argument("-slay","--split_layer",nargs=1,default=['no_split'],help='Split conserved quantities into weather and deep layers at this pressure')
parser.add_argument("-coord","--coordinate_sys",nargs=1,default=['icoh'],help='For KE spectrum, use either icoh grid or llp grid')
parser.add_argument("-ladj","--lmax_adjust",nargs=1,default=[0],help='For KE spectrum, icoh grid, adjust number of wave numbers to fit')
parser.add_argument("-slice","--slice",nargs='+',default=[0,360],help='Plot a long/lat slice or average over all values', type =float)
parser.add_argument("-mt","--maketable",action='store_true',help='Print a table in text file of plot data')
parser.add_argument("-no-log","--no_pressure_log",action='store_true',help='Switch off log coordinates in pressure (vertical)')
parser.add_argument("-llswap","--latlonswap",action='store_true',help='Swap latitude and longitude axes (horizontal plots)')
parser.add_argument("-vc","--vcoord",nargs=1,default=['pressure'],help='Vertical coordinate to use (pressure or height)')
parser.add_argument("-pgrid","--pgrid_ref",nargs=1,default=['auto'],help='Reference file for pressure grid')
args = parser.parse_args()
pview = args.pview

maketable = args.maketable
plog = True
if args.no_pressure_log:
    plog = False
valid = ['uver','ulonver','vver','wver','wlonver','wprof','Tver','Tlonver','Tulev','PTver','PTlonver','ulev','PVver','PVlev',
            'TP','RVlev','cons','stream','pause','tracer','PTP','regrid','KE',
            'SR','uprof','cfl','bvprof','fluxprof','Tsurf','insol','massf']

rg_needed = ['Tver','Tlonver','uver','ulonver','vver','wver','wlonver','Tulev','PTver','PTlonver','ulev','PVver','PVlev',
            'RVlev','stream','tracer','Tsurf','insol','massf']  #these types need regrid

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

if args.vcoord[0] == 'pressure':
    use_p = True
elif args.vcoord[0] == 'height':
    use_p = False
    plog = False  #enforce this on height grid
else:
    raise ValueError('%s not a valid vcoord. Valid options are "pressure" or "height"'%args.vcoord[0])

# regrid function is a special case and interrupts the other stuff
if 'regrid' in pview:
    ham.regrid(resultsf,simulation_ID,ntsi,nts,pressure_vert=use_p,pgrid_ref=args.pgrid_ref[0])
    exit()

outall = ham.GetOutput(resultsf,simulation_ID,ntsi,nts,openrg=openrg,
                        pressure_vert=use_p,pgrid_ref=args.pgrid_ref[0])

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

#--- Vertical plot types-------------------------------
if 'uver' in pview:
    z = {'value':rg.U, 'label':r'Velocity (m s$^{-1}$)', 'name':'u',
         'cmap':'viridis', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    # Averaged zonal winds (latitude vs pressure)
    #ham.u(input,grid,output,rg,sigmaref,slice=args.slice[0])
    ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p,csp=1000)
if 'ulonver' in pview:
    z = {'value':rg.U, 'label':r'Velocity (m s$^{-1}$)', 'name':'u',
         'cmap':'viridis', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    # Averaged zonal winds (latitude vs pressure)
    #ham.u(input,grid,output,rg,sigmaref,slice=args.slice[0])
    ham.vertical_lon(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p,csp=1000)
if 'vver' in pview:
    z = {'value':rg.V, 'label':r'Velocity (m s$^{-1}$)', 'name':'v',
         'cmap':'viridis', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    # Averaged zonal winds (latitude vs pressure)
    #ham.u(input,grid,output,rg,sigmaref,slice=args.slice[0])
    ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p)
if 'wver' in pview:
    z = {'value':rg.W, 'label':r'Velocity (m s$^{-1}$)', 'name':'w',
         'cmap':'viridis', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    # Averaged vertical winds (latitude vs pressure)
    #ham.w_ver(input,grid,output,rg,sigmaref)
    ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p,csp=1)
if 'wlonver' in pview:
    z = {'value':rg.W, 'label':r'Velocity (m s$^{-1}$)', 'name':'w',
         'cmap':'viridis', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    # Averaged zonal winds (latitude vs pressure)
    #ham.u(input,grid,output,rg,sigmaref,slice=args.slice[0])
    ham.vertical_lon(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p,csp=5)
if 'Tver' in pview:
    z = {'value':rg.Temperature, 'label':r'Temperature (K)', 'name':'temperature',
         'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    # Averaged temperature (latitude vs pressure)
    #ham.temperature(input,grid,output,rg,sigmaref)
    ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p,csp=[0])
if 'Tlonver' in pview:
    z = {'value':rg.Temperature, 'label':r'Temperature (K)', 'name':'temperature',
         'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    # Averaged temperature (latitude vs pressure)
    #ham.temperature(input,grid,output,rg,sigmaref)
    ham.vertical_lon(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p,csp=[0])
if 'PTver' in pview:
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient
    pt = rg.Temperature*(rg.Pressure/input.P_Ref)**(-kappa_ad)
    z = {'value':pt, 'label':r'Potential Temperature (K)', 'name':'potential_temp',
         'cmap':'plasma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    # Averaged potential temperature (latitude vs pressure)
    ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p,csp=5000)
if 'PTlonver' in pview:
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient
    pt = rg.Temperature*(rg.Pressure/input.P_Ref)**(-kappa_ad)
    z = {'value':pt, 'label':r'Potential Temperature (K)', 'name':'potential_temp',
         'cmap':'plasma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    # Averaged potential temperature (latitude vs pressure)
    ham.vertical_lon(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p,csp=5000)
if 'PVver' in pview:
    # sigmaref = np.arange(1,0,)
    z = {'value':rg.PV, 'label':r'Potential Vorticity (K m$^2$ kg$^{-1}$ s$^{-1}$)',
         'name':'pot_vort', 'cmap':'viridis', 'lat':rg.lat_lr, 'lon':rg.lon_lr, 'mt':maketable, 'plog':plog}
    sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
    ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p)
    # ham.potential_vort_vert(input,grid,output,sigmaref)
if 'stream' in pview: # RD: needs some work! to adapt to height coordinate
    # strm = ham.calc_moc_streamf(input,grid,output)
    # strm = rg.streamf
    # z = {'value':strm, 'label':r'Eulerian streamfunction (kg s$^{-1}$)', 'name':'streamf2',
    #      'cmap':'viridis', 'lat':rg.lat, 'lon':rg.lon}
    if use_p:
        sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
        # ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice,csp=[0])
        ham.streamf_moc_plot(input,grid,output,rg,sigmaref,mt=maketable,plog=plog)
    else:
        raise ValueError("'stream' plot type requires -vc pressure")
        # no reason to keep this way, just need to fix to use height
if 'massf' in pview:
    if use_p == False:
        # mass flow rate (zonal average)
        massfdl = (input.A+rg.Altitude[None,None,:,:])*np.cos(rg.lat[:,None,None,:]*np.pi/180)/(2*np.pi)*rg.Rho*rg.V
        z = {'value':massfdl, 'label':r'Mass flow', 'name':'massf',
             'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'plog':plog}
        sigmaref = ham.Get_Prange(input,grid,rg,args,xtype='lat',use_p=use_p)
        ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=args.slice,use_p=use_p,csp=[0])
    else:
        raise ValueError("'massf' plot type requires -vc height")

#--- Horizontal plot types-------------------------------
# need to be updated for height coordinates
if 'Tulev' in pview:
    # Averaged temperature and wind field (longitude vs latitude)
    # PR_LV - Pressure level (Pa)
    if use_p:
        PR_LV = np.float(args.horizontal_lev[0])*100
    else:
        PR_LV = np.float(args.horizontal_lev[0])*1000
    z = {'value':rg.Temperature, 'label':r'Temperature (K)', 'name':'temperature-uv',
            'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
if 'ulev' in pview:
    if use_p:
        PR_LV = np.float(args.horizontal_lev[0])*100
    else:
        PR_LV = np.float(args.horizontal_lev[0])*1000
    z = {'value':rg.U, 'label':r'Zonal Velocity (m s$^{-1}$)', 'name':'u',
        'cmap':'viridis', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
    z = {'value':rg.V, 'label':r'Meridional Velocity (m s$^{-1}$)', 'name':'v',
        'cmap':'viridis', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
if 'PVlev' in pview:
    if use_p:
        PR_LV = np.float(args.horizontal_lev[0])*100
    else:
        PR_LV = np.float(args.horizontal_lev[0])*1000
    z = {'value':rg.PV, 'label':r'Potential Vorticity (K m$^2$ kg$^{-1}$ s$^{-1}$)',
        'name':'pot_vort', 'cmap':'viridis', 'lat':rg.lat_lr, 'lon':rg.lon_lr, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
    # ham.potential_vort_lev(input,grid,output,PR_LV)
if 'RVlev' in pview:
    if use_p:
        PR_LV = np.float(args.horizontal_lev[0])*100
    else:
        PR_LV = np.float(args.horizontal_lev[0])*1000
    z = {'value':rg.RV[0], 'label':r'Relative Vorticity (s$^{-1}$)',
        'name':'rela_vort', 'cmap':'viridis', 'lat':rg.lat_lr, 'lon':rg.lon_lr, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
    # ham.rela_vort_lev(input,grid,output,PR_LV)
if 'tracer' in pview:
    if use_p:
        PR_LV = np.float(args.horizontal_lev[0])*100
    else:
        PR_LV = np.float(args.horizontal_lev[0])*1000
    z = {'value':np.log10(rg.ch4), 'label':r'Log(mixing ratio)',
        'name':'chem-ch4-uv1', 'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
    z = {'value':np.log10(rg.co), 'label':r'Log(mixing ratio)',
        'name':'chem-co-uv1', 'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
    z = {'value':np.log10(rg.h2o), 'label':r'Log(mixing ratio)',
        'name':'chem-h2o-uv1', 'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
    z = {'value':np.log10(rg.co2), 'label':r'Log(mixing ratio)',
        'name':'chem-co2-uv1', 'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
    z = {'value':np.log10(rg.nh3), 'label':r'Log(mixing ratio)',
        'name':'chem-nh3-uv1', 'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
if 'insol' in pview:
    PR_LV = np.max(output.Pressure) # not important here
    z = {'value':rg.insol, 'label':r'Insolation (W m$^{-2}$)', 'name':'insol',
        'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)
if 'Tsurf' in pview:
    if not hasattr(rg,"Tsurface"):
        raise ValueError("'Tsurface' not available in regrid file")
    PR_LV = np.max(output.Pressure) # not important here
    z = {'value':rg.Tsurface, 'label':r'Surface Temperature (K)', 'name':'Tsurf',
        'cmap':'magma', 'lat':rg.lat, 'lon':rg.lon, 'mt':maketable, 'llswap':args.latlonswap}
    ham.horizontal_lev(input,grid,output,rg,PR_LV,z,wind_vectors=True,use_p=use_p)

#--- Pressure profile types-------------------------------
if 'TP' in pview:
    z = {'value': output.Pressure/output.Rd/output.Rho, 'label':'Temperature (K)', 'name':'T' }
    #ham.TPprof(input,grid,output,sigmaref,1902)
    ham.profile(input,grid,output,z)
if 'PTP' in pview:
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient
    T = output.Pressure/input.Rd/output.Rho
    pt = T*(output.Pressure/input.P_Ref)**(-kappa_ad)
    z = {'value': pt, 'label':'Potential Temperature (K)', 'name':'PT' }
    ham.profile(input,grid,output,z)
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
if 'hsprof' in pview:
    dpdr = np.gradient(output.Pressure,grid.Altitude,axis=1)
    hseq = (dpdr + output.Rho*input.Gravit)/output.Rho/input.Gravit
    hseq[:,0,:] = np.nan
    hseq[:,-1,:] = np.nan
    z = {'value': hseq, 'label':r'$dP/dr + \rho g$', 'name':'hsprof'}
    ham.profile(input,grid,output,z,stride=20)
if 'bvprof' in pview:
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient
    T = output.Pressure/input.Rd/output.Rho
    pt = T*(output.Pressure/input.P_Ref)**(-kappa_ad)
    dptdr = np.gradient(pt,grid.Altitude,axis=1)
    N = np.sqrt(input.Gravit/pt*dptdr)
    z =  {'value': N, 'label':r'$N$ (s$^{-1}$)', 'name':'BVprof'}
    ham.profile(input,grid,output,z,stride=20)
if 'fluxprof' in pview:
    total_f = output.fnet_up - output.fnet_dn
    fup = total_f[:,:-1,:] + (total_f[:,1:,:]-total_f[:,:-1,:]) *\
          (grid.Altitude[None,:,None]-grid.Altitudeh[None,:-1,None])/\
          (grid.Altitudeh[None,1:,None]-grid.Altitudeh[None,:-1,None])
    z = {'value': fup, 'label':r'Total flux (W m$^{-2}$)', 'name':'ftot'}
    ham.profile(input,grid,output,z,stride=20)

#--- Global diagnostics -----------------------------------
if 'cons' in pview:  # RD: needs some work!
    if args.split_layer[0] == 'no_split':
        split = False
    else:
        split = np.float(args.split_layer[0])*100
    ham.conservation(input,grid,output,split)
if 'KE' in pview:  # RD: needs some work!
    PR_LV = np.float(args.horizontal_lev[0])*100 #not actually used here
    ham.KE_spect(input,grid,output,PR_LV,coord=args.coordinate_sys[0],lmax_adjust=args.lmax_adjust[0])
if 'SR' in pview:
    ham.SRindex(input,grid,output)
if 'RTbalance' in pview:
    ham.RTbalance(input,grid,output)

last = time.time()
print(last-first)
