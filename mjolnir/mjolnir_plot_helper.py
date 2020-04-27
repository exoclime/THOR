import hamarr as ham
import subprocess as spr
import numpy as np
import sys
import traceback

class mjol_args:
    # make_plot args, populate with defaults
    def __init__(self, file):
        self.pview = []
        self.file = [file]
        self.simulation_ID = ['auto']
        self.initial_file = [10]
        self.last_file = ['init']
        self.horizontal_lev = [2.5e2]
        self.vertical_top = ['default']
        self.split_layer = ['no_split']
        self.coordinate_sys = ['icoh']
        self.lmax_adjust = [0]
        self.slice = ['default']
        self.maketable = False
        self.no_pressure_log = False
        self.latlonswap = False
        self.vcoord = ['pressure']
        self.pgrid_ref = ['auto']
        self.clevels = [40]

def call_plot(name, func, *args, **kwargs):
    try:
        pfile = func(*args,**kwargs)
        if pfile is not None:
            print('Created file: ' + str(pfile))
    except:
        print(traceback.format_exc())
        print(f'{name} plot FAILED') 
        
    return pfile

def make_plot(args, save=True, axis=None):
    pview = args.pview

    maketable = args.maketable
    plog = True
    if args.no_pressure_log:
        plog = False

    valid = ['uver', 'ulonver', 'vver', 'wver', 'wlonver', 'wprof', 'Tver', 'Tlonver', 'Tulev', 'PTver', 'PTlonver', 'ulev', 'PVver', 'PVlev',
             'TP', 'RVlev', 'cons', 'stream', 'pause', 'tracer', 'PTP', 'regrid', 'KE',
             'SR', 'uprof', 'cfl', 'bvprof', 'TSfluxprof', 'Tsurf', 'insol', 'massf', 'pause_rg', 'DGfluxprof', 'qheat',
             'DGfutprof', 'DGfdtprof', 'DGfnetprof', 'mustar', 'DGfuptot', 'DGfdowntot', 'DGfnet', 'DGqheat',  # alf stuff
             'TSfutprof', 'TSfdtprof', 'TSfnetprof', 'mustar', 'TSfuptot', 'TSfdowntot', 'TSfnet', 'TSqheat',
             'DGqheatprof', 'TSqheatprof', 'qheatprof',
             'phase','all']

    rg_needed = ['Tver', 'Tlonver', 'uver', 'ulonver', 'vver', 'wver', 'wlonver', 'Tulev', 'PTver', 'PTlonver', 'ulev', 'PVver', 'PVlev',
                 'RVlev', 'stream', 'tracer', 'Tsurf', 'insol', 'massf', 'pause_rg',
                 'mustar', 'qheat',
                 'TSfuptot', 'TSfdowntot', 'TSfnet', 'TSqheat',
                 'DGfuptot', 'DGfdowntot', 'DGfnet', 'DGqheat',
                 'all']  # these types need regrid

    openrg = 0

    for p in pview:
        if openrg == 0 and p in rg_needed:
            openrg = 1
        if p not in valid:
            raise ValueError('%s not a valid plot option. Valid options are ' % p + ', '.join(valid))

    ntsi = args.initial_file[0]  # initial file id number

    if args.last_file[0] == 'init':
        nts = ntsi
    else:
        nts = args.last_file[0]     # last file id number

    if ntsi > nts:
        nts = ntsi

    resultsf = args.file[0]
    if args.simulation_ID[0] == 'auto':
        outname = spr.check_output('ls ' + resultsf + '/esp_output_*_0.h5', shell=True)
        file0 = outname.decode().split(sep='/')[-1]
        simulation_ID = file0.split(sep='_')[2]
    else:
        simulation_ID = args.simulation_ID[0]

    if args.vcoord[0] == 'pressure':
        use_p = True
    elif args.vcoord[0] == 'height':
        use_p = False
        plog = False  # enforce this on height grid
    else:
        raise ValueError('%s not a valid vcoord. Valid options are "pressure" or "height"' % args.vcoord[0])

    # regrid function is a special case and interrupts the other stuff
    if 'regrid' in pview:
        ham.regrid(resultsf, simulation_ID, ntsi, nts, pressure_vert=use_p, pgrid_ref=args.pgrid_ref[0])
        exit()

    outall = ham.GetOutput(resultsf, simulation_ID, ntsi, nts, openrg=openrg,
                           pressure_vert=use_p, pgrid_ref=args.pgrid_ref[0])

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
    if 'pause' in pview or 'pause_rg' in pview:
        import pdb
        pdb.set_trace()

    #########
    # Plots #
    #########

    plots_created = []
    # --- Vertical plot types-------------------------------
    if 'uver' in pview or 'all' in pview:
        rg.load(['U'])  #load these arrays into memory
        z = {'value': rg.U, 'label': r'Velocity (m s$^{-1}$)', 'name': 'u',
             'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('uver',ham.vertical_lat,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, csp=1000, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'ulonver' in pview or 'all' in pview:
        rg.load(['U'])
        z = {'value': rg.U, 'label': r'Velocity (m s$^{-1}$)', 'name': 'u',
             'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('ulonver',ham.vertical_lat,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, csp=1000, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'vver' in pview or 'all' in pview:
        rg.load(['V'])
        z = {'value': rg.V, 'label': r'Velocity (m s$^{-1}$)', 'name': 'v',
             'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('vver',ham.vertical_lat,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'wver' in pview or 'all' in pview:
        rg.load(['W'])
        z = {'value': rg.W, 'label': r'Velocity (m s$^{-1}$)', 'name': 'w',
             'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('wver',ham.vertical_lat,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, csp=1, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'wlonver' in pview or 'all' in pview:
        rg.load(['W'])
        z = {'value': rg.W, 'label': r'Velocity (m s$^{-1}$)', 'name': 'w',
             'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('wlonver',ham.vertical_lon,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, csp=5, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'Tver' in pview or 'all' in pview:
        rg.load(['Temperature'])
        z = {'value': rg.Temperature, 'label': r'Temperature (K)', 'name': 'temperature',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('Tver',ham.vertical_lat,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, csp=[0], clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'Tlonver' in pview or 'all' in pview:
        z = {'value': rg.Temperature, 'label': r'Temperature (K)', 'name': 'temperature',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('Tlonver',ham.vertical_lon,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, csp=[0], clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)
    if 'PTver' in pview or 'all' in pview:
        rg.load(['Temperature','Pressure'])
        kappa_ad = input.Rd / input.Cp  # adiabatic coefficient
        pt = rg.Temperature * (rg.Pressure / input.P_Ref)**(-kappa_ad)
        z = {'value': pt, 'label': r'Potential Temperature (K)', 'name': 'potential_temp',
             'cmap': 'plasma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('PTver',ham.vertical_lat,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, csp=5000, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'PTlonver' in pview or 'all' in pview:
        rg.load(['Temperature','Pressure'])
        kappa_ad = input.Rd / input.Cp  # adiabatic coefficient
        pt = rg.Temperature * (rg.Pressure / input.P_Ref)**(-kappa_ad)
        z = {'value': pt, 'label': r'Potential Temperature (K)', 'name': 'potential_temp',
             'cmap': 'plasma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('PTlonver',ham.vertical_lon,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, csp=5000, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'PVver' in pview or 'all' in pview:
        rg.load(['PV'])
        # sigmaref = np.arange(1,0,)
        z = {'value': rg.PV, 'label': r'Potential Vorticity (K m$^2$ kg$^{-1}$ s$^{-1}$)',
             'name': 'pot_vort', 'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
        sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
        pfile = call_plot('PVver',ham.vertical_lat,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'stream' in pview or 'all' in pview:  # RD: needs some work! to adapt to height coordinate
        # strm = ham.calc_moc_streamf(input,grid,output)
        # strm = rg.streamf
        # z = {'value':strm, 'label':r'Eulerian streamfunction (kg s$^{-1}$)', 'name':'streamf2',
        #      'cmap':'viridis', 'lat':rg.Latitude, 'lon':rg.Longitude}
        if use_p:
            rg.load(['V','W'])
            sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
            pfile = call_plot('stream',ham.streamf_moc_plot,input, grid, output, rg, sigmaref, mt=maketable, plog=plog, clevs=args.clevels, save=save)
        else:
            pfile = "'stream' plot type requires -vc pressure; plot not created"
            print(pfile)
        plots_created.append(pfile)
            # no reason to keep this way, just need to fix to use height

    if 'massf' in pview or 'all' in pview:
        if use_p == False:
            # mass flow rate (zonal average)
            rg.load(['Rho','V'])
            massfdl = (input.A + rg.Altitude[None, None, :, None]) * np.cos(rg.Latitude[:, None, None, None] * np.pi / 180) / (2 * np.pi) * rg.Rho * rg.V
            z = {'value': massfdl, 'label': r'Mass flow', 'name': 'massf',
                 'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'plog': plog}
            sigmaref = ham.Get_Prange(input, grid, rg, args, xtype='lat', use_p=use_p)
            pfile = call_plot('massf',ham.vertical_lat,input, grid, output, rg, sigmaref, z, slice=args.slice, use_p=use_p, csp=[0], clevs=args.clevels, save=save, axis=axis)
        else:
            pfile = "'massf' plot type requires -vc height; plot not created"
            print(pfile)
        plots_created.append(pfile)

    # --- Horizontal plot types-------------------------------
    # need to be updated for height coordinates
    if 'Tulev' in pview or 'all' in pview:
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        rg.load(['Temperature','U','V'])
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        z = {'value': rg.Temperature, 'label': r'Temperature (K)', 'name': 'temperature-uv',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('Tulev',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'ulev' in pview or 'all' in pview:
        rg.load(['U','V'])
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        z = {'value': rg.U, 'label': r'Zonal Velocity (m s$^{-1}$)', 'name': 'u',
             'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('ulev-U',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)
        z = {'value': rg.V, 'label': r'Meridional Velocity (m s$^{-1}$)', 'name': 'v',
             'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('ulev-V',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'PVlev' in pview or 'all' in pview:
        rg.load(['PV','U','V'])
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        z = {'value': rg.PV, 'label': r'Potential Vorticity (K m$^2$ kg$^{-1}$ s$^{-1}$)',
             'name': 'pot_vort', 'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('PVlev',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if 'RVlev' in pview or 'all' in pview:
        rg.load(['RVw','U','V'])
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        z = {'value': rg.RVw, 'label': r'Relative Vorticity (s$^{-1}$)',
             'name': 'rela_vort', 'cmap': 'viridis', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('RVlev',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if ('tracer' in pview or 'all' in pview) and input.chemistry:
        rg.load(['ch4','co','h2o','co2','nh3'])
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        z = {'value': np.log10(rg.ch4), 'label': r'Log(mixing ratio)',
             'name': 'chem-ch4-uv1', 'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('tracer-ch4',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)
        z = {'value': np.log10(rg.co), 'label': r'Log(mixing ratio)',
             'name': 'chem-co-uv1', 'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('tracer-co',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)
        z = {'value': np.log10(rg.h2o), 'label': r'Log(mixing ratio)',
             'name': 'chem-h2o-uv1', 'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('tracer-h2o',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)
        z = {'value': np.log10(rg.co2), 'label': r'Log(mixing ratio)',
             'name': 'chem-co2-uv1', 'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('tracer-co2',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)
        z = {'value': np.log10(rg.nh3), 'label': r'Log(mixing ratio)',
             'name': 'chem-nh3-uv1', 'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('tracer-nh3',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if ('insol' in pview or 'all' in pview) and (input.RT or input.TSRT):
        rg.load(['insol'])
        PR_LV = np.max(output.Pressure)  # not important here
        z = {'value': rg.insol, 'label': r'Insolation (W m$^{-2}$)', 'name': 'insol',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('insol',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if ('Tsurf' in pview or 'all' in pview) and (input.RT or input.TSRT):
        if not hasattr(rg, "Tsurface"):
            raise ValueError("'Tsurface' not available in regrid file")
        rg.load(['Tsurface'])
        PR_LV = np.max(output.Pressure)  # not important here
        z = {'value': rg.Tsurface, 'label': r'Surface Temperature (K)', 'name': 'Tsurf',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('Tsurf',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if ('DGfuptot' in pview or 'all' in pview) and input.RT:
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000

        rg.load(['flw_up'])
        fup = rg.flw_up
        z = {'value': fup, 'label': r'Double Gray Total upward flux (W m^-2)', 'name': 'DGfuptot',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('DGfuptot',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)
        
    if ('TSfuptot' in pview or 'all' in pview) and input.TSRT:
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000

        rg.load(['f_up_tot'])
        fup = rg.f_up_tot
        z = {'value': fup, 'label': r'Two Streams Total upward flux (W m^-2)', 'name': 'TSfuptot',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('TSfuptot',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

        
    if ('DGfdowntot' in pview or 'all' in pview) and input.RT: #add all later
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000

        rg.load(['flw_dn','fsw_dn'])
        fdn = rg.flw_dn + rg.fsw_dn

        z = {'value': fdn, 'label': r'Double Gray Total downward flux (W m^-2)', 'name': 'DGfdowntot',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('DGfdowntot',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if ('TSfdowntot' in pview or 'all' in pview) and input.TSRT: #add all later
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        rg.load(['f_down_tot'])
        fdn = rg.f_down_tot
        z = {'value': fdn, 'label': r'Two Streams Total downward flux (W m^-2)', 'name': 'TSfdowntot',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('TSfdowntot',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)


    if ('DGfnet' in pview or 'all' in pview)  and input.RT:
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        rg.load(['DGf_net'])
        z = {'value': rg.DGf_net, 'label': r'Double Gray Total net flux (W m^-2)', 'name': 'DGfnet',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('DGfnet',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if ('TSfnet' in pview or 'all' in pview)  and input.TSRT:
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        rg.load(['TSf_net'])
        z = {'value': rg.TSf_net, 'label': r'Two Streams Total net flux (W m^-2)', 'name': 'TSfnet',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('TSfnet',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)
        
    if ('qheat' in pview) and (input.RT or input.TSRT):
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        z = {'value': rg.qheat, 'label': r'Q Heat (W m^-3)', 'name': 'qheat',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('qheat',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if ('DGqheat' in pview) and input.RT:
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        z = {'value': rg.DGqheat, 'label': r'Double Gray Q Heat (W m^-3)', 'name': 'DGqheat',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('DGqheat',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    if ('TSqheat' in pview) and input.TSRT:
        # Averaged temperature and wind field (longitude vs latitude)
        # PR_LV - Pressure level (Pa)
        if use_p:
            PR_LV = np.float(args.horizontal_lev[0]) * 100
        else:
            PR_LV = np.float(args.horizontal_lev[0]) * 1000
        z = {'value': rg.TSqheat, 'label': r'Q Heat (W m^-3)', 'name': 'TSqheat',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('TSqheat',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

        
    if ('mustar' in pview) and (input.TSRT):
        PR_LV = np.max(output.Pressure)  # not important here
        z = {'value': rg.mustar, 'label': r'mu_star ', 'name': 'mustar',
             'cmap': 'magma', 'lat': rg.Latitude, 'lon': rg.Longitude, 'mt': maketable, 'llswap': args.latlonswap}
        pfile = call_plot('mustar',ham.horizontal_lev,input, grid, output, rg, PR_LV, z, wind_vectors=True, use_p=use_p, clevs=args.clevels, save=save, axis=axis)
        plots_created.append(pfile)

    # --- Pressure profile types-------------------------------
    if 'TP' in pview or 'all' in pview:
        output.load_reshape(grid,['Pressure','Rd','Rho'])
        z = {'value': output.Pressure / output.Rd / output.Rho, 'label': 'Temperature (K)', 'name': 'T'}
        # ham.TPprof(input,grid,output,sigmaref,1902)
        pfile = call_plot('TP',ham.profile,input, grid, output, z, save=save, axis=axis)
        plots_created.append(pfile)

    if 'PTP' in pview or 'all' in pview:
        output.load_reshape(grid,['Pressure','Rd','Rho','Cp'])
        # kappa_ad = input.Rd/input.Cp  # adiabatic coefficient
        kappa_ad = output.Rd / output.Cp
        T = output.Pressure / input.Rd / output.Rho
        pt = T * (output.Pressure / input.P_Ref)**(-kappa_ad)
        z = {'value': pt, 'label': 'Potential Temperature (K)', 'name': 'PT'}
        pfile = call_plot('PTP',ham.profile,input, grid, output, z, save=save, axis=axis)
        plots_created.append(pfile)

    if 'wprof' in pview or 'all' in pview:  # RD: needs some work!
        output.load_reshape(grid,['Wh','Pressure', 'Rho'])
        z = {'value': output.Wh[:, 1:, :] / output.Rho, 'label': r'Vertical velocity (m s$^{-1}$)', 'name': 'W'}
        pfile = call_plot('wprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

    if 'uprof' in pview or 'all' in pview:  # RD: needs some work!
        output.load_reshape(grid,['Mh','Pressure', 'Rho'])
        u = (-output.Mh[0] * np.sin(grid.lon[:, None, None]) + output.Mh[1] * np.cos(grid.lon[:, None, None])) / output.Rho
        z = {'value': u, 'label': r'Zonal velocity (m s$^{-1}$)', 'name': 'U'}
        pfile = call_plot('uprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

    if 'cfl' in pview or 'all' in pview:
        output.load_reshape(grid,['Pressure','Rho'])
        dt = output.time[0] / output.nstep[0] * 86400
        dx = np.sqrt(np.min(grid.areasT))
        cs = np.sqrt(input.Cp / (input.Cp - input.Rd) * output.Pressure / output.Rho)
        z = {'value': cs * dt / dx, 'label': 'CFL number for (horizontal) acoustic waves', 'name': 'CFL'}
        pfile = call_plot('cfl',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

    if 'bvprof' in pview or 'all' in pview:
        output.load_reshape(grid,['Pressure','Rho'])
        kappa_ad = input.Rd / input.Cp  # adiabatic coefficient
        T = output.Pressure / input.Rd / output.Rho
        pt = T * (output.Pressure / input.P_Ref)**(-kappa_ad)
        dptdr = np.gradient(pt, grid.Altitude, axis=1)
        N = np.sqrt(input.Gravit / pt * dptdr)
        z = {'value': N, 'label': r'$N$ (s$^{-1}$)', 'name': 'BVprof'}
        pfile = call_plot('bvprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

    if ('DGfluxprof' in pview or 'all' in pview) and input.RT:

        output.load_reshape(grid,['flw_up','fsw_dn','flw_dn'])
        fnet_int = output.flw_up - output.flw_dn - output.fsw_dn
        fnet = fnet_int[:, :-1, :] + (fnet_int[:, 1:, :] - fnet_int[:, :-1, :]) *\
            (grid.Altitude[None, :, None] - grid.Altitudeh[None, :-1, None]) /\
            (grid.Altitudeh[None, 1:, None] - grid.Altitudeh[None, :-1, None])
        z = {'value': fnet, 'label': r'Double Gray Net flux (W m$^{-2}$)', 'name': 'DGfnet'}
        pfile = call_plot('DGfluxprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)
        
    if ('TSfluxprof' in pview or 'all' in pview) and input.TSRT:

        output.load_reshape(grid,['f_net'])
        fnet_int = output.f_net
        fnet = fnet_int[:, :-1, :] + (fnet_int[:, 1:, :] - fnet_int[:, :-1, :]) *\
            (grid.Altitude[None, :, None] - grid.Altitudeh[None, :-1, None]) /\
            (grid.Altitudeh[None, 1:, None] - grid.Altitudeh[None, :-1, None])
        z = {'value': fnet, 'label': r'Two Stream Net flux (W m$^{-2}$)', 'name': 'TSfnet'}
        pfile = call_plot('TSfluxprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

    if ('DGfutprof' in pview or 'all' in pview) and input.RT:

        output.load_reshape(grid,['flw_up'])
        fup_int = output.flw_up
        fup = fup_int[:, :-1, :] + (fup_int[:, 1:, :] - fup_int[:, :-1, :]) *\
            (grid.Altitude[None, :, None] - grid.Altitudeh[None, :-1, None]) /\
            (grid.Altitudeh[None, 1:, None] - grid.Altitudeh[None, :-1, None])
        z = {'value': fup, 'label': r'Double Gray Total Upward flux (W m$^{-2}$)', 'name': 'DGfuptot'}
        pfile = call_plot('DGfutprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

        
    if ('TSfutprof' in pview or 'all' in pview) and input.TSRT:
        output.load_reshape(grid,['f_up_tot'])
        fup_int = output.f_up_tot
        fup = fup_int[:, :-1, :] + (fup_int[:, 1:, :] - fup_int[:, :-1, :]) *\
            (grid.Altitude[None, :, None] - grid.Altitudeh[None, :-1, None]) /\
            (grid.Altitudeh[None, 1:, None] - grid.Altitudeh[None, :-1, None])
        z = {'value': fup, 'label': r'Two Stream Total Upward flux (W m$^{-2}$)', 'name': 'TSfuptot'}
        pfile = call_plot('TSfutprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

    if ('DGfdtprof' in pview or 'all' in pview) and input.RT:
        output.load_reshape(grid,['flw_dn','fsw_dn'])
        fdn_int = output.flw_dn + output.fsw_dn
        fdn = fdn_int[:, :-1, :] + (fdn_int[:, 1:, :] - fdn_int[:, :-1, :]) *\
            (grid.Altitude[None, :, None] - grid.Altitudeh[None, :-1, None]) /\
            (grid.Altitudeh[None, 1:, None] - grid.Altitudeh[None, :-1, None])
        z = {'value': fdn, 'label': r'Double Gray Total Downward flux (W m$^{-2}$)', 'name': 'DGfdowntot'}
        pfile = call_plot('DGfdtprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

    if ('TSfdtprof' in pview or 'all' in pview) and input.TSRT:
        output.load_reshape(grid,['f_down_tot'])
        fdn_int = output.f_down_tot
        fdn = fdn_int[:, :-1, :] + (fdn_int[:, 1:, :] - fdn_int[:, :-1, :]) *\
            (grid.Altitude[None, :, None] - grid.Altitudeh[None, :-1, None]) /\
            (grid.Altitudeh[None, 1:, None] - grid.Altitudeh[None, :-1, None])
        z = {'value': fdn, 'label': r'Two Stream Total Downward flux (W m$^{-2}$)', 'name': 'TSfdowntot'}
        pfile = call_plot('TSfdtprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

    if ('TSqheatprof' in pview or 'all' in pview) and input.TSRT:
        output.load_reshape(grid,['TSqheat'])
        qheat = output.TSqheat
        z = {'value': qheat, 'label': r'Two Stream Q heat (W m$^{-2}$)', 'name': 'TSqheatprof'}
        pfile = call_plot('TSqheatprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)

    if ('qheatprof' in pview or 'all' in pview):
        output.load_reshape(grid,['qheat'])
        qheat = output.qheat
        z = {'value': qheat, 'label': r'Q heat (W m$^{-2}$)', 'name': 'qheatprof'}
        pfile = call_plot('qheatprof',ham.profile,input, grid, output, z, stride=20, save=save, axis=axis)
        plots_created.append(pfile)
    # --- Global diagnostics -----------------------------------
    if 'cons' in pview:  # RD: needs some work!
        if args.split_layer[0] == 'no_split':
            split = False
        else:
            split = np.float(args.split_layer[0]) * 100
        ham.conservation(input, grid, output, split)

    if 'KE' in pview:  # RD: needs some work!
        PR_LV = np.float(args.horizontal_lev[0]) * 100  # not actually used here
        ham.KE_spect(input, grid, output, PR_LV, coord=args.coordinate_sys[0], lmax_adjust=args.lmax_adjust[0])

    if 'SR' in pview:
        ham.SRindex(input, grid, output)

    if 'RTbalance' in pview:
        ham.RTbalance(input, grid, output)
        
    if 'phase' in pview and input.RT:
        ham.phase_curve(input,grid,output)

    #some clean up
    output.closeVDS()
    if openrg:
        rg.closeVDS()

    return plots_created
