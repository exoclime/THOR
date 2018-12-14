import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.interpolate as interp
import os
import h5py
import time
import subprocess as spr
import pyshtools as chairs

plt.rcParams['image.cmap'] = 'magma'
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'Helvetica-Normal'

class input:
    def __init__(self,resultsf,simID):
        fileh5 = resultsf+'/esp_output_planet_'+simID+'.h5'
        if os.path.exists(fileh5):
            openh5 = h5py.File(fileh5)
        else:
            fileh5_old = resultsf+'/esp_output_'+simID+'.h5'
            if os.path.exists(fileh5_old):
                openh5 = h5py.File(fileh5_old)
            else:
                raise IOError(fileh5+' or '+fileh5_old+' not found!')
        self.A = openh5['A'][...]  #'A' is the key for this dataset
        self.Rd = openh5['Rd'][...]  #[...] is syntax for "gimme all the data under this key"
        self.Omega = openh5['Omega'][...]
        self.P_Ref = openh5['P_Ref'][...]
        self.Top_altitude = openh5['Top_altitude'][...]
        self.Cp = openh5['Cp'][...]
        self.Gravit = openh5['Gravit'][...]
        if 'spring_beta' in openh5.keys():
            self.spring_beta = openh5['spring_beta'][...]
            self.vlevel = openh5['vlevel'][...]
            self.glevel = openh5['glevel'][...]
            self.spring_dynamics = openh5['spring_dynamics'][...]
        self.resultsf = resultsf
        self.simID = simID

        if 'SpongeLayer' in openh5.keys():
            self.SpongeLayer = openh5['SpongeLayer'][...]
            if self.SpongeLayer[0] == 1:
                self.nlat = openh5['nlat'][...]
                self.ns_sponge = openh5['ns_sponge'][...]
                self.Rv_sponge = openh5['Rv_sponge'][...]
            if 'hstest' in openh5.keys():
                self.core_benchmark = openh5['hstest'][...]
            else:
                self.core_benchmark = openh5['core_benchmark'][...]
            if self.core_benchmark[0] == 0: # need to switch to rt flag
                self.Tstar = openh5['Tstar'][...]
                self.planet_star_dist = openh5['planet_star_dist'][...]
                self.radius_star = openh5['radius_star'][...]
                self.diff_fac = openh5['diff_fac'][...]
                self.Tlow = openh5['Tlow'][...]
                self.albedo = openh5['albedo'][...]
                self.tausw = openh5['tausw'][...]
                self.taulw = openh5['taulw'][...]
        if 'vulcan' in openh5.keys():
            self.chemistry = openh5['vulcan'][...]
        if 'chemistry' in openh5.keys():
            self.chemistry = openh5['chemistry'][...]

        openh5.close()

class grid:
    def __init__(self,resultsf,simID):
        fileh5 = resultsf+'/esp_output_grid_'+simID+'.h5'
        if os.path.exists(fileh5):
            openh5 = h5py.File(fileh5)
        else:
            raise IOError(fileh5+' not found!')
        self.Altitude = openh5['Altitude'][...]
        self.Altitudeh = openh5['Altitudeh'][...]
        self.areasT = openh5['areasT'][...]
        self.lonlat = openh5['lonlat'][...]
        self.point_num = np.int(openh5['point_num'][0])
        self.pntloc = openh5['pntloc'][...]
        self.nv = np.int(openh5['nv'][0])
        openh5.close()
        self.nvi = self.nv+1
        self.lon = (self.lonlat[::2])%(2*np.pi)
        self.lat = self.lonlat[1::2]

class output:
    def __init__(self,resultsf,simID,ntsi,nts,grid,stride=1):
        # Initialize arrays
        self.Rho = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.Pressure = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.Mh = np.zeros((3,grid.point_num,grid.nv,nts-ntsi+1))
        self.Wh = np.zeros((grid.point_num,grid.nvi,nts-ntsi+1))
        self.ntsi = ntsi
        self.nts = nts
        self.time = np.zeros(nts-ntsi+1)
        self.nstep = np.zeros(nts-ntsi+1)
        self.Etotal = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.Mass = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.AngMomx = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.AngMomy = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.AngMomz = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.GlobalE = np.zeros(nts-ntsi+1)
        self.GlobalMass = np.zeros(nts-ntsi+1)
        self.GlobalAMx = np.zeros(nts-ntsi+1)
        self.GlobalAMy = np.zeros(nts-ntsi+1)
        self.GlobalAMz = np.zeros(nts-ntsi+1)
        self.ConvData = np.zeros(nts-ntsi+1)
        self.Insol = np.zeros((grid.point_num,nts-ntsi+1))

        self.ch4 = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.co = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.h2o = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.co2 = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.nh3 = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))

        self.tau_sw = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.tau_lw = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.fnet_up = np.zeros((grid.point_num,grid.nvi,nts-ntsi+1))
        self.fnet_dn = np.zeros((grid.point_num,grid.nvi,nts-ntsi+1))

        # Read model results
        for t in np.arange(ntsi-1,nts,stride):
            fileh5 = resultsf+'/esp_output_'+simID+'_'+np.str(t+1)+'.h5'
            if os.path.exists(fileh5):
                openh5 = h5py.File(fileh5)
            else:
                raise IOError(fileh5+' not found!')

            Rhoi = openh5['Rho'][...]
            Pressurei = openh5['Pressure'][...]
            Mhi = openh5['Mh'][...]
            Whi = openh5['Wh'][...]
            time = openh5['simulation_time'][0]/86400
            nstep = openh5['nstep'][0]
            if 'Etotal' in openh5.keys():
                Etotali = openh5['Etotal'][...]
                Massi = openh5['Mass'][...]
                AngMomxi = openh5['AngMomx'][...]
                AngMomyi = openh5['AngMomy'][...]
                AngMomzi = openh5['AngMomz'][...]
                self.GlobalE[t-ntsi+1] = openh5['GlobalE'][0]
                self.GlobalMass[t-ntsi+1] = openh5['GlobalMass'][0]
                self.GlobalAMx[t-ntsi+1] = openh5['GlobalAMx'][0]
                self.GlobalAMy[t-ntsi+1] = openh5['GlobalAMy'][0]
                self.GlobalAMz[t-ntsi+1] = openh5['GlobalAMz'][0]
                self.ConvData[t-ntsi+1] = True
            else:
                print('Warning: conservation diagnostics not available in file %s'%fileh5)
                self.ConvData[t-ntsi+1] = False
            if 'insol' in openh5.keys():
                self.Insol[:,t-ntsi+1] = openh5['insol'][...]
            if 'tracer' in openh5.keys():
                traceri = openh5['tracer'][...]
            if 'tau' in openh5.keys():
                taui = openh5['tau'][...]
                fupi = openh5['fnet_up'][...]
                fdni = openh5['fnet_dn'][...]
            openh5.close()

            self.Rho[:,:,t-ntsi+1] = np.reshape(Rhoi,(grid.point_num,grid.nv))
            self.Pressure[:,:,t-ntsi+1] = np.reshape(Pressurei,(grid.point_num,grid.nv))
            self.Mh[0,:,:,t-ntsi+1] = np.reshape(Mhi[::3],(grid.point_num,grid.nv))
            self.Mh[1,:,:,t-ntsi+1] = np.reshape(Mhi[1::3],(grid.point_num,grid.nv))
            self.Mh[2,:,:,t-ntsi+1] = np.reshape(Mhi[2::3],(grid.point_num,grid.nv))
            self.Wh[:,:,t-ntsi+1] = np.reshape(Whi,(grid.point_num,grid.nvi))
            self.time[t-ntsi+1] = time
            self.nstep[t-ntsi+1] = nstep
            if 'Etotali' in locals():
                self.Etotal[:,:,t-ntsi+1] = np.reshape(Etotali,(grid.point_num,grid.nv))
                self.Mass[:,:,t-ntsi+1] = np.reshape(Massi,(grid.point_num,grid.nv))
                self.AngMomx[:,:,t-ntsi+1] = np.reshape(AngMomxi,(grid.point_num,grid.nv))
                self.AngMomy[:,:,t-ntsi+1] = np.reshape(AngMomyi,(grid.point_num,grid.nv))
                self.AngMomz[:,:,t-ntsi+1] = np.reshape(AngMomzi,(grid.point_num,grid.nv))

            if 'traceri' in locals():
                self.ch4[:,:,t-ntsi+1] = np.reshape(traceri[::5],(grid.point_num,grid.nv))/self.Rho[:,:,t-ntsi+1]
                self.co[:,:,t-ntsi+1] = np.reshape(traceri[1::5],(grid.point_num,grid.nv))/self.Rho[:,:,t-ntsi+1]
                self.h2o[:,:,t-ntsi+1] = np.reshape(traceri[2::5],(grid.point_num,grid.nv))/self.Rho[:,:,t-ntsi+1]
                self.co2[:,:,t-ntsi+1] = np.reshape(traceri[3::5],(grid.point_num,grid.nv))/self.Rho[:,:,t-ntsi+1]
                self.nh3[:,:,t-ntsi+1] = np.reshape(traceri[4::5],(grid.point_num,grid.nv))/self.Rho[:,:,t-ntsi+1]

            if 'taui' in locals():
                self.tau_sw[:,:,t-ntsi+1] = np.reshape(taui[::2],(grid.point_num,grid.nv))
                self.tau_lw[:,:,t-ntsi+1] = np.reshape(taui[1::2],(grid.point_num,grid.nv))
                self.fnet_dn[:,:,t-ntsi+1] = np.reshape(fdni,(grid.point_num,grid.nvi))
                self.fnet_up[:,:,t-ntsi+1] = np.reshape(fupi,(grid.point_num,grid.nvi))

class rg_out:
    def __init__(self,resultsf,simID,ntsi,nts,input,grid):
        RT = 0
        if input.core_benchmark[0] == 0: # need to switch to 'radiative_tranfer' flag
            RT = 1

        # Read model results
        for t in np.arange(ntsi-1,nts):
            fileh5 = resultsf+'/regrid_'+simID+'_'+np.str(t+1)+'.h5'
            if os.path.exists(fileh5):
                openh5 = h5py.File(fileh5)
            else:
                print(fileh5+' not found, regridding now with default settings...')
                regrid(resultsf,simID,ntsi,nts)
                openh5 = h5py.File(fileh5)

            Rhoi = openh5['Rho'][...]
            Tempi = openh5['Temperature'][...]
            Ui = openh5['U'][...]
            Vi = openh5['V'][...]
            Wi = openh5['W'][...]
            Prei = openh5['Pressure'][...]
            lati = openh5['Latitude'][...]
            loni = openh5['Longitude'][...]

            if RT == 1:
                tau_swi = openh5['tau_sw'][...]
                tau_lwi = openh5['tau_lw'][...]
                fnet_upi = openh5['fnet_up'][...]
                fnet_dni = openh5['fnet_dn'][...]

            openh5.close()

            if t == ntsi-1:
                self.Rho = np.zeros(np.shape(Rhoi)+(nts-ntsi+1,))
                self.U = np.zeros(np.shape(Ui)+(nts-ntsi+1,))
                self.V = np.zeros(np.shape(Vi)+(nts-ntsi+1,))
                self.W = np.zeros(np.shape(Wi)+(nts-ntsi+1,))
                self.Temperature = np.zeros(np.shape(Tempi)+(nts-ntsi+1,))
                self.Pressure = np.zeros(np.shape(Prei)+(nts-ntsi+1,))
                self.lat = np.zeros(np.shape(lati)+(nts-ntsi+1,))
                self.lon = np.zeros(np.shape(loni)+(nts-ntsi+1,))
                if RT == 1:
                    self.tau_sw = np.zeros(np.shape(tau_swi)+(nts-ntsi+1,))
                    self.tau_lw = np.zeros(np.shape(tau_lwi)+(nts-ntsi+1,))
                    self.fnet_up = np.zeros(np.shape(fnet_upi)+(nts-ntsi+1,))
                    self.fnet_dn = np.zeros(np.shape(fnet_dni)+(nts-ntsi+1,))

            self.Rho[:,:,:,t-ntsi+1] = Rhoi
            self.U[:,:,:,t-ntsi+1] = Ui
            self.V[:,:,:,t-ntsi+1] = Vi
            self.W[:,:,:,t-ntsi+1] = Wi
            self.Temperature[:,:,:,t-ntsi+1] = Tempi
            self.Pressure[:,t-ntsi+1] = Prei
            self.lat[:,t-ntsi+1] = lati
            self.lon[:,t-ntsi+1] = loni
            if RT == 1:
                self.tau_sw[:,:,:,t-ntsi+1] = tau_swi
                self.tau_lw[:,:,:,t-ntsi+1] = tau_lwi
                self.fnet_up[:,:,:,t-ntsi+1] = fnet_upi
                self.fnet_dn[:,:,:,t-ntsi+1] = fnet_dni

class GetOutput:
    def __init__(self,resultsf,simID,ntsi,nts,stride=1,openrg=0):
        self.input = input(resultsf,simID)
        self.grid = grid(resultsf,simID)
        self.output = output(resultsf,simID,ntsi,nts,self.grid,stride=stride)
        if openrg == 1:
            self.rg = rg_out(resultsf,simID,ntsi,nts,self.input,self.grid)

def regrid(resultsf,simID,ntsi,nts,res_deg=0.5,nlev=40,pscale='log',overwrite=False,comp=4):
    # runs over files and converts ico-height grid to lat-lon-pr grid
    outall = GetOutput(resultsf,simID,ntsi,nts)
    input = outall.input
    grid = outall.grid
    output = outall.output

    print('Regrid data in folder '+resultsf+'...\n')

    #figure out pressure grid
    pmin = np.min(output.Pressure)
    if pscale == 'log':
        sigmaref = np.logspace(np.log10(input.P_Ref),np.log10(pmin),nlev)/input.P_Ref
    elif pscale == 'lin':
        sigmaref = np.linspace(input.P_Ref,pmin,nlev)/input.P_Ref
    else:
        raise IOError('invalid pressure scale entered! use "lin" or "log"')

    d_sig = np.size(sigmaref)
    Pref = input.P_Ref*sigmaref
    lat_range_tmp = np.arange(-90,90+res_deg,res_deg)
    # recenter so that there are an even number of latitude points
    lat_range = (lat_range_tmp[:-1]+lat_range_tmp[1:])/2
    lon_range = np.arange(0,360,res_deg)
    loni, lati = np.meshgrid(lon_range,lat_range)
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    Temp_icoh = output.Pressure/(input.Rd*output.Rho)
    interpx = (grid.Altitude-grid.Altitudeh[:-1])/(grid.Altitudeh[1:]-grid.Altitudeh[:-1])
    # on even height grid, interpolation is excessive, but wth?
    Wh_icoh = output.Wh[:,:-1,:] + (output.Wh[:,1:,:]-output.Wh[:,:-1,:])*interpx[None,:,None]

    RT = 0
    if hasattr(input,"core_benchmark"):
        if input.core_benchmark[0] == 0: # need to switch to 'radiative_tranfer' flag
            RT = 1
            fnet_up_icoh = output.fnet_up[:,:-1,:]+(output.fnet_up[:,1:,:]-output.fnet_up[:,:-1,:])*interpx[None,:,None]
            fnet_dn_icoh = output.fnet_dn[:,:-1,:]+(output.fnet_dn[:,1:,:]-output.fnet_dn[:,:-1,:])*interpx[None,:,None]

    for t in np.arange(ntsi,nts+1):
        #check for exising h5 files
        proceed = 0
        fileh5 = resultsf+'/regrid_'+simID+'_'+np.str(t)+'.h5'
        if os.path.exists(fileh5):
            if overwrite == True:
                proceed = 1
            else:
                print(fileh5+' already present! Skipping time = %d'%t)
        else:
            proceed = 1

        if proceed == 1:
            print('Regridding time = %d...'%t)
            Temp_icop = np.zeros((grid.point_num,d_sig))
            Temp_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            Rho_icop = np.zeros((grid.point_num,d_sig))
            Rho_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            Mh_icop = np.zeros((3,grid.point_num,d_sig))

            Wh_icop = np.zeros((grid.point_num,d_sig))

            U_llp = np.zeros((d_lon[0],d_lon[1],d_sig))
            V_llp = np.zeros((d_lon[0],d_lon[1],d_sig))
            W_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            if RT == 1:
                tau_sw_icop = np.zeros((grid.point_num,d_sig))
                tau_sw_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                tau_lw_icop = np.zeros((grid.point_num,d_sig))
                tau_lw_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                fnet_up_icop = np.zeros((grid.point_num,d_sig))
                fnet_up_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                fnet_dn_icop = np.zeros((grid.point_num,d_sig))
                fnet_dn_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            for i in np.arange(grid.point_num):
                #interp to pressure grid
                sigma = output.Pressure[i,:,t-ntsi]
                Temp_icop[i,:] = interp.pchip_interpolate(sigma[::-1],Temp_icoh[i,::-1,t-ntsi],Pref[::-1])[::-1]
                Rho_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.Rho[i,::-1,t-ntsi],Pref[::-1])[::-1]
                Mh_icop[0,i,:] = interp.pchip_interpolate(sigma[::-1],output.Mh[0,i,::-1,t-ntsi],Pref[::-1])[::-1]
                Mh_icop[1,i,:] = interp.pchip_interpolate(sigma[::-1],output.Mh[1,i,::-1,t-ntsi],Pref[::-1])[::-1]
                Mh_icop[2,i,:] = interp.pchip_interpolate(sigma[::-1],output.Mh[2,i,::-1,t-ntsi],Pref[::-1])[::-1]
                Wh_icop[i,:] = interp.pchip_interpolate(sigma[::-1],Wh_icoh[i,::-1,t-ntsi],Pref[::-1])[::-1]
                if RT == 1:
                    tau_sw_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.tau_sw[i,::-1,t-ntsi],Pref[::-1])[::-1]
                    tau_lw_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.tau_lw[i,::-1,t-ntsi],Pref[::-1])[::-1]
                    fnet_up_icop[i,:] = interp.pchip_interpolate(sigma[::-1],fnet_up_icoh[i,::-1,t-ntsi],Pref[::-1])[::-1]
                    fnet_dn_icop[i,:] = interp.pchip_interpolate(sigma[::-1],fnet_dn_icoh[i,::-1,t-ntsi],Pref[::-1])[::-1]

            # Convert icosahedral grid into lon-lat grid
            U_icop = (-Mh_icop[0]*np.sin(grid.lon[:,None])+Mh_icop[1]*np.cos(grid.lon[:,None]))/Rho_icop
            V_icop = (-Mh_icop[0]*np.sin(grid.lat[:,None])*np.cos(grid.lon[:,None])\
                      -Mh_icop[1]*np.sin(grid.lat[:,None])*np.sin(grid.lon[:,None])\
                      +Mh_icop[2]*np.cos(grid.lat[:,None]))/Rho_icop

            for lev in np.arange(d_sig):
                Temp_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,Temp_icop[:,lev],(loni,lati),method='nearest')
                Rho_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,Rho_icop[:,lev],(loni,lati),method='nearest')
                U_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,U_icop[:,lev],(loni,lati),method='nearest')
                V_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,V_icop[:,lev],(loni,lati),method='nearest')
                W_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,Wh_icop[:,lev]/Rho_icop[:,lev],(loni,lati),method='nearest')
                if RT == 1:
                    tau_sw_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,tau_sw_icop[:,lev],(loni,lati),method='nearest')
                    tau_lw_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,tau_lw_icop[:,lev],(loni,lati),method='nearest')
                    fnet_up_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,fnet_up_icop[:,lev],(loni,lati),method='nearest')
                    fnet_dn_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,fnet_dn_icop[:,lev],(loni,lati),method='nearest')

            #create h5 files
            openh5 = h5py.File(fileh5,"w")

            print('Writing file '+fileh5+'...')
            # coordinates
            Pre = openh5.create_dataset("Pressure",data=Pref,compression='gzip',compression_opts=comp)
            Lat = openh5.create_dataset("Latitude",data=lat_range,compression='gzip',compression_opts=comp)
            Lon = openh5.create_dataset("Longitude",data=lon_range,compression='gzip',compression_opts=comp)

            # data
            Temp = openh5.create_dataset("Temperature",data=Temp_llp,compression='gzip',compression_opts=comp)
            Rho = openh5.create_dataset("Rho",data=Rho_llp,compression='gzip',compression_opts=comp)
            U = openh5.create_dataset("U",data=U_llp,compression='gzip',compression_opts=comp)
            V = openh5.create_dataset("V",data=V_llp,compression='gzip',compression_opts=comp)
            W = openh5.create_dataset("W",data=W_llp,compression='gzip',compression_opts=comp)

            if RT == 1:
                tau_sw = openh5.create_dataset("tau_sw",data=tau_sw_llp,compression='gzip',compression_opts=comp)
                tau_lw = openh5.create_dataset("tau_lw",data=tau_lw_llp,compression='gzip',compression_opts=comp)
                fnet_up = openh5.create_dataset("fnet_up",data=fnet_up_llp,compression='gzip',compression_opts=comp)
                fnet_dn = openh5.create_dataset("fnet_dn",data=fnet_dn_llp,compression='gzip',compression_opts=comp)

            openh5.close()

def KE_spect(input,output,rg):
    tsp = output.nts-output.ntsi+1
    lon, lat = np.meshgrid(rg.lon[:,0],rg.lat[:,0])
    npre = np.shape(rg.Pressure)[0]

    U = rg.U
    V = rg.V
    W = rg.W

    lmax = np.int(np.shape(lat)[0]/2)
    u_coeffs = np.zeros((2,lmax,lmax,npre,tsp),dtype=complex)
    v_coeffs = np.zeros((2,lmax,lmax,npre,tsp),dtype=complex)
    w_coeffs = np.zeros((2,lmax,lmax,npre,tsp),dtype=complex)
    KE_coeffs = np.zeros((2,lmax,lmax,npre,tsp),dtype=complex)
    KE_power = np.zeros((lmax,npre,tsp))

    waven = np.arange(lmax)  #total spherical wavenumber

    fig, ax = plt.subplots(1, 1)
    for t in np.arange(tsp):
        for p in np.arange(npre):
            # sampling = 2 for a lat-lon grid (lon axis is 2x lat axis)
            u_coeffs[:,:,:,p,t] = chairs.expand.SHExpandDHC(U[:,:,p,t],sampling=2)
            v_coeffs[:,:,:,p,t] = chairs.expand.SHExpandDHC(V[:,:,p,t],sampling=2)
            w_coeffs[:,:,:,p,t] = chairs.expand.SHExpandDHC(W[:,:,p,t],sampling=2)

            KE_coeffs[:,:,:,p,t] = (np.abs(u_coeffs[:,:,:,p,t])**2+np.abs(v_coeffs[:,:,:,p,t])**2)/4.0
            KE_power[:,p,t] = chairs.spectralanalysis.spectrum(KE_coeffs,unit='per_lm')

            ax.plot(waven, KE_power[:,p,t])

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set(ylabel='KE (m$^2$ s$^{-2}$)',xlabel='n')


    plt.figure()
    plt.subplot(2,1,1)
    C = plt.contourf(rg.lon[:,0],rg.lat[:,0],rg.Temp[:,:,0,0],cmap='viridis')
    for cc in C.collections:
        cc.set_edgecolor("face")
    clb = plt.colorbar(C)

    plt.show()

    import pdb; pdb.set_trace()

def temperature(input,grid,output,rg,sigmaref):
    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)

    # Set the latitude-longitude grid.
    loni, lati = np.meshgrid(rg.lon,rg.lat)
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    # Averaging in time and longitude.
    if tsp > 1:
        Temperaturel = np.mean(rg.Temperature[:,:,:,:],axis=1)
        Temperaturelt = np.mean(Temperaturel[:,:,:],axis=2)
        del Temperaturel
    else:
        Temperaturelt = np.mean(rg.Temperature[:,:,:,0],axis=1)

    #################
    # Create figure #
    #################

    # Latitude
    latp = rg.lat[:,0]*np.pi/180

    # need to set desired pressure range
    prange = np.where(np.logical_and(rg.Pressure>=np.min(Pref),rg.Pressure<=np.max(Pref)))

    # Contour plot
    C = plt.contourf(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,Temperaturelt[:,prange[0]].T,40)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.gca().invert_yaxis()
    if np.max(Pref)/np.min(Pref) > 100:
        plt.gca().set_yscale("log")
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (bar)')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    plt.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
    clb = plt.colorbar(C)
    clb.set_label(r'Temperature (K)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/temperature_ver_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def u(input,grid,output,rg,sigmaref):
    # contour spacing
    csp = 500

    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)

    # Set the latitude-longitude grid
    loni, lati = np.meshgrid(rg.lon,rg.lat)
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    ##################
    #    Averages    #
    ##################

    # Averaging in time and longitude
    if tsp > 1:
        ZonalMl = np.mean(rg.U[:,:,:,:],axis=1)
        # Vl = np.mean(rg.V[:,:,:,:],axis=1)
        # Wl = np.mean(rg.W[:,:,:,:],axis=1)
        ZonalMlt = np.mean(ZonalMl[:,:,:],axis=2)
        # Vlt = np.mean(Vl[:,:,:],axis=2)
        # Wlt = np.mean(Wl[:,:,:],axis=2)
        del ZonalMl, Vl, Wl
    else:
        ZonalMlt = np.mean(rg.U[:,:,:,0],axis=1)
        # Vlt = np.mean(rg.V[:,:,:,0],axis=1)
        # Wlt = np.mean(rg.W[:,:,:,0],axis=1)

    #################
    # Create figure #
    #################

    # set up arrows
    # vspacing = 18
    # wspacing = 1
    # Vq = Vlt[::vspacing,::wspacing].ravel()
    # Wq = Wlt[::vspacing,::wspacing].ravel()
    # #preq = rg.Pressure[:,0][::spacing,::spacing].ravel()
    # #latq = lati[::spacing,::spacing].ravel()
    # latq, preq = np.meshgrid(rg.lat[::vspacing,0],rg.Pressure[::wspacing,0])

    # Latitude
    latp = rg.lat[:,0]*np.pi/180

    # need to set desired pressure range (major PITA!)
    prange = np.where(np.logical_and(rg.Pressure>=np.min(Pref),rg.Pressure<=np.max(Pref)))

    # Contour plot
    C = plt.contourf(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,ZonalMlt[:,prange[0]].T,40,cmap='viridis')
    clb = plt.colorbar(C,extend='both')
    clb.set_label(r'Velocity (m s$^{-1}$)')
    levp = np.arange(np.ceil(np.min(ZonalMlt)/csp)*csp,np.floor(np.max(ZonalMlt)/csp)*csp,csp)
    c2 = plt.contour(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,ZonalMlt[:,prange[0]].T,levels=levp,colors='w',linewidths=1)
    plt.clabel(c2,inline=1,fontsize=10)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.gca().invert_yaxis()
    # plt.quiver(latq.ravel(),preq.ravel()/1e5,Vq/np.max(Vq),Wq/np.max(Wq),color='0.5')
    if np.max(Pref)/np.min(Pref) > 100:
        plt.gca().set_yscale("log")
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (bar)')
    plt.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/u_ver_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def w_prof(input,grid,output):
    tsp = output.nts-output.ntsi+1

    w = 0.5*(output.Wh[:,1:,:]+output.Wh[:,:-1,:])/output.Rho

    lon_shift = grid.lon*1.0
    lon_shift[lon_shift>7*np.pi/4] -= 2*np.pi

    day = np.where(np.logical_and(lon_shift>=-np.pi/4,lon_shift<np.pi/4))[0]
    eve = np.where(np.logical_and(lon_shift>=np.pi/4,lon_shift<3*np.pi/4))[0]
    night = np.where(np.logical_and(lon_shift>=3*np.pi/4,lon_shift<5*np.pi/4))[0]
    morn = np.where(np.logical_and(lon_shift>=5*np.pi/4,lon_shift<7*np.pi/4))[0]

    pday = output.Pressure[day,:,:]
    peve = output.Pressure[eve,:,:]
    pnight = output.Pressure[night,:,:]
    pmorn = output.Pressure[morn,:,:]

    pday_avg = np.mean(np.mean(pday,axis=2),axis=0)
    peve_avg = np.mean(np.mean(peve,axis=2),axis=0)
    pnight_avg = np.mean(np.mean(pnight,axis=2),axis=0)
    pmorn_avg = np.mean(np.mean(pmorn,axis=2),axis=0)

    wday_tmp = w[day,:,:]
    weve_tmp = w[eve,:,:]
    wnight_tmp = w[night,:,:]
    wmorn_tmp = w[morn,:,:]

    wday = np.zeros_like(wday_tmp)
    weve = np.zeros_like(weve_tmp)
    wnight = np.zeros_like(wnight_tmp)
    wmorn = np.zeros_like(wmorn_tmp)

    for t in np.arange(tsp):
        for i in np.arange(np.max([len(day),len(eve),len(night),len(morn)])):
            if i < len(day):
                 wday[i,:,t] = interp.pchip_interpolate(pday[i,::-1,t],wday_tmp[i,::-1,t],pday_avg[::-1])[::-1]
            if i < len(eve):
                 weve[i,:,t] = interp.pchip_interpolate(peve[i,::-1,t],weve_tmp[i,::-1,t],peve_avg[::-1])[::-1]
            if i < len(night):
                 wnight[i,:,t] = interp.pchip_interpolate(pnight[i,::-1,t],wnight_tmp[i,::-1,t],pnight_avg[::-1])[::-1]
            if i < len(morn):
                 wmorn[i,:,t] = interp.pchip_interpolate(pmorn[i,::-1,t],wmorn_tmp[i,::-1,t],pmorn_avg[::-1])[::-1]

    wday_avg = np.mean(np.mean(wday,axis=2),axis=0)
    weve_avg = np.mean(np.mean(weve,axis=2),axis=0)
    wnight_avg = np.mean(np.mean(wnight,axis=2),axis=0)
    wmorn_avg = np.mean(np.mean(wmorn,axis=2),axis=0)

    plt.semilogy(wday_avg,pday_avg/1e5,'r-',label='day')
    plt.plot(weve_avg,peve_avg/1e5,'g-',label='evening')
    plt.plot(wnight_avg,pnight_avg/1e5,'b-',label='night')
    plt.plot(wmorn_avg,pmorn_avg/1e5,'c-',label='morning')
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure (bar)')
    plt.xlabel(r'Vertical wind speed [m s$^{-1}$]')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    plt.legend(loc='lower right',fontsize = 8)
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.savefig(input.resultsf+'/figures/Wprofile_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def w_ver(input,grid,output,rg,sigmaref):
    # contour spacing
    csp = 200

    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)

    # Set the latitude-longitude grid
    loni, lati = np.meshgrid(rg.lon,rg.lat)
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    # Averaging in time and longitude
    if tsp > 1:
        VertMl = np.mean(rg.W[:,:,:,:],axis=1)
        VertMeq = np.mean(rg.W[np.logical_and(lati>=-20.0,lati<=20.0)[:,1],:,:,:],axis=0)
        VertMlt = np.mean(VertMl[:,:,:],axis=2)
        VertMeqt = np.mean(VertMeq[:,:,:],axis=2)
        del VertMl, VertMeq
    else:
        VertMlt = np.mean(rg.W[:,:,:,0],axis=1)
        VertMeqt = np.mean(rg.W[np.logical_and(lati>=-20.0,lati<=20.0)[:,1],:,:,0],axis=0)

    #################
    # Create figure #
    #################

    # Latitude
    latp = rg.lat[:,0]*np.pi/180
    lonp = rg.lon[:,0]*np.pi/180

    # need to set desired pressure range
    prange = np.where(np.logical_and(rg.Pressure>=np.min(Pref),rg.Pressure<=np.max(Pref)))

    # Contour plot
    C = plt.contourf(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,VertMlt[:,prange[0]].T,40,cmap='viridis')
    levp = np.arange(np.ceil(np.min(VertMlt)/csp)*csp,np.floor(np.max(VertMlt)/csp)*csp,csp)
    c2 = plt.contour(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,VertMlt[:,prange[0]].T,[0],colors='w',linewidths=1)
    # plt.clabel(c2,inline=1,fontsize=10)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.gca().invert_yaxis()
    if np.max(Pref)/np.min(Pref) > 100:
        plt.gca().set_yscale("log")
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (bar)')
    plt.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Velocity (m s$^{-1}$)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/w_ver_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

    plt.figure()
    C = plt.contourf(lonp*180/np.pi,rg.Pressure[prange[0],0]/1e5,VertMeqt[:,prange[0]].T,40,cmap='viridis')
    levp = np.arange(np.ceil(np.min(VertMeqt)/csp)*csp,np.floor(np.max(VertMeqt)/csp)*csp,csp)
    c2 = plt.contour(lonp*180/np.pi,rg.Pressure[prange[0],0]/1e5,VertMeqt[:,prange[0]].T,[0],colors='w',linewidths=1)
    # plt.clabel(c2,inline=1,fontsize=10)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.gca().invert_yaxis()
    if np.max(Pref)/np.min(Pref) > 100:
        plt.gca().set_yscale("log")
    plt.xlabel('Longitude (deg)')
    plt.ylabel('Pressure (bar)')
    plt.plot(lonp*180/np.pi,np.zeros_like(lonp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Velocity (m s$^{-1}$)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/weq_ver_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def uv_lev(input,grid,output,rg,Plev):
    loni, lati = np.meshgrid(rg.lon[:,0],rg.lat[:,0])

    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    Uii = np.zeros(np.shape(loni)+(tsp,))
    Vii = np.zeros(np.shape(loni)+(tsp,))

    for t in np.arange(tsp):
        # interpolate
        # above is index of higher pressure bound (lower index)
        # below is index of lower pressure bound (higher index)
        if np.size(rg.Pressure[rg.Pressure>Plev]) == 0:
            above = np.where(rg.Pressure==np.max(rg.Pressure))[0][0]
            below = above + 1
        elif np.size(rg.Pressure[rg.Pressure<Plev]) == 0:
            below = np.where(rg.Pressure==np.min(rg.Pressure))[0][0]
            above = below - 1
        else:
            above = np.where(rg.Pressure==np.min(rg.Pressure[rg.Pressure>Plev]))[0][0]
            below = np.where(rg.Pressure==np.max(rg.Pressure[rg.Pressure<Plev]))[0][0]

        Uii[:,:,t] = (rg.U[:,:,below,t]*(rg.Pressure[above,t] - Plev)\
            + rg.U[:,:,above,t]*(Plev - rg.Pressure[below,t]))\
            / (rg.Pressure[above,t]-rg.Pressure[below,t])
        Vii[:,:,t] = (rg.V[:,:,below,t]*(rg.Pressure[above,t] - Plev)\
            + rg.V[:,:,above,t]*(Plev - rg.Pressure[below,t]))\
            / (rg.Pressure[above,t]-rg.Pressure[below,t])
        # Wii[:,:,t] = (rg.W[:,:,below,t]*(rg.Pressure[above,t] - Plev)\
        #     + rg.W[:,:,above,t]*(Plev - rg.Pressure[below,t]))\
        #     / (rg.Pressure[above,t]-rg.Pressure[below,t])

    # Averaging in time
    if tsp > 1:
        Uiii = np.mean(Uii,axis=2)
        Viii = np.mean(Vii,axis=2)
        del Uii, Vii
    else:
        Uiii = Uii[:,:,0]
        Viii = Vii[:,:,0]

    #################
    # Create Figure #
    #################

    plt.figure()
    lonp = rg.lon[:,0]
    latp = rg.lat[:,0]
    C = plt.contourf(lonp,latp,Uiii,50,cmap='viridis')
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Velocity (m s$^{-1}$)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/u_lev%#.3fmbar_i%d_l%d.pdf'%((Plev/100),output.ntsi,output.nts))
    plt.close()

    plt.figure()
    C = plt.contourf(lonp,latp,Viii,50,cmap='viridis')
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Velocity (m s$^{-1}$)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/grid.figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/v_lev%#.3fmbar_i%d_l%d.pdf'%((Plev/100),output.ntsi,output.nts))
    plt.close()

def temperature_u_lev(input,grid,output,rg,Plev):
    # Set the latitude-longitude grid.
    loni, lati = np.meshgrid(rg.lon[:,0],rg.lat[:,0])

    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    Temperatureii = np.zeros(np.shape(loni)+(tsp,))
    Uii = np.zeros(np.shape(loni)+(tsp,))
    Vii = np.zeros(np.shape(loni)+(tsp,))

    for t in np.arange(tsp):
        # interpolate
        # above is index of higher pressure bound (lower index)
        # below is index of lower pressure bound (higher index)
        if np.size(rg.Pressure[rg.Pressure>Plev]) == 0:
            above = np.where(rg.Pressure==np.max(rg.Pressure))[0][0]
            below = above + 1
        elif np.size(rg.Pressure[rg.Pressure<Plev]) == 0:
            below = np.where(rg.Pressure==np.min(rg.Pressure))[0][0]
            above = below - 1
        else:
            above = np.where(rg.Pressure==np.min(rg.Pressure[rg.Pressure>Plev]))[0][0]
            below = np.where(rg.Pressure==np.max(rg.Pressure[rg.Pressure<Plev]))[0][0]

        Temperatureii[:,:,t] = (rg.Temperature[:,:,below,t]*(rg.Pressure[above,t] - Plev)\
            + rg.Temperature[:,:,above,t]*(Plev - rg.Pressure[below,t]))\
            / (rg.Pressure[above,t]-rg.Pressure[below,t])
        Uii[:,:,t] = (rg.U[:,:,below,t]*(rg.Pressure[above,t] - Plev)\
            + rg.U[:,:,above,t]*(Plev - rg.Pressure[below,t]))\
            / (rg.Pressure[above,t]-rg.Pressure[below,t])
        Vii[:,:,t] = (rg.V[:,:,below,t]*(rg.Pressure[above,t] - Plev)\
            + rg.V[:,:,above,t]*(Plev - rg.Pressure[below,t]))\
            / (rg.Pressure[above,t]-rg.Pressure[below,t])

    # Averaging in time
    if tsp > 1:
        Temperature = np.mean(Temperatureii,axis=2)
        Uiii = np.mean(Uii,axis=2)
        Viii = np.mean(Vii,axis=2)
        del Temperatureii, Uii, Vii
    else:
        Temperature = Temperatureii[:,:,0]
        Uiii = Uii[:,:,0]
        Viii = Vii[:,:,0]

    # Wind arrays
    d_z = np.shape(Uiii)
    spacing = 40
    U = Uiii[::spacing,::spacing].ravel()
    V = Viii[::spacing,::spacing].ravel()
    lonq = loni[::spacing,::spacing].ravel()
    latq = lati[::spacing,::spacing].ravel()
    del Uiii, Viii

    #################
    # Create Figure #
    #################

    lonp = rg.lon[:,0]
    latp = rg.lat[:,0]
    C = plt.contourf(lonp,latp,Temperature,50)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')
    plt.quiver(lonq,latq,U,V,color='0.5')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Temperature (K)')

    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/temperature-uv_lev%#.3fmbar_i%d_l%d.pdf'%(Plev/100,output.ntsi,output.nts))
    plt.close()

def tracer_u_lev(input,grid,output,Plev,trace):
    # Set the latitude-longitude grid.
    res_deg = 0.005
    loni, lati = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2,np.pi/2,res_deg))
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    #######################
    # Winds & Tracer      #
    #######################

    # id tracer
    tracer = getattr(output,trace)

    # Initialize arrays
    tr = np.zeros((grid.nv,1))
    Pr = np.zeros((grid.nv,1))
    Mx = np.zeros((grid.nv,1))
    My = np.zeros((grid.nv,1))
    Mz = np.zeros((grid.nv,1))
    tri = np.zeros((grid.point_num,1))
    Rhot = np.zeros((grid.point_num,1))
    Mxf = np.zeros((grid.point_num,1))
    Myf = np.zeros((grid.point_num,1))
    Mzf = np.zeros((grid.point_num,1))
    Ui = np.zeros((grid.point_num,1))
    Vi = np.zeros((grid.point_num,1))
    trii = np.zeros((d_lon[0],d_lon[1],tsp))
    Uii = np.zeros((d_lon[0],d_lon[1],tsp))
    Vii = np.zeros((d_lon[0],d_lon[1],tsp))

    # Compute winds and temperatures
    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            for lev in np.arange(grid.nv):
                Pr[lev] = output.Pressure[i,lev,t]
                tr[lev] = tracer[i,lev,t]
                Mx[lev] = output.Mh[0,i,lev,t]
                My[lev] = output.Mh[1,i,lev,t]
                Mz[lev] = output.Mh[2,i,lev,t]
            # Interpolate in pressure
            Rhot[i] = interp.interp1d(Pr.T[0],output.Rho[i,:,t],kind='linear',fill_value='extrapolate')(Plev)
            Mxf[i] = interp.interp1d(Pr.T[0],Mx.T[0],kind='linear',fill_value='extrapolate')(Plev)
            Myf[i] = interp.interp1d(Pr.T[0],My.T[0],kind='linear',fill_value='extrapolate')(Plev)
            Mzf[i] = interp.interp1d(Pr.T[0],Mz.T[0],kind='linear',fill_value='extrapolate')(Plev)
            tri[i] = interp.interp1d(Pr.T[0],tr.T[0],kind='linear',fill_value='extrapolate')(Plev)
        for i in np.arange(grid.point_num):
            Ui[i] = (Mxf[i]*(-np.sin(grid.lon[i])) + \
                     Myf[i]*np.cos(grid.lon[i]) + \
                     Mzf[i]*(0))/Rhot[i]
            Vi[i] = (Mxf[i]*(-np.sin(grid.lat[i])*np.cos(grid.lon[i])) + \
                     Myf[i]*(-np.sin(grid.lat[i])*np.sin(grid.lon[i])) + \
                     Mzf[i]*np.cos(grid.lat[i]))/Rhot[i]
        # Convert icosahedral grid into lon-lat grid
        # import pdb; pdb.set_trace()
        trii[:,:,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,tri.T[0],(loni,lati),method='linear')
        Uii[:,:,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Ui.T[0],(loni,lati),method='cubic')
        Vii[:,:,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Vi.T[0],(loni,lati),method='cubic')

    # Averaging in time
    if tsp > 1:
        tr_plot = np.mean(trii,axis=2)
        Uiii = np.mean(Uii,axis=2)
        Viii = np.mean(Vii,axis=2)
        del trii, Uii, Vii
    else:
        tr_plot = trii[:,:,0]
        Uiii = Uii[:,:,0]
        Viii = Vii[:,:,0]

    # Wind arrays
    d_z = np.shape(Uiii)
    spacing = 80
    U = Uiii[::spacing,::spacing].ravel()
    V = Viii[::spacing,::spacing].ravel()
    lonq = loni[::spacing,::spacing].ravel()
    latq = lati[::spacing,::spacing].ravel()
    del Uiii, Viii

    #################
    # Create Figure #
    #################

    lonp = np.arange(0,2*np.pi,res_deg)
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)
    C = plt.contourf(lonp*180/np.pi,latp*180/np.pi,np.log10(tr_plot),50)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')
    plt.quiver(lonq*180/np.pi,latq*180/np.pi,U,V,color='0.5')
    #plt.plot(grid.lon*180/np.pi,grid.lat*180/np.pi,'w.')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Log(mixing ratio)')

    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/chem-'+trace+'-uv_lev%#.3fmbar_i%d_l%d.pdf'%(Plev/100,output.ntsi,output.nts))
    plt.close()

def potential_temp(input,grid,output,rg,sigmaref):
    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient

    # Set the latitude-longitude grid.
    res_deg = 0.005
    loni, lati = np.meshgrid(rg.lon,rg.lat)
    tsp = output.nts-output.ntsi+1

    #########################
    # Potential temperature #
    #########################

    Theta = rg.Pressure[None,None,:,:]/(input.Rd*rg.Rho) * \
             (rg.Pressure[None,None,:,:]/input.P_Ref)**(-kappa_ad)

    # Averaging in time and longitude.
    if tsp > 1:
        Thetal = np.mean(Theta[:,:,:,:],axis=1)
        del Theta
        Thetalt = np.mean(Thetal[:,:,:],axis=2)
        del Thetal
    else:
        Thetalt = np.mean(Theta[:,:,:,0],axis=1)
        del Theta

    #################
    # Create figure #
    #################

    # Latitude
    latp = rg.lat[:,0]*np.pi/180

    # need to set desired pressure range
    prange = np.where(np.logical_and(rg.Pressure>=np.min(Pref),rg.Pressure<=np.max(Pref)))

    # Contour plot
    C = plt.contourf(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,Thetalt[:,prange[0]].T,40,linewidths=None)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.gca().invert_yaxis()
    if np.max(Pref)/np.min(Pref) > 100:
        plt.gca().set_yscale("log")
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (bar)')
    plt.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')

    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Potential temperature (K)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/potential_temp_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def dFunc_dZ(f, alti, nv):
    dfdz = np.zeros_like(f)
    dz = alti[0,0,1,0]-alti[0,0,0,0]

    # Over z (vertical direction)
    for t in np.arange(np.shape(f)[-1]):
        for i in np.arange(nv):
            if i == 0:
                #lower edge, first order gradient
                dfdz[:,:,i,t] = (f[:,:,i+1,t]-f[:,:,i,t])/dz
            elif i == 1 or i == (nv-2):
                #second order gradient at penultimate layers
                dfdz[:,:,i,t] = (f[:,:,i+1,t]-f[:,:,i-1,t])/(2*dz)
            elif i == (nv-1):
                #highest edge, first order gradient
                dfdz[:,:,i,t] = (f[:,:,i,t]-f[:,:,i-1,t])/dz
            else:
                #five point stencil
                dfdz[:,:,i,t] = (-f[:,:,i+2,t]+8*f[:,:,i+1,t]\
                                    -8*f[:,:,i-1,t]+f[:,:,i-2,t])/(12*dz)
    return dfdz

def dFunc_dLat(f, dLat):
    dfdlat = np.zeros_like(f)

    # Over lat (vertical direction)
    for t in np.arange(np.shape(f)[-1]):
        for i in np.arange(np.shape(f)[0]):
            if i == 0:
                #lower edge, first order gradient
                dfdlat[i,:,:,t] = (f[i+1,:,:,t]-f[i,:,:,t])/dLat
            elif i == 1 or i == (np.shape(f)[0]-2):
                #second order gradient at penultimate layers
                dfdlat[i,:,:,t] = (f[i+1,:,:,t]-f[i-1,:,:,t])/(2*dLat)
            elif i == (np.shape(f)[0]-1):
                #highest edge, first order gradient
                dfdlat[i,:,:,t] = (f[i,:,:,t]-f[i-1,:,:,t])/dLat
            else:
                #five point stencil
                dfdlat[i,:,:,t] = (-f[i+2,:,:,t]+8*f[i+1,:,:,t]\
                                    -8*f[i-1,:,:,t]+f[i-2,:,:,t])/(12*dLat)
    return dfdlat

def dFunc_dLon(f, dLon):
    dfdlon = np.zeros_like(f)

    # Over lon (vertical direction)
    for t in np.arange(np.shape(f)[-1]):
        for i in np.arange(np.shape(f)[1]):
            if i == 0:
                #lower edge, first order gradient
                dfdlon[:,i,:,t] = (f[:,i+1,:,t]-f[:,i,:,t])/dLon
            elif i == 1 or i == (np.shape(f)[1]-2):
                #second order gradient at penultimate layers
                dfdlon[:,i,:,t] = (f[:,i+1,:,t]-f[:,i-1,:,t])/(2*dLon)
            elif i == (np.shape(f)[1]-1):
                #highest edge, first order gradient
                dfdlon[:,i,:,t] = (f[:,i,:,t]-f[:,i-1,:,t])/dLon
            else:
                #five point stencil
                dfdlon[:,i,:,t] = (-f[:,i+2,:,t]+8*f[:,i+1,:,t]\
                                    -8*f[:,i-1,:,t]+f[:,i-2,:,t])/(12*dLon)
    return dfdlon

def CurlF(fr,flat,flon,lati,alti,res_deg,nv,A):
    curlFz = np.zeros_like(fr)
    curlFlat = np.zeros_like(fr)
    curlFlon = np.zeros_like(fr)

    curlFz = -1.0*(dFunc_dLat(flon*np.cos(lati),res_deg)-dFunc_dLon(flat,res_deg))
    curlFz /= (np.cos(lati)*(A+alti))

    curlFlat = -1.0*(dFunc_dLon(fr,res_deg)/np.cos(lati)-dFunc_dZ((A+alti)*flon,alti,nv))
    curlFlat /= (A+alti)

    curlFlon = -1.0*(dFunc_dZ((A+alti)*flat,alti,nv)-dFunc_dLat(fr,res_deg))
    curlFlon /= (A+alti)

    return curlFz, curlFlat, curlFlon

def potential_vort_vert(input,grid,output,sigmaref):
    # Set the reference pressure

    Pref = input.P_Ref*sigmaref
    #Pref = np.array([sigmaref])
    d_sig = np.size(sigmaref)
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient

    # Set the latitude-longitude grid.
    res_deg = 0.1
    tsp = output.nts-output.ntsi+1
    lonm, altm = np.meshgrid(grid.lon,grid.Altitude)
    latm, altm = np.meshgrid(grid.lat,grid.Altitude)
    loni, lati, alti, ti = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2+2*res_deg,np.pi/2-2*res_deg,res_deg),grid.Altitude,np.arange(tsp))
    d_lon = np.shape(loni)

    #########################
    # Potential temperature #
    #########################

    # Initialize arrays
    Thetai = np.zeros((grid.point_num,grid.nv))
    Theta = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Ui = np.zeros((grid.point_num,grid.nv))
    U = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Vi = np.zeros((grid.point_num,grid.nv))
    V = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Wi = np.zeros((grid.point_num,grid.nv))
    W = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Rholl = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Pressll = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))

    # Compute temperatures
    for t in np.arange(tsp):
        Thetai = output.Pressure[:,:,t]/(input.Rd*output.Rho[:,:,t]) * \
                    (output.Pressure[:,:,t]/input.P_Ref)**(-kappa_ad)
        Ui = (output.Mh[0,:,:,t]*(-np.sin(lonm.T)) + \
                     output.Mh[1,:,:,t]*np.cos(lonm.T) + \
                     output.Mh[2,:,:,t]*(0))/output.Rho[:,:,t]
        Vi = (output.Mh[0,:,:,t]*(-np.sin(latm.T)*np.cos(lonm.T)) + \
                     output.Mh[1,:,:,t]*(-np.sin(latm.T)*np.sin(lonm.T)) + \
                     output.Mh[2,:,:,t]*np.cos(latm.T))/output.Rho[:,:,t]
        Wi = (output.Mh[0,:,:,t]*(np.cos(latm.T)*np.cos(lonm.T)) + \
                     output.Mh[1,:,:,t]*(np.cos(latm.T)*np.sin(lonm.T)) + \
                     output.Mh[2,:,:,t]*np.sin(latm.T))/output.Rho[:,:,t]
        for lev in np.arange(grid.nv):
            Theta[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Thetai[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            U[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Ui[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            V[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Vi[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            W[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Wi[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            Rholl[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,output.Rho[:,lev,t],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            Pressll[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,output.Pressure[:,lev,t],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')

    # Compute gradient of potential temp
    dThdz = dFunc_dZ(Theta,alti,grid.nv)

    dThdlat = dFunc_dLat(Theta,res_deg)/(input.A+alti)

    dThdlon = dFunc_dLon(Theta,res_deg)/(input.A+alti)/np.cos(lati)

    # Compute curl of wind speeds
    curlVz, curlVlat, curlVlon = CurlF(W,V,U,lati,alti,res_deg,grid.nv,input.A)

    # ok, no idea if I got the curl right, but now I need to add 2*Omega and dot into gradTheta
    curlVlon += 2*input.Omega

    curlVz /= Rholl
    curlVlat /= Rholl
    curlVlon /= Rholl

    Phi_z = curlVz * dThdz + curlVlat * dThdlat + curlVlon * dThdlon

    # Now, I need to interpolate to a pressure grid
    Phi = np.zeros((d_lon[0],d_lon[1], d_sig, tsp))
    for t in np.arange(tsp):
        for ilat in np.arange(d_lon[0]):
            for ilon in np.arange(d_lon[1]):
                Phi[ilat,ilon,:,t] = interp.pchip_interpolate(Pressll[ilat,ilon,:,t][::-1],
                                            Phi_z[ilat,ilon,:,t][::-1],Pref[::-1])[::-1]
    # Averaging in time and longitude.
    if tsp > 1:
        Phil = np.mean(Phi,axis=1)
        Philt = np.mean(Phil,axis=2)
        del Phil
    else:
        Philt = np.mean(Phi[:,:,:,0],axis=1)

    #################
    # Create figure #
    #################

    # Latitude
    latp = lati[:,0,0,0]
    lonp = loni[0,:,0,0]

    # Contour plot
    C = plt.contourf(latp*180/np.pi,Pref/1e5,Philt.T,40,linewidths=None,cmap='plasma')
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0

    plt.gca().invert_yaxis()
    if np.max(Pref)/np.min(Pref) > 100:
        plt.gca().set_yscale("log")
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (bar)')
    plt.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')

    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    #plt.ylabel('Altitude (m)')
    clb = plt.colorbar(C)
    clb.set_label(r'Potential vorticity (K m$^2$ kg$^{-1}$ s$^{-1}$)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/pot_vort_z_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def potential_vort_lev(input,grid,output,sigmaref):
    # Set the reference pressure

    #Pref = ps0*sigmaref
    Pref = np.array([sigmaref])
    d_sig = np.size(sigmaref)
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient

    # Set the latitude-longitude grid.
    res_deg = 0.1
    tsp = output.nts-output.ntsi+1
    lonm, altm = np.meshgrid(grid.lon,grid.Altitude)
    latm, altm = np.meshgrid(grid.lat,grid.Altitude)
    loni, lati, alti, ti = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2+2*res_deg,np.pi/2-2*res_deg,res_deg),grid.Altitude,np.arange(tsp))
    d_lon = np.shape(loni)

    #########################
    # Potential temperature #
    #########################

    # Initialize arrays
    Thetai = np.zeros((grid.point_num,grid.nv))
    Theta = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Ui = np.zeros((grid.point_num,grid.nv))
    U = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Vi = np.zeros((grid.point_num,grid.nv))
    V = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Wi = np.zeros((grid.point_num,grid.nv))
    W = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Rholl = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Pressll = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))

    # Compute temperatures
    for t in np.arange(tsp):
        Thetai = output.Pressure[:,:,t]/(input.Rd*output.Rho[:,:,t]) * \
                    (output.Pressure[:,:,t]/input.P_Ref)**(-kappa_ad)
        Ui = (output.Mh[0,:,:,t]*(-np.sin(lonm.T)) + \
                     output.Mh[1,:,:,t]*np.cos(lonm.T) + \
                     output.Mh[2,:,:,t]*(0))/output.Rho[:,:,t]
        Vi = (output.Mh[0,:,:,t]*(-np.sin(latm.T)*np.cos(lonm.T)) + \
                     output.Mh[1,:,:,t]*(-np.sin(latm.T)*np.sin(lonm.T)) + \
                     output.Mh[2,:,:,t]*np.cos(latm.T))/output.Rho[:,:,t]
        Wi = (output.Mh[0,:,:,t]*(np.cos(latm.T)*np.cos(lonm.T)) + \
                     output.Mh[1,:,:,t]*(np.cos(latm.T)*np.sin(lonm.T)) + \
                     output.Mh[2,:,:,t]*np.sin(latm.T))/output.Rho[:,:,t]
        for lev in np.arange(grid.nv):
            Theta[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Thetai[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            U[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Ui[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            V[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Vi[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            W[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Wi[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            Rholl[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,output.Rho[:,lev,t],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            Pressll[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,output.Pressure[:,lev,t],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')

    # Compute gradient of potential temp
    dThdz = dFunc_dZ(Theta,alti,grid.nv)

    dThdlat = dFunc_dLat(Theta,res_deg)/(input.A+alti)

    dThdlon = dFunc_dLon(Theta,res_deg)/(input.A+alti)/np.cos(lati)

    # Compute curl of wind speeds
    curlVz, curlVlat, curlVlon = CurlF(W,V,U,lati,alti,res_deg,grid.nv,input.A)

    # ok, no idea if I got the curl right, but now I need to add 2*Omega and dot into gradTheta
    curlVlon += 2*input.Omega

    curlVz /= Rholl
    curlVlat /= Rholl
    curlVlon /= Rholl

    Phi_z = curlVz * dThdz + curlVlat * dThdlat + curlVlon * dThdlon

    # Now, I need to interpolate to a pressure grid (or should I?)
    Phi = np.zeros((d_lon[0],d_lon[1], d_sig, tsp))
    for t in np.arange(tsp):
        for ilat in np.arange(d_lon[0]):
            for ilon in np.arange(d_lon[1]):
                Phi[ilat,ilon,:,t] = interp.pchip_interpolate(Pressll[ilat,ilon,:,t][::-1],
                                            Phi_z[ilat,ilon,:,t][::-1],Pref[::-1])[::-1]
    #################
    # Create figure #
    #################

    # Latitude
    latp = lati[:,0,0,0]
    lonp = loni[0,:,0,0]

    # Contour plot
    C = plt.contourf(lonp*180/np.pi,latp*180/np.pi,Phi[:,:,0,0],40,linewidths=None,cmap='plasma')
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    # plt.gca().invert_yaxis()
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Potential vorticity (K m$^2$ kg$^{-1}$ s$^{-1}$)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/pot_vort_lev%d_i%d_l%d.pdf'%(sigmaref/input.P_Ref*1000,output.ntsi,output.nts))
    plt.close()

def rela_vort_lev(input,grid,output,sigmaref):
    # Set the reference pressure

    #Pref = ps0*sigmaref
    Pref = np.array([sigmaref])
    d_sig = np.size(sigmaref)
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient

    # Set the latitude-longitude grid.
    res_deg = 0.1
    tsp = output.nts-output.ntsi+1
    lonm, altm = np.meshgrid(grid.lon,grid.Altitude)
    latm, altm = np.meshgrid(grid.lat,grid.Altitude)
    loni, lati, alti, ti = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2+2*res_deg,np.pi/2-2*res_deg,res_deg),grid.Altitude,np.arange(tsp))
    d_lon = np.shape(loni)

    #########################
    # Potential temperature #
    #########################

    # Initialize arrays
    # Thetai = np.zeros((grid.point_num,grid.nv))
    # Theta = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Ui = np.zeros((grid.point_num,grid.nv))
    U = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Vi = np.zeros((grid.point_num,grid.nv))
    V = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Wi = np.zeros((grid.point_num,grid.nv))
    W = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Rholl = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))
    Pressll = np.zeros((d_lon[0],d_lon[1],grid.nv,tsp))

    # Compute temperatures
    for t in np.arange(tsp):
        # Thetai = output.Pressure[:,:,t]/(input.Rd*output.Rho[:,:,t]) * \
                    # (output.Pressure[:,:,t]/input.P_Ref)**(-kappa_ad)
        Ui = (output.Mh[0,:,:,t]*(-np.sin(lonm.T)) + \
                     output.Mh[1,:,:,t]*np.cos(lonm.T) + \
                     output.Mh[2,:,:,t]*(0))/output.Rho[:,:,t]
        Vi = (output.Mh[0,:,:,t]*(-np.sin(latm.T)*np.cos(lonm.T)) + \
                     output.Mh[1,:,:,t]*(-np.sin(latm.T)*np.sin(lonm.T)) + \
                     output.Mh[2,:,:,t]*np.cos(latm.T))/output.Rho[:,:,t]
        Wi = (output.Mh[0,:,:,t]*(np.cos(latm.T)*np.cos(lonm.T)) + \
                     output.Mh[1,:,:,t]*(np.cos(latm.T)*np.sin(lonm.T)) + \
                     output.Mh[2,:,:,t]*np.sin(latm.T))/output.Rho[:,:,t]
        for lev in np.arange(grid.nv):
            # Theta[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Thetai[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            U[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Ui[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            V[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Vi[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            W[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Wi[:,lev],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            Rholl[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,output.Rho[:,lev,t],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')
            Pressll[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,output.Pressure[:,lev,t],(loni[:,:,0,0],lati[:,:,0,0]),method='nearest')

    # Compute gradient of potential temp
    # dThdz = dFunc_dZ(Theta,alti,grid.nv)
    #
    # dThdlat = dFunc_dLat(Theta,res_deg)/(input.A+alti)
    #
    # dThdlon = dFunc_dLon(Theta,res_deg)/(input.A+alti)/np.cos(lati)

    # Compute curl of wind speeds
    curlVz, curlVlat, curlVlon = CurlF(W,V,U,lati,alti,res_deg,grid.nv,input.A)

    # ok, no idea if I got the curl right, but now I need to add 2*Omega and dot into gradTheta
    # curlVlon += 2*input.Omega
    #
    # curlVz /= Rholl
    # curlVlat /= Rholl
    # curlVlon /= Rholl
    #
    # Phi_z = curlVz * dThdz + curlVlat * dThdlat + curlVlon * dThdlon

    # Now, I need to interpolate to a pressure grid (or should I?)
    curlVz2 = np.zeros((d_lon[0],d_lon[1], d_sig, tsp))
    for t in np.arange(tsp):
        for ilat in np.arange(d_lon[0]):
            for ilon in np.arange(d_lon[1]):
                curlVz2[ilat,ilon,:,t] = interp.pchip_interpolate(Pressll[ilat,ilon,:,t][::-1],
                                            curlVz[ilat,ilon,:,t][::-1],Pref[::-1])[::-1]
    #################
    # Create figure #
    #################

    # Latitude
    latp = lati[:,0,0,0]
    lonp = loni[0,:,0,0]

    # Contour plot
    C = plt.contourf(lonp*180/np.pi,latp*180/np.pi,curlVz2[:,:,0,0],40,linewidths=None,cmap='viridis')
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    # plt.gca().invert_yaxis()
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Relative vorticity (s$^{-1}$)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/rela_vort_lev%d_i%d_l%d.pdf'%(sigmaref/input.P_Ref*1000,output.ntsi,output.nts))
    plt.close()

def TPprof(input,grid,output,sigmaref,column):
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)

    tsp = output.nts-output.ntsi+1

    for column in np.arange(0,grid.point_num,50):
        if tsp > 1:
            P = np.mean(output.Pressure[column,:,:],axis=2)
            T = np.mean(output.Pressure[column,:,:]/input.Rd/output.Rho[column,:,:],axis=2)
        else:
            P = output.Pressure[column,:,0]
            T = output.Pressure[column,:,0]/input.Rd/output.Rho[column,:,0]

        kappa = input.Rd/input.Cp

        plt.semilogy(T,P/1e5,'k-',alpha= 0.5,lw=1)
        plt.plot(T[np.int(np.floor(grid.nv/2))],P[np.int(np.floor(grid.nv/2))]/100000,'r+',ms =5,alpha=0.5)
        plt.plot(T[np.int(np.floor(grid.nv*0.75))],P[np.int(np.floor(grid.nv*0.75))]/100000,'g+',ms =5,alpha=0.5)

    Tad = T[15]*(P/P[15])**kappa # need to come up with better way of visualizing the adiabatic conditions

    # plt.plot(Tad,P/100,'r--')
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure (bar)')
    plt.xlabel('Temperature [K]')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.savefig(input.resultsf+'/figures/TPprofile_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def PTPprof(input,grid,output,sigmaref,column):
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient

    tsp = output.nts-output.ntsi+1

    for column in np.arange(0,grid.point_num,50):
        if tsp > 1:
            P = np.mean(output.Pressure[column,:,:],axis=2)
            PT = np.mean(output.Pressure[column,:,:]/(input.Rd*output.Rho[column,:,:]) * \
                                    (output.Pressure[column,:,:]/input.P_Ref)**(-kappa_ad),axis=2)
        else:
            P = output.Pressure[column,:,0]
            PT = output.Pressure[column,:,0]/(input.Rd*output.Rho[column,:,0]) * \
                                    (output.Pressure[column,:,0]/input.P_Ref)**(-kappa_ad)

        plt.semilogy(PT,P/1e5,'k-',alpha= 0.5,lw=1)
        plt.plot(PT[np.int(np.floor(grid.nv/2))],P[np.int(np.floor(grid.nv/2))]/100000,'r+',ms =5,alpha=0.5)
        plt.plot(PT[np.int(np.floor(grid.nv*0.75))],P[np.int(np.floor(grid.nv*0.75))]/100000,'g+',ms =5,alpha=0.5)

    plt.gca().invert_yaxis()
    plt.ylabel('Pressure (bar)')
    plt.xlabel('Potential Temperature [K]')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.savefig(input.resultsf+'/figures/PTPprofile_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def streamf(input,grid,output,sigmaref):
    # eulerian stream function or mass stream function
    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)

    # Set the latitude-longitude grid.
    res_deg = 0.005
    loni, lati = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2,np.pi/2,res_deg))
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1
    ###################
    # Stream function #
    ###################

    # Initialize arrays
    vii = (output.Mh[0]*(-np.sin(grid.lat[:,None,None])*np.cos(grid.lon[:,None,None])) + \
                output.Mh[1]*(-np.sin(grid.lat[:,None,None])*np.sin(grid.lon[:,None,None])) + \
                output.Mh[2]*np.cos(grid.lat[:,None,None]))/output.Rho
    # vi = np.zeros((grid.point_num,d_sig))
    # v = np.zeros((d_lon[0],d_lon[1],d_sig,tsp))

    # Compute meridional wind
    Strii = np.zeros_like(vii)
    Stri = np.zeros((grid.point_num,d_sig))
    Stream = np.zeros((d_lon[0],d_lon[1],d_sig,tsp))

    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            sigma = output.Pressure[i,:,t]
            for j in np.arange(0,len(sigma)):
                # import pdb; pdb.set_trace()
                Strii[i,j,t] = np.trapz(vii[i,j:,t][::-1],x=sigma[j:][::-1])*np.cos(grid.lat[i])
            # Interpolate atmospheric column ot the reference pressure.
            # Need to reverse index all arrays because pchip only works if x is increasing
            Stri[i,:] = interp.pchip_interpolate(sigma[::-1],Strii[i,::-1,t],Pref[::-1])[::-1]
        # Convert icosahedral grid into lon-lat grid
        for lev in np.arange(d_sig):
            Stream[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Stri[:,lev],(loni,lati),method='nearest')

    # del vii, vi
    del Strii, Stri

    # Averaging in time and longitude.
    if tsp > 1:
        Streaml = np.mean(Stream[:,:,:,:],axis=1)
        del Stream
        Streamlt = np.mean(Streaml[:,:,:],axis=2)
        del Streaml
    else:
        Streamlt = np.mean(Stream[:,:,:,0],axis=1)
        del Stream

    # Latitude
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)

    # # Now calculate stream fxn !!!
    # Stream = np.zeros_like(vlt)
    # for i in np.arange(np.shape(vlt)[0]):
    #     for j in np.arange(1,np.shape(vlt)[1]):
    #         Stream[i,j] = np.trapz(vlt[i,::-1][-j-1:],x=Pref[::-1][-j-1:])*np.cos(latp[i])

    Streamlt *= 2*np.pi*input.A/input.Gravit

    #################
    # Create figure #
    #################

    # Contour plot
    C = plt.contourf(latp*180/np.pi,Pref/100000,Streamlt.T,40,cmap='viridis')
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0

    c2 = plt.contour(latp*180/np.pi,Pref/1e5,Streamlt.T,[0],colors='w',linewidths=1)

    plt.gca().invert_yaxis()
    if np.max(Pref)/np.min(Pref) > 100:
        plt.gca().set_yscale("log")
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (bar)')
    plt.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')

    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Stream function (kg s$^{-1}$)')
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/streamf_ver_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()


def CalcE_M_AM(input,grid,output,split):
    temperature = output.Pressure/(input.Rd*output.Rho)
    # dz = grid.Altitude[1]-grid.Altitude[0]
    # Vol = grid.areasT*dz

    ## testing @@
    tsp = output.nts-output.ntsi+1

    Atot = input.A**2
    solid_ang = grid.areasT/Atot

    rint = ((input.A+grid.Altitudeh[1:])**3-(input.A+grid.Altitudeh[:-1])**3)/3.0
    Vol0 = solid_ang[:,None]*rint[None,:]

    Eint = output.Rho*(input.Cp-input.Rd)*temperature
    Eg = output.Rho*input.Gravit*grid.Altitude[None,:,None]
    W = 0.5*(output.Wh[:,1:,:]+output.Wh[:,:-1,:])
    Wx = W*np.cos(grid.lat[:,None,None])*np.cos(grid.lon[:,None,None])
    Wy = W*np.cos(grid.lat[:,None,None])*np.sin(grid.lon[:,None,None])
    Wz = W*np.sin(grid.lat[:,None,None])

    Mtotx = output.Mh[0]+Wx
    Mtoty = output.Mh[1]+Wy
    Mtotz = output.Mh[2]+Wz
    Mtot2 = (Mtotx)**2 + (Mtoty)**2 + (Mtotz)**2

    Ek = 0.5*Mtot2/output.Rho

    output.Etotal = (Eint+Eg+Ek)*Vol0[:,:,None]

    output.Mass = output.Rho*Vol0[:,:,None]

    r = input.A+grid.Altitude

    rx = r[None,:,None]*np.cos(grid.lat[:,None,None])*np.cos(grid.lon[:,None,None])
    ry = r[None,:,None]*np.cos(grid.lat[:,None,None])*np.sin(grid.lon[:,None,None])
    rz = r[None,:,None]*np.sin(grid.lat[:,None,None])

    output.AngMomx = (ry*Mtotz-rz*Mtoty-output.Rho*input.Omega[0]*\
                     rz*r[None,:,None]*np.cos(grid.lat[:,None,None])*\
                     np.cos(grid.lon[:,None,None]))*Vol0[:,:,None]
    output.AngMomy = (-rx*Mtotz+rz*Mtotx-output.Rho*input.Omega[0]*\
                     rz*r[None,:,None]*np.cos(grid.lat[:,None,None])*\
                     np.sin(grid.lon[:,None,None]))*Vol0[:,:,None]
    output.AngMomz = (rx*Mtoty-ry*Mtotx + output.Rho*input.Omega[0]*\
                     r[None,:,None]**2*np.cos(grid.lat[:,None,None])*\
                     np.cos(grid.lat[:,None,None]))*Vol0[:,:,None]

    if split == False:
        output.GlobalE = np.sum(np.sum(output.Etotal,0),0)
        output.GlobalMass = np.sum(np.sum(output.Mass,0),0)
        output.GlobalAMx = np.sum(np.sum(output.AngMomx,0),0)
        output.GlobalAMy = np.sum(np.sum(output.AngMomy,0),0)
        output.GlobalAMz = np.sum(np.sum(output.AngMomz,0),0)
    else:
        output.WeatherE = np.zeros(np.shape(output.Pressure)[2])
        output.DeepE = np.zeros(np.shape(output.Pressure)[2])
        output.WeatherMass = np.zeros(np.shape(output.Pressure)[2])
        output.DeepMass = np.zeros(np.shape(output.Pressure)[2])
        output.WeatherAMx = np.zeros(np.shape(output.Pressure)[2])
        output.DeepAMx = np.zeros(np.shape(output.Pressure)[2])
        output.WeatherAMy = np.zeros(np.shape(output.Pressure)[2])
        output.DeepAMy = np.zeros(np.shape(output.Pressure)[2])
        output.WeatherAMz = np.zeros(np.shape(output.Pressure)[2])
        output.DeepAMz = np.zeros(np.shape(output.Pressure)[2])

        for ii in np.arange(np.shape(output.Pressure)[2]):
            output.WeatherE[ii] = np.sum(output.Etotal[output.Pressure[:,:,ii]<=split][:,ii])
            output.DeepE[ii] = np.sum(output.Etotal[output.Pressure[:,:,ii]>split][:,ii])
            output.WeatherMass[ii] = np.sum(output.Mass[output.Pressure[:,:,ii]<=split][:,ii])
            output.DeepMass[ii] = np.sum(output.Mass[output.Pressure[:,:,ii]>split][:,ii])
            output.WeatherAMx[ii] = np.sum(output.AngMomx[output.Pressure[:,:,ii]<=split][:,ii])
            output.DeepAMx[ii] = np.sum(output.AngMomx[output.Pressure[:,:,ii]>split][:,ii])
            output.WeatherAMy[ii] = np.sum(output.AngMomy[output.Pressure[:,:,ii]<=split][:,ii])
            output.DeepAMy[ii] = np.sum(output.AngMomy[output.Pressure[:,:,ii]>split][:,ii])
            output.WeatherAMz[ii] = np.sum(output.AngMomz[output.Pressure[:,:,ii]<=split][:,ii])
            output.DeepAMz[ii] = np.sum(output.AngMomz[output.Pressure[:,:,ii]>split][:,ii])

def CalcEntropy(input,grid,output,split):
    temperature = output.Pressure/(input.Rd*output.Rho)
    # dz = grid.Altitude[1]-grid.Altitude[0]
    # Vol = grid.areasT*dz

    Atot = input.A**2
    solid_ang = grid.areasT/Atot

    rint = ((input.A+grid.Altitudeh[1:])**3-(input.A+grid.Altitudeh[:-1])**3)/3.0
    Vol0 = solid_ang[:,None]*rint[None,:]

    kappa = input.Rd/input.Cp

    potT = temperature*(input.P_Ref/output.Pressure)**kappa
    S = input.Cp * np.log(potT)
    output.Entropy = S*Vol0[:,:,None]
    if split == False:
        output.GlobalEnt = np.sum(np.sum(output.Entropy,0),0)
    else:
        output.WeatherEnt = np.zeros(np.shape(output.Pressure)[2])
        output.DeepEnt = np.zeros(np.shape(output.Pressure)[2])
        for ii in np.arange(np.shape(output.Pressure)[2]):
            output.WeatherEnt[ii] = np.sum(output.Entropy[output.Pressure[:,:,ii]<=split][:,ii])
            output.DeepEnt[ii] = np.sum(output.Entropy[output.Pressure[:,:,ii]>split][:,ii])

def conservation(input,grid,output,split):
    # plot quantities that are interesting for conservation
    if split == False:
        # if (output.ConvData == False).any():
        print('Calculating energy, mass, angular momentum...')
        CalcE_M_AM(input,grid,output,split)
        plots = ['global']

    else:
        CalcE_M_AM(input,grid,output,split)
        plots = ['weather', 'deep']

    CalcEntropy(input,grid,output,split) #should put in Thor at some point?

    for ii in np.arange(len(plots)):
        fig = plt.figure(figsize=(12,8))
        fig.suptitle('Time = %#.3f - %#.3f days, split = %e bar'%(output.time[0],output.time[-1],split/1e5))
        fig.subplots_adjust(wspace=0.25,left=0.07,right=0.98,top=0.94,bottom=0.07)
        plt.subplot(2,3,1)
        if split == False:
            plt.plot(output.time,output.GlobalE,'ko',linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time,output.WeatherE,'bo',linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time,output.DeepE,'ro',linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel('Total energy of atmosphere (J)')

        plt.subplot(2,3,2)
        if split == False:
            plt.plot(output.time,output.GlobalMass,'ko',linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time,output.WeatherMass,'bo',linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time,output.DeepMass,'ro',linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel('Total mass of atmosphere (kg)')

        plt.subplot(2,3,3)
        if split == False:
            plt.plot(output.time,output.GlobalAMz,'ko',linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time,output.WeatherAMz,'bo',linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time,output.DeepAMz,'ro',linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel(r'Z angular momentum of atmosphere (kg m$^2$ s$^{-1}$)')

        plt.subplot(2,3,4)
        if split == False:
            plt.plot(output.time,output.GlobalEnt,'ko',linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time,output.WeatherEnt,'bo',linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time,output.DeepEnt,'ro',linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel('Total entropy of atmosphere (J K$^{-1}$)')

        plt.subplot(2,3,5)
        if split == False:
            plt.plot(output.time,output.GlobalAMx,'ko',linestyle='--')
            plt.plot(output.time,output.GlobalAMy,'o',color='0.5',linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time,output.WeatherAMx,'bo',linestyle='--')
                plt.plot(output.time,output.WeatherAMy,'co',linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time,output.DeepAMx,'ro',linestyle='--')
                plt.plot(output.time,output.DeepAMy,'o',color='orange',linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel(r'X, Y angular momentum of atmosphere (kg m$^2$ s$^{-1}$)')

        plt.subplot(2,3,6)
        if split == False:
            plt.plot(output.time,np.sqrt(output.GlobalAMx**2+output.GlobalAMy**2),'ko',linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time,np.sqrt(output.WeatherAMx**2+output.WeatherAMy**2),'bo',linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time,np.sqrt(output.DeepAMx**2+output.DeepAMy**2),'ro',linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel(r'Horizontal angular momentum of atmosphere (kg m$^2$ s$^{-1}$)')

        if not os.path.exists(input.resultsf+'/figures'):
            os.mkdir(input.resultsf+'/figures')
        plt.savefig(input.resultsf+'/figures/conservation_s%e_i%d_l%d_%s.pdf'%(split/100,output.ntsi,output.nts,plots[ii]))
        plt.close()

def SRindex(input,grid,output):
    tsp = output.nts-output.ntsi+1
    #total surface area
    Atot = 4*np.pi*input.A**2
    solid_ang = grid.areasT/Atot

    U = (-output.Mh[0]*np.sin(grid.lon[:,None,None])+output.Mh[1]*np.cos(grid.lon[:,None,None]))/output.Rho

    # should I adjust radius for height?? yes!!
    mom = (input.A+grid.Altitude[None,:,None])*np.cos(grid.lat[:,None,None])*\
          ((input.A+grid.Altitude[None,:,None])*input.Omega*np.cos(grid.lat[:,None,None])+U)
    mom0 = (input.A+grid.Altitude[None,:,None])*np.cos(grid.lat[:,None,None])*\
          ((input.A+grid.Altitude[None,:,None])*input.Omega*np.cos(grid.lat[:,None,None])+np.zeros(np.shape(U[:,:,:])))

    rint = ((input.A+grid.Altitudeh[1:])**3-(input.A+grid.Altitudeh[:-1])**3)/3.0
    mom_grid = mom*rint[None,:,None]*output.Rho
    mom_grid0 = mom0*rint[None,:,None]*output.Rho

    Mtot = np.sum(np.sum(mom_grid*solid_ang[:,None,None],axis=0),axis=0)
    Mtot0 = np.sum(np.sum(mom_grid0*solid_ang[:,None,None],axis=0),axis=0)

    S = Mtot/Mtot0 - 1

    plt.semilogy(output.time, S)
    plt.xlabel('Time (days)')
    plt.ylabel('Super-rotation index')
    plt.tight_layout()

    if not os.path.exists(input.resultsf+'/figures'):
            os.mkdir(input.resultsf+'/figures')
    plt.savefig(input.resultsf+'/figures/SRindex_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()
