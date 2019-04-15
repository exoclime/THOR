import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.axes as axes
import matplotlib.cm as cm
import scipy.interpolate as interp
import os
import h5py
import time
import subprocess as spr
import pyshtools as chairs
import pdb

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
        self.pntloc = np.reshape(openh5['pntloc'][...],(self.point_num,6)) # indices of nearest neighbors
        self.nv = np.int(openh5['nv'][0])
        if 'grad' in openh5.keys():
            self.grad = np.reshape(openh5['grad'][...],(self.point_num,7,3))
            self.div = np.reshape(openh5['div'][...],(self.point_num,7,3))
        if 'curlz' in openh5.keys():
            self.curlz = np.reshape(openh5['curlz'][...],(self.point_num,7,3))
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
        self.Entropy = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.Mass = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.AngMomx = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.AngMomy = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.AngMomz = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.GlobalE = np.zeros(nts-ntsi+1)
        self.GlobalEnt = np.zeros(nts-ntsi+1)
        self.GlobalMass = np.zeros(nts-ntsi+1)
        self.GlobalAMx = np.zeros(nts-ntsi+1)
        self.GlobalAMy = np.zeros(nts-ntsi+1)
        self.GlobalAMz = np.zeros(nts-ntsi+1)
        self.ConvData = np.zeros(nts-ntsi+1)
        self.EntData = np.zeros(nts-ntsi+1)
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
            if 'Entropy' in openh5.keys():
                Entropyi = openh5['Entropy'][...]
                self.GlobalEnt[t-ntsi+1] = openh5['GlobalEnt'][0]
                self.EntData[t-ntsi+1] = True
            else:
                print('Warning: entropy not available in file %s'%fileh5)
                self.EntData[t-ntsi+1] = False
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
            if 'Entropyi' in locals():
                self.Entropy[:,:,t-ntsi+1] = np.reshape(Entropyi,(grid.point_num,grid.nv))

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
    def __init__(self,resultsf,simID,ntsi,nts,input,output,grid,pressure_vert=True):
        RT = 0
        if "core_benchmark" in dir(input):
            if input.core_benchmark[0] == 0: # need to switch to 'radiative_tranfer' flag
                RT = 1

        chem = 0
        if hasattr(input,"chemistry"):
            if input.chemistry == 1:
                chem = 1

        # Read model results
        for t in np.arange(ntsi-1,nts):
            if pressure_vert == True:
                fileh5 = resultsf+'/regrid_'+simID+'_'+np.str(t+1)+'.h5'
            else:
                fileh5 = resultsf+'/regrid_height_'+simID+'_'+np.str(t+1)+'.h5'
            if os.path.exists(fileh5):
                openh5 = h5py.File(fileh5)
            else:
                print(fileh5+' not found, regridding now with default settings...')
                regrid(resultsf,simID,ntsi,nts,pressure_vert=pressure_vert)
                openh5 = h5py.File(fileh5)

            # if not 'del_hseq' in openh5.keys():
            #     openh5.close()
            #     regrid(resultsf,simID,ntsi,nts,pressure_vert=pressure_vert)
            #     openh5 = h5py.File(fileh5)

            Rhoi = openh5['Rho'][...]
            Tempi = openh5['Temperature'][...]
            Ui = openh5['U'][...]
            Vi = openh5['V'][...]
            Wi = openh5['W'][...]
            # del_hseqi = openh5['del_hseq'][...]
            if pressure_vert == True:
                Prei = openh5['Pressure'][...]
            else:
                Alti = openh5['Altitude'][...]
            lati = openh5['Latitude'][...]
            loni = openh5['Longitude'][...]

            if RT == 1:
                tau_swi = openh5['tau_sw'][...]
                tau_lwi = openh5['tau_lw'][...]
                fnet_upi = openh5['fnet_up'][...]
                fnet_dni = openh5['fnet_dn'][...]

            if chem == 1:
                ch4i = openh5['ch4'][...]
                coi = openh5['co'][...]
                h2oi = openh5['h2o'][...]
                co2i = openh5['co2'][...]
                nh3i = openh5['nh3'][...]

            if not 'PV' in openh5.keys():
                calc_RV_PV(grid,output,input,loni,lati,Prei,t-ntsi+1,fileh5,pressure_vert=pressure_vert)

            PVi = openh5['PV'][...]  #these are regridded at low resolution because they are a pita to calculate
            RVi = openh5['RV'][...]
            lati_lr = openh5['Lat_lowres'][...]
            loni_lr = openh5['Lon_lowres'][...]

            # if not 'streamf' in openh5.keys():
            #     calc_moc_streamf(grid,output,input,loni,lati,Prei,t-ntsi+1,fileh5,pressure_vert=pressure_vert)

            # streamfi = openh5['streamf'][...]

            openh5.close()

            if t == ntsi-1:
                self.Rho = np.zeros(np.shape(Rhoi)+(nts-ntsi+1,))
                self.U = np.zeros(np.shape(Ui)+(nts-ntsi+1,))
                self.V = np.zeros(np.shape(Vi)+(nts-ntsi+1,))
                self.W = np.zeros(np.shape(Wi)+(nts-ntsi+1,))
                self.Temperature = np.zeros(np.shape(Tempi)+(nts-ntsi+1,))
                if pressure_vert == True:
                    self.Pressure = np.zeros(np.shape(Prei)+(nts-ntsi+1,))
                else:
                    self.Altitude = np.zeros(np.shape(Alti)+(nts-ntsi+1,))
                self.lat = np.zeros(np.shape(lati)+(nts-ntsi+1,))
                self.lon = np.zeros(np.shape(loni)+(nts-ntsi+1,))
                self.PV = np.zeros(np.shape(PVi)+(nts-ntsi+1,))
                self.RV = np.zeros(np.shape(RVi)+(nts-ntsi+1,))
                self.lat_lr = np.zeros(np.shape(lati_lr)+(nts-ntsi+1,))
                self.lon_lr = np.zeros(np.shape(loni_lr)+(nts-ntsi+1,))
                # self.streamf = np.zeros(np.shape(streamfi)+(nts-ntsi+1,))

                # self.del_hseq = np.zeros(np.shape(del_hseqi)+(nts-ntsi+1,))
                if RT == 1:
                    self.tau_sw = np.zeros(np.shape(tau_swi)+(nts-ntsi+1,))
                    self.tau_lw = np.zeros(np.shape(tau_lwi)+(nts-ntsi+1,))
                    self.fnet_up = np.zeros(np.shape(fnet_upi)+(nts-ntsi+1,))
                    self.fnet_dn = np.zeros(np.shape(fnet_dni)+(nts-ntsi+1,))
                if chem == 1:
                    self.ch4 = np.zeros(np.shape(ch4i)+(nts-ntsi+1,))
                    self.co = np.zeros(np.shape(coi)+(nts-ntsi+1,))
                    self.h2o = np.zeros(np.shape(h2oi)+(nts-ntsi+1,))
                    self.co2 = np.zeros(np.shape(co2i)+(nts-ntsi+1,))
                    self.nh3 = np.zeros(np.shape(nh3i)+(nts-ntsi+1,))

            self.Rho[:,:,:,t-ntsi+1] = Rhoi
            self.U[:,:,:,t-ntsi+1] = Ui
            self.V[:,:,:,t-ntsi+1] = Vi
            self.W[:,:,:,t-ntsi+1] = Wi
            self.Temperature[:,:,:,t-ntsi+1] = Tempi
            # self.del_hseq[:,:,:,t-ntsi+1] = del_hseqi
            if pressure_vert == True:
                self.Pressure[:,t-ntsi+1] = Prei
            else:
                self.Altitude[:,t-ntsi+1] = Alti
            self.lat[:,t-ntsi+1] = lati
            self.lon[:,t-ntsi+1] = loni
            self.PV[:,:,:,t-ntsi+1] = PVi
            self.RV[:,:,:,:,t-ntsi+1] = RVi
            self.lat_lr[:,t-ntsi+1] = lati_lr
            self.lon_lr[:,t-ntsi+1] = loni_lr
            # self.streamf[:,:,:,t-ntsi+1] = streamfi

            if RT == 1:
                self.tau_sw[:,:,:,t-ntsi+1] = tau_swi
                self.tau_lw[:,:,:,t-ntsi+1] = tau_lwi
                self.fnet_up[:,:,:,t-ntsi+1] = fnet_upi
                self.fnet_dn[:,:,:,t-ntsi+1] = fnet_dni
            if chem == 1:
                self.ch4[:,:,:,t-ntsi+1] = ch4i
                self.co[:,:,:,t-ntsi+1] = coi
                self.h2o[:,:,:,t-ntsi+1] = h2oi
                self.co2[:,:,:,t-ntsi+1] = co2i
                self.nh3[:,:,:,t-ntsi+1] = nh3i


class GetOutput:
    def __init__(self,resultsf,simID,ntsi,nts,stride=1,openrg=0,pressure_vert=True):
        self.input = input(resultsf,simID)
        self.grid = grid(resultsf,simID)
        self.output = output(resultsf,simID,ntsi,nts,self.grid,stride=stride)
        if openrg == 1:
            self.rg = rg_out(resultsf,simID,ntsi,nts,self.input,self.output,self.grid,pressure_vert=pressure_vert)

def calc_RV_PV(grid,output,input,lons,lats,sigma,t_ind,fileh5,comp=4,pressure_vert=True):
    #Calculates relative and potential vorticity on height levels, then interpolates to pressure.
    #It is easiest to calculate these on a lat-lon-altitude grid since conversion
    #to pressure requires transformation of vertical velocity and computing 3D curls
    #at the grid level would require heavy interpolation and invocation of stokes thm
    ang_res = 4.0/2**(input.glevel-4)
    # first interpolate to near thor grid resolution

    lat_range_tmp = np.arange(-90,90+ang_res,ang_res)
    lat_range = (lat_range_tmp[:-1]+lat_range_tmp[1:])/2
    lon_range = np.arange(0,360,ang_res)

    loni, lati, alti = np.meshgrid(lon_range,lat_range,grid.Altitude)
    d_lon = np.shape(loni)
    res_deg = ang_res*np.pi/180

    U = (-output.Mh[0,:,:,t_ind]*np.sin(grid.lon[:,None])+output.Mh[1,:,:,t_ind]*np.cos(grid.lon[:,None]))/output.Rho[:,:,t_ind]
    V = (-output.Mh[0,:,:,t_ind]*np.sin(grid.lat[:,None])*np.cos(grid.lon[:,None])\
                      -output.Mh[1,:,:,t_ind]*np.sin(grid.lat[:,None])*np.sin(grid.lon[:,None])\
                      +output.Mh[2,:,:,t_ind]*np.cos(grid.lat[:,None]))/output.Rho[:,:,t_ind]
    interpx = (grid.Altitude-grid.Altitudeh[:-1])/(grid.Altitudeh[1:]-grid.Altitudeh[:-1])
    W = (output.Wh[:,:-1,t_ind] + (output.Wh[:,1:,t_ind]-output.Wh[:,:-1,t_ind])*interpx[None,:])/output.Rho[:,:,t_ind]
    Temp_icoh = output.Pressure[:,:,t_ind]/(input.Rd*output.Rho[:,:,t_ind])
    pt = Temp_icoh*(output.Pressure[:,:,t_ind]/input.P_Ref)**(-input.Rd/input.Cp)

    #full meshgrid
    lonf, latf = np.meshgrid(lons,lats)

    pt_llz = np.zeros((d_lon[0],d_lon[1],len(grid.Altitude)))
    rho_llz = np.zeros((d_lon[0],d_lon[1],len(grid.Altitude)))
    p_llz = np.zeros((d_lon[0],d_lon[1],len(grid.Altitude)))
    U_llz = np.zeros((d_lon[0],d_lon[1],len(grid.Altitude)))
    V_llz = np.zeros((d_lon[0],d_lon[1],len(grid.Altitude)))
    W_llz = np.zeros((d_lon[0],d_lon[1],len(grid.Altitude)))

    for lev in np.arange(len(grid.Altitude)):
        # interpolate onto lat-lon grid
        pt_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,pt[:,lev],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
        rho_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,output.Rho[:,lev,t_ind],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
        p_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,output.Pressure[:,lev,t_ind],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
        U_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,U[:,lev],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
        V_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,V[:,lev],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
        W_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,W[:,lev],(loni[:,:,lev],lati[:,:,lev]),method='nearest')

    # Compute gradient of potential temp
    dz = np.abs(grid.Altitude[1] - grid.Altitude[0])
    dptdz = dFunc_dZ(pt_llz,alti,grid.nv,dz)

    dptdlat = dFunc_dLat(pt_llz,res_deg)/(input.A+alti)

    dptdlon = dFunc_dLon(pt_llz,res_deg)/(input.A+alti)/np.cos(lati)

    # Compute curl of wind speeds
    curlVz, curlVlat, curlVlon = CurlF(W_llz,V_llz,U_llz,lati*np.pi/180,alti,res_deg,grid.nv,input.A,dz)
    RV_llz = np.zeros((3,d_lon[0],d_lon[1],len(grid.Altitude)))
    RV_llz[0] = curlVz
    RV_llz[1] = curlVlat
    RV_llz[2] = curlVlon

    curlVlon += 2*input.Omega

    curlVz /= rho_llz
    curlVlat /= rho_llz
    curlVlon /= rho_llz

    PV_llz = curlVz*dptdz + curlVlat*dptdlat + curlVlon*dptdlon

    if pressure_vert == True:
        PV_llp = np.zeros((d_lon[0],d_lon[1], len(sigma)))
        RV_llp = np.zeros((3,d_lon[0],d_lon[1], len(sigma)))

        for ilat in np.arange(d_lon[0]):
            for ilon in np.arange(d_lon[1]):
                # interpolate to pressure grid
                PV_llp[ilat,ilon,:] = interp.pchip_interpolate(p_llz[ilat,ilon,:][::-1],
                                                PV_llz[ilat,ilon,:][::-1],sigma[::-1])[::-1]
                RV_llp[0,ilat,ilon,:] = interp.pchip_interpolate(p_llz[ilat,ilon,:][::-1],
                                                RV_llz[0,ilat,ilon,:][::-1],sigma[::-1])[::-1]
                RV_llp[1,ilat,ilon,:] = interp.pchip_interpolate(p_llz[ilat,ilon,:][::-1],
                                                RV_llz[1,ilat,ilon,:][::-1],sigma[::-1])[::-1]
                RV_llp[2,ilat,ilon,:] = interp.pchip_interpolate(p_llz[ilat,ilon,:][::-1],
                                                RV_llz[2,ilat,ilon,:][::-1],sigma[::-1])[::-1]

        #testing something...
        PV_llp_f = PV_llp
        RV_llp_f = RV_llp
    else:
        PV_llp_f = PV_llz
        RV_llp_f = PV_llz

    # # grid now to higher resolution
    # PV_llp_f = np.zeros(np.shape(lonf)+(len(sigma),))
    # RV_llp_f = np.zeros((3,)+np.shape(lonf)+(len(sigma),))
    #
    # for lev in np.arange(len(sigma)):
    #     # interpolate onto higher resolution grid
    #     locs = np.vstack([np.ravel(loni[:,:,lev]),np.ravel(lati[:,:,lev])]).T
    #     PV_llp_f[:,:,lev] = interp.griddata(locs,np.ravel(PV_llp[:,:,lev]),(lonf,latf),method='nearest')
    #     RV_llp_f[0,:,:,lev] = interp.griddata(locs,np.ravel(RV_llp[0,:,:,lev]),(lonf,latf),method='nearest')
    #     RV_llp_f[1,:,:,lev] = interp.griddata(locs,np.ravel(RV_llp[1,:,:,lev]),(lonf,latf),method='nearest')
    #     RV_llp_f[2,:,:,lev] = interp.griddata(locs,np.ravel(RV_llp[2,:,:,lev]),(lonf,latf),method='nearest')

    openh5 = h5py.File(fileh5,"r+")
    print('Adding vorticities to regrid file...')
    PV = openh5.create_dataset("PV",data=PV_llp_f,compression='gzip',compression_opts=comp)
    RV = openh5.create_dataset("RV",data=RV_llp_f,compression='gzip',compression_opts=comp)
    Latlr = openh5.create_dataset("Lat_lowres",data=lat_range,compression='gzip',compression_opts=comp)
    Lonlr = openh5.create_dataset("Lon_lowres",data=lon_range,compression='gzip',compression_opts=comp)
    openh5.close()

def regrid(resultsf,simID,ntsi,nts,nlev=40,pscale='log',overwrite=False,comp=4,pressure_vert=True):
    # runs over files and converts ico-height grid to lat-lon-pr grid
    outall = GetOutput(resultsf,simID,ntsi,nts)
    input = outall.input
    grid = outall.grid
    output = outall.output
    res_deg = 4.0/2**(input.glevel-4)

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
    Pref = input.P_Ref*sigmaref[:,0]
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
    del_hseq = np.gradient(output.Pressure,grid.Altitude,axis=1) + output.Rho*input.Gravit

    RT = 0
    if hasattr(input,"core_benchmark"):
        if input.core_benchmark[0] == 0: # need to switch to 'radiative_tranfer' flag
            RT = 1
            fnet_up_icoh = output.fnet_up[:,:-1,:]+(output.fnet_up[:,1:,:]-output.fnet_up[:,:-1,:])*interpx[None,:,None]
            fnet_dn_icoh = output.fnet_dn[:,:-1,:]+(output.fnet_dn[:,1:,:]-output.fnet_dn[:,:-1,:])*interpx[None,:,None]

    chem = 0
    if hasattr(input,"chemistry"):
        if input.chemistry == 1:
            chem = 1

    for t in np.arange(ntsi,nts+1):
        #check for exising h5 files
        proceed = 0
        if pressure_vert == True:
            fileh5 = resultsf+'/regrid_'+simID+'_'+np.str(t)+'.h5'
        else:
            fileh5 = resultsf+'/regrid_height_'+simID+'_'+np.str(t)+'.h5'
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

            del_hseq_icop = np.zeros((grid.point_num,d_sig))
            del_hseq_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            if RT == 1:
                tau_sw_icop = np.zeros((grid.point_num,d_sig))
                tau_sw_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                tau_lw_icop = np.zeros((grid.point_num,d_sig))
                tau_lw_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                fnet_up_icop = np.zeros((grid.point_num,d_sig))
                fnet_up_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                fnet_dn_icop = np.zeros((grid.point_num,d_sig))
                fnet_dn_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                insol_ll = np.zeros((d_lon[0],d_lon[1]))


            if chem == 1:
                ch4_icop = np.zeros((grid.point_num,d_sig))
                ch4_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                co_icop = np.zeros((grid.point_num,d_sig))
                co_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                h2o_icop = np.zeros((grid.point_num,d_sig))
                h2o_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                co2_icop = np.zeros((grid.point_num,d_sig))
                co2_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                nh3_icop = np.zeros((grid.point_num,d_sig))
                nh3_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            if pressure_vert == True: # use pressure as vertical coordinate
                for i in np.arange(grid.point_num):
                    #interp to pressure grid
                    sigma = output.Pressure[i,:,t-ntsi]
                    Temp_icop[i,:] = interp.pchip_interpolate(sigma[::-1],Temp_icoh[i,::-1,t-ntsi],Pref[::-1])[::-1]
                    Rho_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.Rho[i,::-1,t-ntsi],Pref[::-1])[::-1]
                    Mh_icop[0,i,:] = interp.pchip_interpolate(sigma[::-1],output.Mh[0,i,::-1,t-ntsi],Pref[::-1])[::-1]
                    Mh_icop[1,i,:] = interp.pchip_interpolate(sigma[::-1],output.Mh[1,i,::-1,t-ntsi],Pref[::-1])[::-1]
                    Mh_icop[2,i,:] = interp.pchip_interpolate(sigma[::-1],output.Mh[2,i,::-1,t-ntsi],Pref[::-1])[::-1]
                    Wh_icop[i,:] = interp.pchip_interpolate(sigma[::-1],Wh_icoh[i,::-1,t-ntsi],Pref[::-1])[::-1]
                    del_hseq_icop[i,:] = interp.pchip_interpolate(sigma[::-1],del_hseq[i,::-1,t-ntsi],Pref[::-1])[::-1]
                    if RT == 1:
                        tau_sw_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.tau_sw[i,::-1,t-ntsi],Pref[::-1])[::-1]
                        tau_lw_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.tau_lw[i,::-1,t-ntsi],Pref[::-1])[::-1]
                        fnet_up_icop[i,:] = interp.pchip_interpolate(sigma[::-1],fnet_up_icoh[i,::-1,t-ntsi],Pref[::-1])[::-1]
                        fnet_dn_icop[i,:] = interp.pchip_interpolate(sigma[::-1],fnet_dn_icoh[i,::-1,t-ntsi],Pref[::-1])[::-1]
                    if chem == 1:
                        ch4_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.ch4[i,::-1,t-ntsi],Pref[::-1])[::-1]
                        co_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.co[i,::-1,t-ntsi],Pref[::-1])[::-1]
                        h2o_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.h2o[i,::-1,t-ntsi],Pref[::-1])[::-1]
                        co2_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.co2[i,::-1,t-ntsi],Pref[::-1])[::-1]
                        nh3_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.nh3[i,::-1,t-ntsi],Pref[::-1])[::-1]
            else:  # keep height at vertical coordinate (sometimes useful)
                Temp_icop[:,:] = Temp_icoh[:,:,t-ntsi]
                Rho_icop[:,:] = output.Rho[:,:,t-ntsi]
                Mh_icop[:,:,:] = output.Mh[:,:,:,t-ntsi]
                Wh_icop[:,:] = Wh_icoh[:,:,t-ntsi]
                del_hseq_icop[:,:] = del_hseq[:,:,t-ntsi]
                if RT == 1:
                    tau_sw_icop[:,:] = output.tau_sw[:,:,t-ntsi]
                    tau_lw_icop[:,:] = output.tau_lw[:,:,t-ntsi]
                    fnet_up_icop[:,:] = fnet_up_icoh[:,:,t-ntsi]
                    fnet_dn_icop[:,:] = fnet_dn_icoh[:,:,t-ntsi]
                if chem == 1:
                    ch4_icop[:,:] = output.ch4[:,:,t-ntsi]
                    co_icop[:,:] = output.co[:,:,t-ntsi]
                    h2o_icop[:,:] = output.h2o[:,:,t-ntsi]
                    co2_icop[:,:] = output.co2[:,:,t-ntsi]
                    nh3_icop[:,:] = output.nh3[:,:,t-ntsi]

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
                del_hseq_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,del_hseq_icop[:,lev],(loni,lati),method='nearest')
                if RT == 1:
                    tau_sw_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,tau_sw_icop[:,lev],(loni,lati),method='nearest')
                    tau_lw_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,tau_lw_icop[:,lev],(loni,lati),method='nearest')
                    fnet_up_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,fnet_up_icop[:,lev],(loni,lati),method='nearest')
                    fnet_dn_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,fnet_dn_icop[:,lev],(loni,lati),method='nearest')
                    insol_ll[:,:] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,output.Insol[:,t-ntsi],(loni,lati),method='nearest')
                if chem == 1:
                    ch4_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,ch4_icop[:,lev],(loni,lati),method='nearest')
                    co_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,co_icop[:,lev],(loni,lati),method='nearest')
                    h2o_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,h2o_icop[:,lev],(loni,lati),method='nearest')
                    co2_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,co2_icop[:,lev],(loni,lati),method='nearest')
                    nh3_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,nh3_icop[:,lev],(loni,lati),method='nearest')

            #create h5 files
            openh5 = h5py.File(fileh5,"w")

            print('Writing file '+fileh5+'...')
            # coordinates
            if pressure_vert == True:
                Pre = openh5.create_dataset("Pressure",data=Pref,compression='gzip',compression_opts=comp)
            else:
                Alt = openh5.create_dataset("Altitude",data=grid.Altitude,compression='gzip',compression_opts=comp)
            Lat = openh5.create_dataset("Latitude",data=lat_range,compression='gzip',compression_opts=comp)
            Lon = openh5.create_dataset("Longitude",data=lon_range,compression='gzip',compression_opts=comp)
            # data
            Temp = openh5.create_dataset("Temperature",data=Temp_llp,compression='gzip',compression_opts=comp)

            Rho = openh5.create_dataset("Rho",data=Rho_llp,compression='gzip',compression_opts=comp)
            U = openh5.create_dataset("U",data=U_llp,compression='gzip',compression_opts=comp)
            V = openh5.create_dataset("V",data=V_llp,compression='gzip',compression_opts=comp)
            W = openh5.create_dataset("W",data=W_llp,compression='gzip',compression_opts=comp)
            del_hseqf = openh5.create_dataset("del_hseq",data=del_hseq_llp,compression='gzip',compression_opts=comp)
            if RT == 1:
                tau_sw = openh5.create_dataset("tau_sw",data=tau_sw_llp,compression='gzip',compression_opts=comp)
                tau_lw = openh5.create_dataset("tau_lw",data=tau_lw_llp,compression='gzip',compression_opts=comp)
                fnet_up = openh5.create_dataset("fnet_up",data=fnet_up_llp,compression='gzip',compression_opts=comp)
                fnet_dn = openh5.create_dataset("fnet_dn",data=fnet_dn_llp,compression='gzip',compression_opts=comp)
                insol = openh5.create_dataset("insol",data=insol_ll,compression='gzip',compression_opts=comp)
            if chem == 1:
                ch4 = openh5.create_dataset("ch4",data=ch4_llp,compression='gzip',compression_opts=comp)
                co = openh5.create_dataset("co",data=co_llp,compression='gzip',compression_opts=comp)
                h2o = openh5.create_dataset("h2o",data=h2o_llp,compression='gzip',compression_opts=comp)
                co2 = openh5.create_dataset("co2",data=co2_llp,compression='gzip',compression_opts=comp)
                nh3 = openh5.create_dataset("nh3",data=nh3_llp,compression='gzip',compression_opts=comp)

            openh5.close()

            ## calculate relative and potential vorticity and add to regrid file
            calc_RV_PV(grid,output,input,lon_range,lat_range,Pref,t-ntsi,fileh5,pressure_vert)
            # calc_moc_streamf(grid,output,input,lon_range,lat_range,Pref,t-ntsi,fileh5,pressure_vert)


def Get_Prange(input,grid,output,args,xtype='lat'):
    # Sigma (normalized pressure) values for the plotting
    if not isinstance(args.slice,list):
        raise IOError("'slice' argument must be a list")

    if (args.pressure_min[0]=='default'):
        if xtype == 'lat':
            if len(args.slice) == 2:
                grid_mask = np.logical_and(grid.lon*180/np.pi>=args.slice[0],grid.lon*180/np.pi<=args.slice[1])
                if (grid_mask==False).all():
                    #there were no grid points in the range
                    grid_mask = np.argmin(np.abs(grid.lon*180/np.pi-0.5*(args.slice[0]+args.slice[1])))
            elif len(args.slice) == 1:
                grid_mask = np.argmin(np.abs(grid.lon*180/np.pi-args.slice[0]))
        elif xtype == 'lon':
            if len(args.slice) == 2:
                grid_mask = np.logical_and(grid.lat*180/np.pi>=args.slice[0],grid.lat*180/np.pi<=args.slice[1])
                if (grid_mask==False).all():
                    #there were no grid points in the range
                    grid_mask = np.argmin(np.abs(grid.lat*180/np.pi-0.5*(args.slice[0]+args.slice[1])))
            elif len(args.slice) == 1:
                grid_mask = np.argmin(np.abs(grid.lat*180/np.pi-args.slice[0]))
        args.pressure_min[0] = np.max(output.Pressure[grid_mask,grid.nv-1,:])/100

    if np.max(input.P_Ref)/np.float(args.pressure_min[0]) > 1000:
        sigmaref = np.logspace(np.log10(input.P_Ref),np.log10(np.float(args.pressure_min[0])*100),20)/input.P_Ref
    else:
        sigmaref = np.linspace(input.P_Ref,np.float(args.pressure_min[0])*100,20)/input.P_Ref

    return sigmaref


def KE_spect(input,grid,output,sigmaref,coord = 'icoh',lmax_adjust = 0):
    tsp = output.nts-output.ntsi+1
    lmax_grid = np.int(np.floor(np.sqrt(grid.point_num))/2-1)

    if coord == 'icoh':
        W = 0.5*(output.Wh[:,1:,:]+output.Wh[:,:-1,:])
        Wx = W*np.cos(grid.lat[:,None,None])*np.cos(grid.lon[:,None,None])
        Wy = W*np.cos(grid.lat[:,None,None])*np.sin(grid.lon[:,None,None])
        Wz = W*np.sin(grid.lat[:,None,None])

        Vx = (output.Mh[0] + Wx)
        Vy = (output.Mh[1] + Wy)
        Vz = (output.Mh[2] + Wz)

        KE = 0.5*(Vx**2+Vy**2+Vz**2)
        lmax = np.int(lmax_grid+lmax_adjust) #sets lmax based on grid size

        x_coeffs = np.zeros((2,lmax+1,lmax+1,grid.nv,tsp),dtype=complex)
        y_coeffs = np.zeros((2,lmax+1,lmax+1,grid.nv,tsp),dtype=complex)
        z_coeffs = np.zeros((2,lmax+1,lmax+1,grid.nv,tsp),dtype=complex)
        KE_coeffs = np.zeros((2,lmax+1,lmax+1,grid.nv,tsp),dtype=complex)
        KE_power = np.zeros((lmax+1,grid.nv,tsp))
        waven = np.arange(lmax+1)  #total spherical wavenumber

        cmap = cm.get_cmap('cividis')
        fig, ax = plt.subplots(1, 1)

        if tsp == 1:
            for lev in np.arange(grid.nv):
                KE_coeffs[:,:,:,lev,0], chiz = chairs.expand.SHExpandLSQ(KE[:,lev,0],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                KE_power[:,lev,0] = chairs.spectralanalysis.spectrum(KE_coeffs,unit='per_lm')
                ax.plot(waven, KE_power[:,lev,0],'k-',c=cmap(lev/grid.nv),lw=1)
        else:
            for t in np.arange(tsp):
                KE_coeffs[:,:,:,grid.nv-1,t], chiz = chairs.expand.SHExpandLSQ(KE[:,grid.nv-1,t],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                KE_power[:,grid.nv-1,t] = chairs.spectralanalysis.spectrum(KE_coeffs,unit='per_lm')
                ax.plot(waven, KE_power[:,grid.nv-1,t],'k-',c=cmap(t/tsp),lw=1)

    else:
        raise IOError("Invalid coord option! Valid options are 'icoh'")

    # ax.plot(waven,np.mean(KE_power[:,:,t],axis=1),'k-',lw=2)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.vlines(lmax_grid,ax.get_ylim()[0],ax.get_ylim()[1],zorder=1000,linestyle='--')
    ax.set(ylabel='KE density (kg m$^{-1}$ s$^{-2}$)',xlabel='n')
    # ke3line = ax.get_ylim()[0] + waven**(-3.0)
    # ke53line = ax.get_ylim()[0] + waven**(-5.0/3)
    # ax.plot(waven,ke3line,'r:')
    # ax.plot(waven,ke53line,':',color='r')

    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/KEspectrum_%i_%i_%s.pdf'%(output.ntsi,output.nts,coord))
    plt.close()

    norm = colors.Normalize(vmin = np.min(KE[:,0,0]),vmax=np.max(KE[:,0,0]))

    fig = plt.figure()
    fig.subplots_adjust(left=0.1,right=0.97)
    plt.subplot(2,1,1)
    if coord == 'llp':
        plt.imshow(KE[:,:,0,0],origin='lower',extent=(0,360,-90,90),norm=norm,aspect='auto')
    if coord == 'icoh':
        plt.tricontourf(grid.lon*180/np.pi,grid.lat*180/np.pi,KE[:,0,0],levels=30)
    plt.xlabel('Longitude ($^{\circ}$)')
    plt.ylabel('Latitude ($^{\circ}$)')
    clb = plt.colorbar()
    clb.set_label('Kinetic energy density (kg m$^{-1}$ s$^{-2}$)')

    KEcomp = chairs.expand.MakeGridDH(KE_coeffs[:,:,:,0,0],sampling=2)
    lat = np.linspace(-90,90,np.shape(KEcomp)[0])
    lon = np.linspace(0,360,np.shape(KEcomp)[1])

    plt.subplot(2,1,2)
    plt.imshow(np.real(KEcomp),origin='lower',extent=(0,360,-90,90),norm=norm,aspect='auto')
    plt.xlabel('Latitude ($^{\circ}$)')
    plt.ylabel('Longitude ($^{\circ}$)')
    plt.colorbar()
    clb.set_label('Kinetic energy density (kg m$^{-1}$ s$^{-2}$)')

    plt.savefig(input.resultsf+'/figures/KEmap_lowest_%i_%s.pdf'%(output.ntsi,coord))
    plt.close()

def vertical_lat(input,grid,output,rg,sigmaref,z,slice=[0,360],save=True,axis=False,csp=500,wind_vectors=False):
    # generic pressure/latitude plot function

    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)
    tsp = output.nts-output.ntsi+1

    if not isinstance(slice,list):
        raise IOError("'slice' argument must be a list")

    lat = z['lat']
    lon = z['lon']
    if len(slice) == 2:
        # Set the latitude-longitude grid
        if slice[1]-slice[0] > 360:
            raise IOError("'slice' cannot exceed a range of 360 degrees")
        if slice[1] < slice[0]:
            raise IOError("'slice' values must be arranged (small, large) because I am too lazy to code it better")
        if slice[0] < 0:
            mask_ind = np.logical_or((lon-360)>=slice[0],lon<=slice[1])
        else:
            mask_ind = np.logical_and(lon>=slice[0],lon<=slice[1])
        loni, lati = np.meshgrid(lon[mask_ind],lat)
        d_lon = np.shape(loni)

        ##################
        #    Averages    #
        ##################

        # Averaging in time and longitude
        if tsp > 1:
            Zonall = np.mean(z['value'][:,mask_ind[:,0],:,:],axis=1)
            # Vl = np.mean(rg.V[:,:,:,:],axis=1)
            # Wl = np.mean(rg.W[:,:,:,:],axis=1)
            Zonallt = np.mean(Zonall[:,:,:],axis=2)
            # Vlt = np.mean(Vl[:,:,:],axis=2)
            # Wlt = np.mean(Wl[:,:,:],axis=2)
            del Zonall
            if wind_vectors == True:
                Vl = np.mean(rg.V[:,mask_ind[:,0],:,:],axis=1)
                Wl = np.mean(rg.W[:,mask_ind[:,0],:,:],axis=1)
                Vlt = np.mean(Vl[:,:,:],axis=2)
                Wlt = np.mean(Wl[:,:,:],axis=2)
                del Vl, Wl
        else:
            Zonallt = np.mean(z['value'][:,:,:,0][:,mask_ind[:,0],:],axis=1)
            # Vlt = np.mean(rg.V[:,:,:,0],axis=1)
            # Wlt = np.mean(rg.W[:,:,:,0],axis=1)
            if wind_vectors == True:
                Vlt = np.mean(rg.V[:,:,:,0][:,mask_ind[:,0],:],axis=1)
                Wlt = np.mean(rg.W[:,:,:,0][:,mask_ind[:,0],:],axis=1)

    elif len(slice) == 1:
        if slice[0] in lon:
            Zonall = z['value'][:,lon[:,0]==slice[0],:,:]
            if wind_vectors == True:
                Vl = rg.V[:,lon[:,0]==slice[0],:,:]
                Wl = rg.W[:,lon[:,0]==slice[0],:,:]
        else:
            Zonall = np.zeros((len(lat),1,d_sig,tsp))
            if wind_vectors == True:
                Vl = np.zeros((len(lat),1,d_sig,tsp))
                Wl = np.zeros((len(lat),1,d_sig,tsp))
            # interpolate to slice given
            for t in tsp:
                for lev in np.arange(d_sig):
                    Zonall[:,0,lev,tsp] = interp.griddata(np.vstack([lon,lat]).T,z['value'][:,:,lev,tsp],(slice[0],lat))
                    if wind_vectors == True:
                        Vl[:,0,lev,tsp] = interp.griddata(np.vstack([lon,lat]).T,rg.V[:,:,lev,tsp],(slice[0],lat))
                        Wl[:,0,lev,tsp] = interp.griddata(np.vstack([lon,lat]).T,rg.W[:,:,lev,tsp],(slice[0],lat))


        # Averaging in time
        if tsp > 1:
            Zonallt = np.mean(Zonall[:,0,:,:],axis=2)
            del Zonall
            if wind_vectors == True:
                Vlt = np.mean(Vl[:,0,:,:],axis=2)
                Wlt = np.mean(Wl[:,0,:,:],axis=2)
                del Vl, Wl
        else:
            Zonallt = Zonall[:,0,:,0]
            if wind_vectors == True:
                Vlt = Vl[:,0,:,0]
                Wlt = Wl[:,0,:,0]

    else:
        raise IOError("'slice' must have 1 or 2 values")

    #################
    # Create figure #
    #################
    # Latitude
    latp = lat[:,0]*np.pi/180

    # need to set desired pressure range (major PITA!)
    prange = np.where(np.logical_and(rg.Pressure>=np.min(Pref),rg.Pressure<=np.max(Pref)))

    # Contour plot
    if isinstance(axis,axes.SubplotBase):
        C = axis.contourf(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,Zonallt[:,prange[0]].T,40,cmap=z['cmap'])
        ax = axis
    elif axis == False:
        C = plt.contourf(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,Zonallt[:,prange[0]].T,40,cmap=z['cmap'])
        ax = plt.gca()
    else:
        raise IOError("'axis = {}' but {} is not an axes.SubplotBase instance".format(axis,axis))

    if wind_vectors == True:
        vspacing = np.int(np.shape(rg.lat)[0]/10)
        wspacing = np.int(np.shape(rg.Pressure)[0]/10)
        Vlt = Vlt[:,prange[0]]
        Wlt = Wlt[:,prange[0]]
        Vq = Vlt[::vspacing,::wspacing].ravel()
        Wq = Wlt[::vspacing,::wspacing].ravel()
        #preq = rg.Pressure[:,0][::spacing,::spacing].ravel()
        #latq = lati[::spacing,::spacing].ravel()
        latq, preq = np.meshgrid(rg.lat[::vspacing,0],rg.Pressure[::wspacing,0][prange[0]])
        del Vlt, Wlt
        plt.quiver(latq.ravel(),preq.ravel()/1e5,Vq,Wq,color='0.5')

    clb = plt.colorbar(C,extend='both',ax=ax)
    clb.set_label(z['label'])
    if isinstance(csp,list) or isinstance(csp,tuple):
        levp = csp
    else:
        levp = np.arange(np.ceil(np.min(Zonallt[:,prange[0]])/csp)*csp,np.floor(np.max(Zonallt[:,prange[0]])/csp)*csp,csp)
    c2 = ax.contour(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,Zonallt[:,prange[0]].T,levels=levp,colors='w',linewidths=1)
    plt.clabel(c2,inline=1,fontsize=10)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    ax.invert_yaxis()
    # plt.quiver(latq.ravel(),preq.ravel()/1e5,Vq/np.max(Vq),Wq/np.max(Wq),color='0.5')
    if np.max(Pref)/np.min(Pref) > 100:
        ax.set_yscale("log")
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Pressure (bar)')
    ax.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
    #    plt.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
    if np.min(rg.Pressure[prange[0],0]) < np.max(output.Pressure[:,grid.nv-1,:]):
        ax.set_ylim(np.max(rg.Pressure[prange[0],0])/1e5,np.min(rg.Pressure[prange[0],0])/1e5)
    else:
        ax.set_ylim(np.max(rg.Pressure[prange[0],0])/1e5,np.max(output.Pressure[:,grid.nv-1,:])/1e5)

    if len(slice) == 2:
        ax.set_title('Time = %#.3f-%#.3f days, Lon = (%#.3f,%#.3f)'%(output.time[0],output.time[-1],slice[0],slice[1]))
    else:
        ax.set_title('Time = %#.3f-%#.3f days, Lon = (%#.3f,)'%(output.time[0],output.time[-1],slice[0]))

    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    if save == True:
        # save the plot to file designated by z
        if len(slice)==2:
            plt.savefig(input.resultsf+'/figures/%s_ver_i%d_l%d_lon%#.2f-%#.2f.pdf'%(z['name'],output.ntsi,output.nts,slice[0],slice[1]))
        else:
            plt.savefig(input.resultsf+'/figures/%s_ver_i%d_l%d_lon%#.2f.pdf'%(z['name'],output.ntsi,output.nts,slice[0]))
        plt.close()

def horizontal_lev(input,grid,output,rg,Plev,z,save=True,axis=False,wind_vectors=False):
    # Set the latitude-longitude grid.
    loni, lati = np.meshgrid(rg.lon[:,0],rg.lat[:,0])

    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    zlev = np.zeros(np.shape(loni)+(tsp,))
    if wind_vectors == True:
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

        zlev[:,:,t] = (z['value'][:,:,below,t]*(rg.Pressure[above,t] - Plev)\
            + z['value'][:,:,above,t]*(Plev - rg.Pressure[below,t]))\
            / (rg.Pressure[above,t]-rg.Pressure[below,t])

        if wind_vectors == True:
            Uii[:,:,t] = (rg.U[:,:,below,t]*(rg.Pressure[above,t] - Plev)\
                + rg.U[:,:,above,t]*(Plev - rg.Pressure[below,t]))\
                / (rg.Pressure[above,t]-rg.Pressure[below,t])
            Vii[:,:,t] = (rg.V[:,:,below,t]*(rg.Pressure[above,t] - Plev)\
                + rg.V[:,:,above,t]*(Plev - rg.Pressure[below,t]))\
                / (rg.Pressure[above,t]-rg.Pressure[below,t])

    # Averaging in time
    if tsp > 1:
        zlevt = np.mean(zlev,axis=2)
        del zlev
        if wind_vectors == True:
            Uiii = np.mean(Uii,axis=2)
            Viii = np.mean(Vii,axis=2)
            del Uii, Vii
    else:
        zlevt = zlev[:,:,0]
        del zlev
        if wind_vectors == True:
            Uiii = Uii[:,:,0]
            Viii = Vii[:,:,0]
            del Uii, Vii


    #################
    # Create Figure #
    #################

    lonp = rg.lon[:,0]
    latp = rg.lat[:,0]

    if isinstance(axis,axes.SubplotBase):
        C = axis.contourf(lonp,latp,zlevt,50,cmap=z['cmap'])
        ax = axis
    elif axis == False:
        C = plt.contourf(lonp,latp,zlevt,50,cmap=z['cmap'])
        ax = plt.gca()
    else:
        raise IOError("'axis = {}' but {} is not an axes.SubplotBase instance".format(axis,axis))

    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    ax.set_ylabel('Latitude (deg)')
    ax.set_xlabel('Longitude (deg)')

    ax.set_title('Time = %#.3f-%#.3f days, Plev = %#.3f bar'%(output.time[0],output.time[-1],Plev/1e5))

    if wind_vectors == True:
        d_z = np.shape(Uiii)
        spacing = np.int(np.shape(Uiii)[0]/10)
        U = Uiii[::spacing,::spacing].ravel()
        V = Viii[::spacing,::spacing].ravel()
        lonq = loni[::spacing,::spacing].ravel()
        latq = lati[::spacing,::spacing].ravel()
        del Uiii, Viii
        plt.quiver(lonq,latq,U,V,color='0.5')

    clb = plt.colorbar(C)
    clb.set_label(z['label'])

    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    if save == True:
        plt.savefig(input.resultsf+'/figures/%s_lev%#.3fmbar_i%d_l%d.pdf'%(z['name'],Plev/100,output.ntsi,output.nts))
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

def dFunc_dZ(f, alti, nv, dz):
    dfdz = np.zeros_like(f)
    # dz = alti[0,0,1]-alti[0,0,0]

    # Over z (vertical direction)
    for i in np.arange(nv):
        if i == 0:
            #lower edge, first order gradient
            dfdz[:,:,i] = (f[:,:,i+1]-f[:,:,i])/dz
        elif i == 1 or i == (nv-2):
            #second order gradient at penultimate layers
            dfdz[:,:,i] = (f[:,:,i+1]-f[:,:,i-1])/(2*dz)
        elif i == (nv-1):
            #highest edge, first order gradient
            dfdz[:,:,i] = (f[:,:,i]-f[:,:,i-1])/dz
        else:
            #five point stencil
            dfdz[:,:,i] = (-f[:,:,i+2]+8*f[:,:,i+1]\
                                -8*f[:,:,i-1]+f[:,:,i-2])/(12*dz)

    return dfdz

def dFunc_dLat(f, dLat):
    dfdlat = np.zeros_like(f)

    # Over lat (vertical direction)
    # for t in np.arange(np.shape(f)[-1]):
    for i in np.arange(np.shape(f)[0]):
        if i == 0:
            #lower edge, first order gradient
            dfdlat[i,:,:] = (f[i+1,:,:]-f[i,:,:])/dLat
        elif i == 1 or i == (np.shape(f)[0]-2):
            #second order gradient at penultimate layers
            dfdlat[i,:,:] = (f[i+1,:,:]-f[i-1,:,:])/(2*dLat)
        elif i == (np.shape(f)[0]-1):
            #highest edge, first order gradient
            dfdlat[i,:,:] = (f[i,:,:]-f[i-1,:,:])/dLat
        else:
            #five point stencil
            dfdlat[i,:,:] = (-f[i+2,:,:]+8*f[i+1,:,:]\
                                -8*f[i-1,:,:]+f[i-2,:,:])/(12*dLat)
    return dfdlat

def dFunc_dLon(f, dLon):
    dfdlon = np.zeros_like(f)

    # Over lon (vertical direction)
    for i in np.arange(np.shape(f)[1]):
        if i == 0:
            #lower edge, first order gradient
            dfdlon[:,i,:] = (f[:,i+1,:]-f[:,i,:])/dLon
        elif i == 1 or i == (np.shape(f)[1]-2):
            #second order gradient at penultimate layers
            dfdlon[:,i,:] = (f[:,i+1,:]-f[:,i-1,:])/(2*dLon)
        elif i == (np.shape(f)[1]-1):
            #highest edge, first order gradient
            dfdlon[:,i,:] = (f[:,i,:]-f[:,i-1,:])/dLon
        else:
            #five point stencil
            dfdlon[:,i,:] = (-f[:,i+2,:]+8*f[:,i+1,:]\
                                -8*f[:,i-1,:]+f[:,i-2,:])/(12*dLon)
    return dfdlon

def CurlF(fr,flat,flon,lati,alti,res_deg,nv,A,dz):
    curlFz = np.zeros_like(fr)
    curlFlat = np.zeros_like(fr)
    curlFlon = np.zeros_like(fr)

    curlFz = -1.0*(dFunc_dLat(flon*np.cos(lati),res_deg)-dFunc_dLon(flat,res_deg))
    curlFz /= (np.cos(lati)*(A+alti))

    curlFlat = -1.0*(dFunc_dLon(fr,res_deg)/np.cos(lati)-dFunc_dZ((A+alti)*flon,alti,nv,dz))
    curlFlat /= (A+alti)

    curlFlon = -1.0*(dFunc_dZ((A+alti)*flat,alti,nv,dz)-dFunc_dLat(fr,res_deg))
    curlFlon /= (A+alti)

    return curlFz, curlFlat, curlFlon

def calc_moc_streamf(grid,output,input,lons,lats,Pref,t_ind,fileh5,comp=4,pressure_vert=True):
    nlev = len(grid.Altitude)
    res_deg = 4.0/2**(input.glevel-4)
    #figure out pressure grid
    # pmin = np.min(output.Pressure)
    # pscale = 'log'
    # if pscale == 'log':
    #     sigmaref = np.logspace(np.log10(input.P_Ref),np.log10(pmin),nlev)/input.P_Ref
    # elif pscale == 'lin':
    #     sigmaref = np.linspace(input.P_Ref,pmin,nlev)/input.P_Ref
    # else:
    #     raise IOError('invalid pressure scale entered! use "lin" or "log"')
    #
    d_sig = np.size(Pref)
    # Pref = input.P_Ref*sigmaref[:,0]
    lat_range_tmp = np.arange(-90,90+res_deg,res_deg)
    # recenter so that there are an even number of latitude points
    lat_range = (lat_range_tmp[:-1]+lat_range_tmp[1:])/2
    lon_range = np.arange(0,360,res_deg)
    loni, lati = np.meshgrid(lon_range,lat_range)
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    rv_icoh = (output.Mh[0,:,:,t_ind]*(-np.sin(grid.lat[:,None])*np.cos(grid.lon[:,None])) + \
                output.Mh[1,:,:,t_ind]*(-np.sin(grid.lat[:,None])*np.sin(grid.lon[:,None])) + \
                output.Mh[2,:,:,t_ind]*np.cos(grid.lat[:,None]))

    # Pref = input.P_Ref*sigmaref[:,0]
    SF_llp = np.zeros((d_lon[0],d_lon[1],d_sig))
    # for t in np.arange(tsp):
    SF_icop = np.zeros((grid.point_num,d_sig))
    for inum in np.arange(grid.point_num):
        SF_icoh = np.zeros(input.vlevel)
        for ilev in np.arange(input.vlevel):
            if ilev == 0:
                z_tmp = np.linspace(0,grid.Altitude[ilev],10)
            else:
                z_tmp = np.linspace(grid.Altitude[ilev-1],grid.Altitude[ilev],10)
            rv_tmp = interp.pchip_interpolate(grid.Altitude,rv_icoh[inum,:],z_tmp)
            SF_icoh[ilev] = 2*np.pi*np.trapz((input.A+z_tmp)*rv_tmp,x=z_tmp)*np.cos(grid.lat[inum])
            if ilev > 0:
                SF_icoh[ilev] += SF_icoh[ilev-1]

        SF_icop[inum,:] = interp.pchip_interpolate(output.Pressure[inum,:,t_ind][::-1],SF_icoh[::-1],Pref[::-1])[::-1]

    for lev in np.arange(d_sig):
        SF_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,SF_icop[:,lev],(loni,lati),method='nearest')

    openh5 = h5py.File(fileh5,"r+")
    print('Adding streamfunction to regrid file...')
    stream = openh5.create_dataset("streamf",data=SF_llp,compression='gzip',compression_opts=comp)
    openh5.close()

def streamf_moc_plot(input,grid,output,rg,sigmaref,save=True,axis=False,wind_vectors=False):
    # special plotting function for the mass streamfunction

    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)
    tsp = output.nts-output.ntsi+1

    if tsp > 1:
        Vavgl = np.mean(rg.V,axis=1)
        Vavglt = np.mean(Vavgl,axis=2)
        if wind_vectors == True:
            Wavgl = np.mean(rg.W,axis=1)
            Wavglt = np.mean(Wavgl,axis=2)
    else:
        Vavglt = np.mean(rg.V[:,:,:,0],axis=1)
        if wind_vectors == True:
            Wavglt = np.mean(rg.V[:,:,:,0],axis=1)

    sf = np.zeros_like(Vavglt)
    arg = 2*np.pi*input.A*np.cos(rg.lat[:,0][:,None]*np.pi/180)/input.Gravit*Vavglt
    for ilat in np.arange(np.shape(Vavglt)[0]):
        for ilev in np.arange(np.shape(Vavglt)[1]):
            if ilev == 0:
                sf[ilat,-1] = arg[ilat,-1]*rg.Pressure[-1,0]
            else:
                sf[ilat,-(ilev+1)] = np.trapz(arg[ilat,-1:-(ilev+2):-1],x=rg.Pressure[-1:-(ilev+2):-1][:,0])

    # need to set desired pressure range (major PITA!)
    prange = np.where(np.logical_and(rg.Pressure>=np.min(Pref),rg.Pressure<=np.max(Pref)))

    # Contour plot
    if isinstance(axis,axes.SubplotBase):
        C = axis.contourf(rg.lat[:,0],rg.Pressure[prange[0],0]/1e5,sf[:,prange[0]].T,40,cmap = 'viridis')
        ax = axis
    elif axis == False:
        C = plt.contourf(rg.lat[:,0],rg.Pressure[prange[0],0]/1e5,sf[:,prange[0]].T,40,cmap = 'viridis')
        ax = plt.gca()
    else:
        raise IOError("'axis = {}' but {} is not an axes.SubplotBase instance".format(axis,axis))

    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0

    if wind_vectors == True:
        vspacing = np.int(np.shape(rg.lat)[0]/10)
        wspacing = np.int(np.shape(rg.Pressure)[0]/10)
        Vlt = Vavglt[:,prange[0]]
        Wlt = Wavglt[:,prange[0]]
        Vq = Vlt[::vspacing,::wspacing].ravel()
        Wq = Wlt[::vspacing,::wspacing].ravel()
        #preq = rg.Pressure[:,0][::spacing,::spacing].ravel()
        #latq = lati[::spacing,::spacing].ravel()
        latq, preq = np.meshgrid(rg.lat[::vspacing,0],rg.Pressure[::wspacing,0][prange[0]])
        del Vlt, Wlt
        plt.quiver(latq.ravel(),preq.ravel()/1e5,Vq,Wq,color='0.5')

    ax.invert_yaxis()
    ax.set_yscale("log")
    c2 = ax.contour(rg.lat[:,0],rg.Pressure[:,0]/1e5,sf.T,levels=[0.0],colors='w',linewidths=1)
    clb = plt.colorbar(C)
    clb.set_label(r'Eulerian streamfunction (kg s$^{-1}$)')
    if np.max(Pref)/np.min(Pref) > 100:
        ax.set_yscale("log")
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Pressure (bar)')
    ax.plot(rg.lat[:,0],np.zeros_like(rg.lat[:,0])+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
    if np.min(rg.Pressure[prange[0],0]) < np.max(output.Pressure[:,grid.nv-1,:]):
        ax.set_ylim(np.max(rg.Pressure[prange[0],0])/1e5,np.min(rg.Pressure[prange[0],0])/1e5)
    else:
        ax.set_ylim(np.max(rg.Pressure[prange[0],0])/1e5,np.max(output.Pressure[:,grid.nv-1,:])/1e5)
    ax.set_title('Time = %#.3f-%#.3f days, Lon = (0,360)'%(output.time[0],output.time[-1]))
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    if save == True:
        plt.savefig(input.resultsf+'/figures/streamf_ver_i%d_l%d.pdf'%(output.ntsi,output.nts))
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
    dz = np.abs(grid.Altitude[1] - grid.Altitude[0])
    dThdz = dFunc_dZ(Theta,alti,grid.nv,dz)

    dThdlat = dFunc_dLat(Theta,res_deg)/(input.A+alti)

    dThdlon = dFunc_dLon(Theta,res_deg)/(input.A+alti)/np.cos(lati)

    # Compute curl of wind speeds
    curlVz, curlVlat, curlVlon = CurlF(W,V,U,lati,alti,res_deg,grid.nv,input.A,dz)

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
    dz = np.abs(grid.Altitude[1] - grid.Altitude[0])

    curlVz, curlVlat, curlVlon = CurlF(W,V,U,lati,alti,res_deg,grid.nv,input.A,dz)

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

def profile(input,grid,output,z,stride=50):
    # Pref = input.P_Ref*sigmaref
    # d_sig = np.size(sigmaref)

    tsp = output.nts-output.ntsi+1

    for column in np.arange(0,grid.point_num,stride):
        if tsp > 1:
            P = np.mean(output.Pressure[column,:,:],axis=2)
            x = np.mean(z['value'][column,:,:],axis=2)
        else:
            P = output.Pressure[column,:,0]
            x = z['value'][column,:,0]

        plt.semilogy(x,P/1e5,'k-',alpha= 0.5,lw=1)
        plt.plot(x[np.int(np.floor(grid.nv/2))],P[np.int(np.floor(grid.nv/2))]/100000,'r+',ms =5,alpha=0.5)
        plt.plot(x[np.int(np.floor(grid.nv*0.75))],P[np.int(np.floor(grid.nv*0.75))]/100000,'g+',ms =5,alpha=0.5)

    # plt.plot(Tad,P/100,'r--')
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure (bar)')
    plt.xlabel(z['label'])
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.savefig(input.resultsf+'/figures/%sPprofile_i%d_l%d.pdf'%(z['name'],output.ntsi,output.nts))
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
    Pref = input.P_Ref*sigmaref[:,0]
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
        if (output.ConvData == False).any():
            print('Calculating energy, mass, angular momentum...')
            CalcE_M_AM(input,grid,output,split)
        if (output.EntData == False).any():
            print('Calculating entropy...')
            CalcEntropy(input,grid,output,split)
        plots = ['global']

    else:
        CalcE_M_AM(input,grid,output,split)
        CalcEntropy(input,grid,output,split)
        plots = ['weather', 'deep']

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
