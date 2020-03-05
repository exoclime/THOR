import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.axes as axes
import matplotlib.cm as cm
import scipy.interpolate as interp
import scipy.ndimage as ndimage
import os
import h5py
import time
import subprocess as spr
import pyshtools as chairs
import pdb
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import pycuda.gpuarray as gpuarray

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
        if 'radiative_transfer' in openh5.keys():
            self.RT = openh5['radiative_transfer'][0]
            self.Tstar = openh5['Tstar'][...]
            self.planet_star_dist = openh5['planet_star_dist'][...]
            self.radius_star = openh5['radius_star'][...]
            if 'diff_fac' in openh5.keys(): #called diff_fac in old versions
                self.diff_ang = openh5['diff_fac'][...]
            else:
                self.diff_ang = openh5['diff_ang'][...]
            if 'Tint' in openh5.keys():
                self.Tint = openh5['Tint'][...]
            self.albedo = openh5['albedo'][...]
            self.tausw = openh5['tausw'][...]
            self.taulw = openh5['taulw'][...]
            if 'surface' in openh5.keys():
                self.surface = openh5['surface'][0]
            else:
                self.surface = 0
        else:
            self.RT = 0
        if 'vulcan' in openh5.keys():
            self.chemistry = openh5['vulcan'][...]
        if 'chemistry' in openh5.keys():
            self.chemistry = openh5['chemistry'][...]

        openh5.close()

def LatLong2Cart(lat,lon):
    x = np.cos(lat)*np.cos(lon)
    y = np.cos(lat)*np.sin(lon)
    z = np.sin(lat)
    return x, y, z

def Cart2LatLong(x,y,z):
    lon = np.arctan2(y,x)
    lat = np.arctan2(z,np.sqrt(x**2+y**2))
    return lat, lon

def Rot_y(theta,x,y,z):
    xnew = np.cos(theta)*x + np.sin(theta)*z
    znew = -np.sin(theta)*x + np.cos(theta)*z
    return xnew, y, znew

def Rot_z(phi,x,y,z):
    xnew = np.cos(phi)*x - np.sin(phi)*y
    ynew = np.sin(phi)*x + np.cos(phi)*y
    return xnew, ynew, z

class grid:
    def __init__(self,resultsf,simID,rotation=False,theta_y=0,theta_z=0):
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
        if rotation == True:
            xtmp, ytmp, ztmp = LatLong2Cart(self.lat, self.lon)
            xtmp1, ytmp1, ztmp1 = Rot_z(theta_z,xtmp,ytmp,ztmp)
            xtmp2, ytmp2, ztmp2 = Rot_y(theta_y,xtmp1,ytmp1,ztmp1)
            self.lat, self.lon = Cart2LatLong(xtmp2,ytmp2,ztmp2)
            self.lon = self.lon%(2*np.pi)

class output:
    def __init__(self,resultsf,simID,ntsi,nts,input,grid,stride=1):
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
        self.Tsurface = np.zeros((grid.point_num,nts-ntsi+1))

        self.ch4 = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.co = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.h2o = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.co2 = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.nh3 = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))

        self.tau_sw = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.tau_lw = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.flw_up = np.zeros((grid.point_num,grid.nvi,nts-ntsi+1))
        self.flw_dn = np.zeros((grid.point_num,grid.nvi,nts-ntsi+1))
        self.fsw_dn = np.zeros((grid.point_num,grid.nvi,nts-ntsi+1))

        self.Rd = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.Cp = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))

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

            if 'Rd' in openh5.keys():
                Rdi = openh5['Rd'][...]
                Cpi = openh5['Cp'][...]
            else:
                Rdi = np.zeros_like(Rhoi)+input.Rd
                Cpi = np.zeros_like(Rhoi)+input.Cp

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
                if 'fnet_up' in openh5.keys(): #old version of LW output
                    fupi = openh5['fnet_up'][...]
                    fdni = openh5['fnet_dn'][...]
                    sdni = np.zeros_like(fdni)
                else:
                    fupi = openh5['flw_up'][...]
                    fdni = openh5['flw_dn'][...]
                    sdni = openh5['fsw_dn'][...]
            if 'Tsurface' in openh5.keys():
                self.Tsurface[:,t-ntsi+1] = openh5['Tsurface'][...]
            openh5.close()

            self.Rho[:,:,t-ntsi+1] = np.reshape(Rhoi,(grid.point_num,grid.nv))
            self.Pressure[:,:,t-ntsi+1] = np.reshape(Pressurei,(grid.point_num,grid.nv))
            self.Mh[0,:,:,t-ntsi+1] = np.reshape(Mhi[::3],(grid.point_num,grid.nv))
            self.Mh[1,:,:,t-ntsi+1] = np.reshape(Mhi[1::3],(grid.point_num,grid.nv))
            self.Mh[2,:,:,t-ntsi+1] = np.reshape(Mhi[2::3],(grid.point_num,grid.nv))
            self.Wh[:,:,t-ntsi+1] = np.reshape(Whi,(grid.point_num,grid.nvi))
            self.Rd[:,:,t-ntsi+1] = np.reshape(Rdi,(grid.point_num,grid.nv))
            self.Cp[:,:,t-ntsi+1] = np.reshape(Cpi,(grid.point_num,grid.nv))

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
                self.flw_dn[:,:,t-ntsi+1] = np.reshape(fdni,(grid.point_num,grid.nvi))
                self.flw_up[:,:,t-ntsi+1] = np.reshape(fupi,(grid.point_num,grid.nvi))
                self.fsw_dn[:,:,t-ntsi+1] = np.reshape(sdni,(grid.point_num,grid.nvi))

class rg_out:
    def __init__(self,resultsf,simID,ntsi,nts,input,output,grid,pressure_vert=True,pgrid_ref='auto'):

        if input.RT == 1:
            surf = 0
            if input.surface:
                surf = 1
        else:
            surf = 0

        chem = 0
        if hasattr(input,"chemistry"):
            if input.chemistry == 1:
                chem = 1

        # Read model results
        for t in np.arange(ntsi-1,nts):
            if pressure_vert == True:
                fileh5 = resultsf+'/regrid_'+simID+'_'+np.str(t+1)
            else:
                fileh5 = resultsf+'/regrid_height_'+simID+'_'+np.str(t+1)
            fileh5 += '.h5'
            if os.path.exists(fileh5):
                openh5 = h5py.File(fileh5)
            else:
                print(fileh5+' not found, regridding now with default settings...')
                regrid(resultsf,simID,ntsi,nts,pgrid_ref=pgrid_ref)
                openh5 = h5py.File(fileh5)

            Rhoi = openh5['Rho'][...]
            Tempi = openh5['Temperature'][...]
            Ui = openh5['U'][...]
            Vi = openh5['V'][...]
            Wi = openh5['W'][...]
            if pressure_vert == True:
                Prei = openh5['Pressure'][...]
            else:
                Alti = openh5['Altitude'][...]
                Prei = np.zeros_like(Alti)
            lati = openh5['Latitude'][...]
            loni = openh5['Longitude'][...]

            if input.RT == 1:
                tau_swi = openh5['tau_sw'][...]
                tau_lwi = openh5['tau_lw'][...]
                if 'fnet_up' in openh5.keys():
                    flw_upi = openh5['fnet_up'][...]
                    flw_dni = openh5['fnet_dn'][...]
                else:
                    flw_upi = openh5['flw_up'][...]
                    flw_dni = openh5['flw_dn'][...]
                insoli = openh5['insol'][...]
                if 'Tsurface' in openh5.keys():
                    Tsurfi = openh5['Tsurface'][...]
                else:
                    surf = 0

            if chem == 1:
                ch4i = openh5['ch4'][...]
                coi = openh5['co'][...]
                h2oi = openh5['h2o'][...]
                co2i = openh5['co2'][...]
                nh3i = openh5['nh3'][...]

            PVi = openh5['PV'][...]
            if 'RV' in openh5.keys():
                RVi = openh5['RV'][...][0]
            else:
                RVi = openh5['RVw'][...]

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

                if input.RT == 1:
                    self.tau_sw = np.zeros(np.shape(tau_swi)+(nts-ntsi+1,))
                    self.tau_lw = np.zeros(np.shape(tau_lwi)+(nts-ntsi+1,))
                    self.flw_up = np.zeros(np.shape(flw_upi)+(nts-ntsi+1,))
                    self.flw_dn = np.zeros(np.shape(flw_dni)+(nts-ntsi+1,))
                    self.insol = np.zeros(np.shape(insoli)+(nts-ntsi+1,))
                    if surf:
                        self.Tsurface = np.zeros(np.shape(Tsurfi)+(nts-ntsi+1,))
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
            if pressure_vert == True:
                self.Pressure[:,t-ntsi+1] = Prei
                if t > ntsi-1:
                    if (Prei == self.Pressure[:,0]).all() == False:
                        print('Warning: Pressure grid in file %d does not match file %d'%(t,ntsi))
            else:
                self.Altitude[:,t-ntsi+1] = Alti
            self.lat[:,t-ntsi+1] = lati
            self.lon[:,t-ntsi+1] = loni
            self.PV[:,:,:,t-ntsi+1] = PVi
            self.RV[:,:,:,t-ntsi+1] = RVi

            if input.RT == 1:
                self.tau_sw[:,:,:,t-ntsi+1] = tau_swi
                self.tau_lw[:,:,:,t-ntsi+1] = tau_lwi
                self.flw_up[:,:,:,t-ntsi+1] = flw_upi
                self.flw_dn[:,:,:,t-ntsi+1] = flw_dni
                self.insol[:,:,t-ntsi+1] = insoli
                if surf:
                    self.Tsurface[:,:,t-ntsi+1] = Tsurfi
            if chem == 1:
                self.ch4[:,:,:,t-ntsi+1] = ch4i
                self.co[:,:,:,t-ntsi+1] = coi
                self.h2o[:,:,:,t-ntsi+1] = h2oi
                self.co2[:,:,:,t-ntsi+1] = co2i
                self.nh3[:,:,:,t-ntsi+1] = nh3i


class GetOutput:
    def __init__(self,resultsf,simID,ntsi,nts,stride=1,openrg=0,pressure_vert=True,rotation=False,theta_y=0,theta_z=0,pgrid_ref='auto'):
        self.input = input(resultsf,simID)
        self.grid = grid(resultsf,simID,rotation=rotation,theta_y=theta_y,theta_z=theta_z)
        self.output = output(resultsf,simID,ntsi,nts,self.input,self.grid,stride=stride)
        if openrg == 1:
            self.rg = rg_out(resultsf,simID,ntsi,nts,self.input,self.output,self.grid,
                                pressure_vert=pressure_vert,pgrid_ref=pgrid_ref)

def define_Pgrid(resultsf,simID,ntsi,nts,stride,overwrite=False):
    # first we need the grid size
    fileh5 = resultsf+'/esp_output_grid_'+simID+'.h5'
    if os.path.exists(fileh5):
        openh5 = h5py.File(fileh5)
    else:
        raise IOError(fileh5+' not found!')
    nv = np.int(openh5['nv'][0])
    point_num = np.int(openh5['point_num'][0])
    Altitude = openh5['Altitude'][...]
    Altitudeh = openh5['Altitudeh'][...]
    openh5.close()

    # we also need to know if we included a surface
    fileh5 = resultsf+'/esp_output_planet_'+simID+'.h5'
    if os.path.exists(fileh5):
        openh5 = h5py.File(fileh5)
    else:
        raise IOError(fileh5+' not found!')
    if 'surface' in openh5.keys():
        surf = openh5['surface'][0]
    else:
        surf = 0
    openh5.close()

    # now we'll loop over all the files to get the pressure_mean
    num_out = np.int((nts-ntsi)/stride) + 1
    pressure_mean = np.zeros((point_num,nv,num_out))
    for i in np.arange(num_out):
        t = ntsi + stride*i
        fileh5 = resultsf+'/esp_output_'+simID+'_'+np.str(t)+'.h5'
        if os.path.exists(fileh5):
            openh5 = h5py.File(fileh5)
        else:
            raise IOError(fileh5+' not found!')

        pressure_mean[:,:,i] = np.reshape(openh5['Pressure_mean'][...],(point_num,nv))
        openh5.close()

    pgrid = np.mean(np.mean(pressure_mean,axis=0),axis=1)
    if surf == 1:
        extrap_low = (Altitudeh[0]-Altitude[1])/(Altitude[0]-Altitude[1])
        Psurf = pressure_mean[:,1,:]+extrap_low*(pressure_mean[:,0,:]-pressure_mean[:,1,:])
        Psurf_mean = np.mean(np.mean(Psurf,axis=0),axis=0)
        pgrid = np.concatenate((np.array([Psurf_mean]),pgrid))

    pfile = 'pgrid_%d_%d_%d.txt'%(ntsi,nts,stride)
    if not os.path.exists(resultsf+'/'+pfile) or overwrite==True:
        f = open(resultsf+'/'+pfile,'w')
        for j in np.arange(len(pgrid)):
            f.write('%d %#.6e\n'%(j,pgrid[j]))
        f.close()
    else:
        raise IOError(pfile+' already exists in this directory!')

regrid_tools = SourceModule("""
    __device__ double dot_product(double *a, double *b) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }

    __global__ void find_nearest(int *near, double *a, double *a_ico, int num_ll, int num_ico) {
        int id = blockIdx.x * blockDim.x + threadIdx.x;

        if (id < num_ll) {
            double dot_max = 0.0, dot_new;
            int closest = -1;
            for (int i_ico = 0; i_ico < num_ico; i_ico++) {
                dot_new = dot_product(&a[id*3],&a_ico[i_ico*3]);
                if (dot_new > dot_max) {
                    dot_max = dot_new;
                    closest = i_ico;
                }
            }
            near[id] = closest;
        }
    }

    __device__ double determinant(double *a, double *b, double *c) {
        return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]);
    }

    __device__ void weights_tuv(double *a, double *b, double *c, double *d,
                                double *t, double *u, double *v) {
        // Simple code-up of Moller-Trumbore algorithm
        // Ref: [Fast, Minimum Storage Ray/Triangle Intersection.
        //       Möller & Trumbore. Journal of Graphics Tools, 1997.]

        double *OD = new double[3]();
        double *OA = new double[3]();
        double *BA = new double[3]();
        double *CA = new double[3]();

        for (int i = 0; i < 3; i++){
            OD[i] = -d[i];
            OA[i] = -a[i];
            BA[i] = b[i] - a[i];
            CA[i] = c[i] - a[i];
        }

        double det0, det1, det2, det3;

        det0 = determinant(OD,BA,CA);
        det1 = determinant(OA,BA,CA);
        det2 = determinant(OD,OA,CA);
        det3 = determinant(OD,BA,OA);

        *t = det1/det0;
        *u = det2/det0;
        *v = det3/det0;

        delete[] OD;
        delete[] OA;
        delete[] BA;
        delete[] CA;

    }

    __global__ void calc_weights(double *weight3, int *near3, double *a, double *a_ico,
                                 int *pntloc, int num_ll, int num_ico) {
        int id = blockIdx.x * blockDim.x + threadIdx.x;

        if (id < num_ll) {
            int i_ico = near3[id*3+0], j = 0, jp1, ii1, ii2;
            double t, u, v;
            bool correct_triangle = false;
            while (correct_triangle == false){
                jp1 = (j+1)%6;
                ii1 = pntloc[i_ico*6+j];
                ii2 = pntloc[i_ico*6+jp1];
                weights_tuv(&a_ico[i_ico*3],&a_ico[ii1*3],&a_ico[ii2*3],&a[id*3],&t,&u,&v);

                if ((u<0) || (u>1)) {
                    correct_triangle = false;
                } else if ((v<0) || (u+v>1)) {
                    correct_triangle = false;
                } else {
                    correct_triangle = true;
                    near3[id*3+1] = ii1;
                    near3[id*3+2] = ii2;
                    weight3[id*3+0] = 1.0-u-v;
                    weight3[id*3+1] = u;
                    weight3[id*3+2] = v;
                }
                j += 1;
            }
        }
    }

    __global__ void vert_lin_interp(double *x, double *y, double *xnew,
                                    double *ynew, int nxnew, int point_num, int nv,
                                    int *x_non_mono_check, int *xnew_non_mono_check) {
        // linear interpolation for vertical columns (pressure or height)
        // x and xnew must be monotonically increasing
        // in the case that they are not, the mono_check array entries will = 1
        int id = blockIdx.x * blockDim.x + threadIdx.x;

        if (id < point_num) {
            for (int lev = 1; lev < nv; lev++) {
                //check for monotonicity of column (x -> increasing order)
                if (x[id*nv+lev] < x[id*nv+lev-1]) {
                    x_non_mono_check[id] = 1;
                }
            }

            for (int ilev = 1; ilev < nxnew; ilev++) {
                //check for monotonicity of new vert grid
                if (xnew[ilev] < xnew[ilev-1]) {
                    xnew_non_mono_check[0] = 1;
                }
            }

            if ((x_non_mono_check[id]==0) && (xnew_non_mono_check[0]==0)){
                //monotonicity ok, do the interpolation
                for (int ilev = 0; ilev < nxnew; ilev++) {
                    int lev;
                    if (xnew[ilev] < x[id*nv+0]) {
                        //extrapolate below
                        lev = 0;
                    } else if (xnew[ilev] > x[id*nv+nv-1]) {
                        //extrapolate above
                        lev = nv - 2;
                    } else {
                        //interpolation
                        lev = 0;
                        while (xnew[ilev] > x[id*nv+lev+1]) {
                            lev += 1;
                        }
                    }
                    ynew[id*nxnew+ilev] = y[id*nv+lev] + (xnew[ilev]-x[id*nv+lev])*
                                        (y[id*nv+lev+1]-y[id*nv+lev])/
                                        (x[id*nv+lev+1]-x[id*nv+lev]);
                }
            }
        }
    }

""")

def create_rg_map(resultsf,simID,rotation=False,theta_z=0,theta_y=0):
    outall = GetOutput(resultsf,simID,0,0,rotation=rotation,theta_z=theta_z,theta_y=theta_y)
    input = outall.input
    grid = outall.grid

    #vector positions of ico grid points
    v_ico = np.vstack([np.cos(grid.lon)*np.cos(grid.lat),\
                       np.sin(grid.lon)*np.cos(grid.lat),\
                       np.sin(grid.lat)]).T

    #set horizontal angular resolution of latlon grid roughly same as ico grid
    ang_res = (4.0/2**(input.glevel-4))*np.pi/180

    #1D lat and lon arrays
    lat_range = np.arange(-np.pi/2+ang_res/2,np.pi/2,ang_res)
    lon_range = np.arange(0,2*np.pi,ang_res)
    loni, lati = np.meshgrid(lon_range,lat_range)
    num_ll = len(lat_range)*len(lon_range)

    #vector positions of latitude-longitude points
    v_ll = np.vstack([(np.cos(loni)*np.cos(lati)).ravel(),\
                      (np.sin(loni)*np.cos(lati)).ravel(),\
                      (np.sin(lati)).ravel()]).T

    #arrays for nearest neighbors and weights
    near3 = np.zeros((num_ll,3),dtype=np.int32)
    weight3 = np.zeros((num_ll,3))
    near = np.zeros(num_ll,dtype=np.int32)

    #cuda functions and sizes
    find_nearest = regrid_tools.get_function("find_nearest")
    calc_weights =  regrid_tools.get_function("calc_weights")
    gridGPU = (np.int(np.floor(num_ll/256))+1,1,1)
    blockGPU = (256,1,1)

    #nearest neighbor search
    find_nearest(cuda.Out(near),cuda.In(v_ll.ravel()),cuda.In(v_ico.ravel()),\
             np.int32(num_ll),np.int32(grid.point_num),\
             block=blockGPU,grid=gridGPU)
    near3[:,0] = near

    #next 2 neighbors and weights from Moller-Trumbore algorithm
    calc_weights(cuda.Out(weight3.ravel()),cuda.InOut(near3.ravel()),\
            cuda.In(v_ll.ravel()),cuda.In(v_ico.ravel()), \
            cuda.In(grid.pntloc.ravel()),np.int32(num_ll),\
            np.int32(grid.point_num),block=blockGPU,grid=gridGPU)

    #save into numpy zip file
    np.savez(resultsf+'/regrid_map.npz',lon=lon_range,lat=lat_range,\
             near3=near3,weight3=weight3)

def vertical_regrid_field(source_array,nv,x,xnew):
    #handles set up and running of vertical interpolation in pycuda
    vert_lin_interp = regrid_tools.get_function("vert_lin_interp")
    x_non_mono_check = np.zeros(np.shape(source_array)[:-1],dtype=np.int32)
    xnew_non_mono_check = np.zeros(1,dtype=np.int32)
    gridgpu = (np.int(np.floor(len(x_non_mono_check.ravel())/256))+1,1,1)
    blockgpu = (256,1,1)

    y = np.ascontiguousarray(source_array.ravel())
    ynew = np.zeros(np.shape(source_array)[:-1]+(len(xnew),))
    vert_lin_interp(cuda.In(x),cuda.In(y),cuda.In(xnew),\
                    cuda.Out(ynew.ravel()),np.int32(len(xnew)),\
                    np.int32(len(x_non_mono_check.ravel())),np.int32(nv),\
                    cuda.InOut(x_non_mono_check.ravel()),\
                    cuda.InOut(xnew_non_mono_check),grid=gridgpu,block=blockgpu)
    if (x_non_mono_check==1).any():
        raise ValueError('Pressure interpolation failed! Pressure not monotonically decreasing with height')
    if (xnew_non_mono_check[0] == 1):
        raise ValueError('Pressure interpolation failed! Destination pgrid not monotonically decreasing with height')
    return ynew

def regrid(resultsf,simID,ntsi,nts,pgrid_ref='auto',overwrite=False,comp=4,
            rotation=False,theta_z=0,theta_y=0,mask_surf=True):
    #New regridding process, smoother and faster than previous one
    #For each point on destination grid (lat-lon), find 3 points (near3) on the ico grid
    #that define a triangle containing the destination point. Then use the
    #Moller-Trumbore algorithm to calculate the weights (weight3) that interpolate the 3
    #points to the destination point. PyCuda is used to accelerate the process.
    #Interpolating each field is then a vectorized process using Numpy, taking
    #field_on_ll_grid = weight3 * field_on_ico_grid[near3]
    #Remapping from height to pressure is done on all columns via pycuda

    if not os.path.exists(resultsf+'/regrid_map.npz'):
        #if numpy zip file containing weights d.n.e., then create it
        create_rg_map(resultsf,simID,rotation=rotation,theta_z=theta_z,theta_y=theta_y)

    #open numpy zip file containing weights
    rg_map = np.load(resultsf+'/regrid_map.npz')
    lon_range = rg_map['lon']
    lat_range = rg_map['lat']
    near3 = rg_map['near3']
    weight3 = rg_map['weight3']

    #open first file to get some basic properties
    outall = GetOutput(resultsf,simID,ntsi,ntsi,rotation=rotation,theta_z=theta_z,theta_y=theta_y)
    input = outall.input
    grid = outall.grid
    output = outall.output

    print('Regrid data in folder '+resultsf+'...\n')

    #handling of vertical coordinate
    if pgrid_ref == 'auto':
        try:
            files_found = spr.check_output('ls '+resultsf+'/pgrid_*.txt',shell=True).split()
        except:
            raise IOError('No pgrid file found in "%s/", please run "pgrid -i <first> -l <last> %s"'%(resultsf,resultsf))
        if len(files_found) > 1:
            raise IOError('Multiple pgrid files found in "%s/", please specify with -pgrid flag'%resultsf)
        else:
            pgrid_file = files_found[0].decode()
    else:
        if os.path.exists(resultsf+'/'+pgrid_ref):
            pgrid_file = resultsf+'/'+ pgrid_ref
        else:
            raise IOError('pgrid file %s not found!'%(resultsf+'/'+pgrid_ref))
    print('Vertical coordinate = pressure from file %s'%pgrid_file)
    Pref = np.loadtxt(pgrid_file,unpack=True,usecols=1)

    #destination grid and set some dimensions
    loni, lati = np.meshgrid(lon_range,lat_range)
    d_lon = np.shape(loni)
    tsp = nts-ntsi+1
    d_sig = np.size(Pref)

    #handling of module properties
    if input.RT == 1:#gray radiative transfer module
        surf = 0
        if input.surface:
            surf = 1
            extrap_low = (grid.Altitudeh[0]-grid.Altitude[1])/(grid.Altitude[0]-grid.Altitude[1])
            Psurf = output.Pressure[:,1,:]+extrap_low*(output.Pressure[:,0,:]-output.Pressure[:,1,:])
    else:
        surf = 0

    chem = 0        #chemistry module
    if hasattr(input,"chemistry"):
        if input.chemistry == 1:
            chem = 1

    #begin regrid loop over times
    for t in np.arange(ntsi,nts+1):
        #check for existing h5 files
        proceed = 0
        fileh5p = resultsf+'/regrid_'+simID+'_'+np.str(t)

        fileh5h = resultsf+'/regrid_height_'+simID+'_'+np.str(t)
        if rotation == True:        #tell the user they asked for a rotated grid
            print('Applied rotation (theta_z,theta_y) = (%f,%f) to grid\n'\
                  %(theta_z*180/np.pi,theta_y*180/np.pi))
        fileh5p += '.h5'; fileh5h += '.h5'
        if os.path.exists(fileh5p) or os.path.exists(fileh5h):  #regrid files exist
            if overwrite == True:    #overwrite existing files
                proceed = 1
            else:                    #skip existing files
                print(fileh5p+'or '+fileh5h+' already present! Skipping time = %d'%t)
        else:                        #no existing regrid files, all clear
            proceed = 1

        #begin regridding
        if proceed == 1:
            print('Regridding time = %d...'%t)

            #read in one file at a time to prevent mem overflow
            outall = GetOutput(resultsf,simID,t,t,rotation=rotation,theta_z=theta_z,theta_y=theta_y)
            output = outall.output

            #interpolate in z to middle of layers using this
            interpz = (grid.Altitude-grid.Altitudeh[:-1])/(grid.Altitudeh[1:]-grid.Altitudeh[:-1])

            #set up dictionary with source arrays (add new quantities here!)
            source = {'Temperature':(output.Pressure/(output.Rd*output.Rho))[:,:,0],
                      'W':(output.Wh[:,:-1,0] + (output.Wh[:,1:,0]-output.Wh[:,:-1,0])*interpz[None,:])/output.Rho[:,:,0],
                      'Rho':output.Rho[:,:,0],
                      'Mh':output.Mh[:,:,:,0],
                      'Pressure':output.Pressure[:,:,0],
                      'Rd':output.Rd[:,:,0],
                      'Cp':output.Cp[:,:,0]}

            if input.RT == 1:
                source['flw_up'] = output.flw_up[:,:-1,0]+(output.flw_up[:,1:,0]-output.flw_up[:,:-1,0])*interpz[None,:]
                source['flw_dn'] = output.flw_dn[:,:-1,0]+(output.flw_dn[:,1:,0]-output.flw_dn[:,:-1,0])*interpz[None,:]
                source['tau_sw'] = output.tau_sw[:,:,0]
                source['tau_lw'] = output.tau_lw[:,:,0]
                source['insol'] = output.Insol[:,0]
                if surf == 1:
                    source['Tsurface'] = output.Tsurface[:,0]
                    source['Psurf'] = Psurf[:,0]

            if chem == 1:
                source['ch4'] = output.ch4[:,:,0]
                source['co'] = output.co[:,:,0]
                source['h2o'] = output.h2o[:,:,0]
                source['co2'] = output.co2[:,:,0]
                source['nh3'] = output.nh3[:,:,0]

            #calculate zonal and meridional velocity (special step for Mh)
            source['U'] = (-source['Mh'][0]*np.sin(grid.lon[:,None])+\
                           source['Mh'][1]*np.cos(grid.lon[:,None]))/source['Rho']
            source['V'] = (-source['Mh'][0]*np.sin(grid.lat[:,None])*np.cos(grid.lon[:,None])\
                      -source['Mh'][1]*np.sin(grid.lat[:,None])*np.sin(grid.lon[:,None])\
                      +source['Mh'][2]*np.cos(grid.lat[:,None]))/source['Rho']

            #set up intermediate arrays (icogrid and pressure)
            interm = {}
            for key in source.keys():
                if key == 'Mh':
                    pass
                elif np.shape(source[key])==(grid.point_num,):
                    #2D field (e.g., insolation) -> not needed
                    interm[key] = np.zeros((d_lon[0],d_lon[1]))
                else:
                    interm[key] = np.zeros((d_lon[0],d_lon[1],grid.nv))

            for key in interm.keys():
                if np.shape(interm[key]) == (d_lon[0],d_lon[1]):
                    #2D field (e.g., insolation)
                    tmp = np.sum(weight3[:,:]*source[key][near3],axis=1)
                    interm[key][:,:] = tmp.reshape((d_lon[0],d_lon[1]))
                else:
                    tmp = np.sum(weight3[:,:,None]*source[key][near3],axis=1)
                    interm[key][:,:] = tmp.reshape((d_lon[0],d_lon[1],grid.nv))

            #dealing with RV and PV
            #potential temp (only for constant Rd and Cp, currently)
            pt = interm['Temperature']*(interm['Pressure']/input.P_Ref)**(-interm['Rd']/interm['Cp'])
            dptdz = np.gradient(pt,grid.Altitude,axis=2)
            # these are not very accurate at the 'edges'
            dptdlat = np.gradient(pt,lat_range,axis=0)/(input.A+grid.Altitude[None,None,:])
            dptdlon = np.gradient(pt,lon_range,axis=1)/(input.A+grid.Altitude[None,None,:])/np.cos(lat_range[:,None,None])

            #relative vorticity
            curlVz, curlVlat, curlVlon = CurlF(interm['W'],interm['V'],interm['U'],\
                                               lat_range,lon_range,grid.Altitude,input.A)
            interm['RVw'] = curlVz    #vertical component
            interm['RVv'] = curlVlat  #horizontal components
            interm['RVu'] = curlVlon  # these are probably small in most cases

            curlVlon += 2*input.Omega #add solid body rotation
            interm['PV'] = (curlVz*dptdz + curlVlat*dptdlat + curlVlon*dptdlon)/interm['Rho']

            #set up destination arrays (lat-lon and pressure/height)
            dest = {}
            for key in interm.keys():
                if key == 'Mh' or key == 'Pressure':
                    pass  #don't need these any further
                elif np.shape(interm[key])==(d_lon[0],d_lon[1]):
                    #2D field (e.g., insolation)
                    dest[key] = np.zeros((d_lon[0],d_lon[1]))
                else:
                    dest[key] = np.zeros((d_lon[0],d_lon[1],d_sig))

            #regrid to pressure coordinate
            #need to make sure the arrays are contiguous
            x = np.ascontiguousarray(interm['Pressure'][:,:,::-1].ravel())
            xnew = np.ascontiguousarray(Pref[::-1])

            for key in dest.keys():
                if np.shape(interm[key])==(d_lon[0],d_lon[1]):
                    #2D field (e.g., insolation)
                    dest[key] = interm[key]

                else:
                    dest[key][:,:,:] = vertical_regrid_field(interm[key][:,:,::-1],grid.nv,x,xnew)[:,:,::-1]
                    if surf and mask_surf:
                        dest[key][(Pref[None,None,:]>interm['Psurf'][:,:,None])] = np.nan

            #create h5 files (pressure grid)
            openh5 = h5py.File(fileh5p,"w")
            print('Writing file '+fileh5p+'...')
            # coordinates
            Pre = openh5.create_dataset("Pressure",data=Pref,\
                                            compression='gzip',compression_opts=comp)
            Lat = openh5.create_dataset("Latitude",data=lat_range*180/np.pi,\
                                        compression='gzip',compression_opts=comp)
            Lon = openh5.create_dataset("Longitude",data=lon_range*180/np.pi,\
                                        compression='gzip',compression_opts=comp)
            # data
            for key in dest.keys():
                tmp_data = openh5.create_dataset(key,data=dest[key],\
                                                 compression='gzip',compression_opts=comp)
            openh5.close()

            #create h5 files (height grid)
            openh5 = h5py.File(fileh5h,"w")
            print('Writing file '+fileh5h+'...')
            Alt = openh5.create_dataset("Altitude",data=grid.Altitude,\
                                            compression='gzip',compression_opts=comp)
            Lat = openh5.create_dataset("Latitude",data=lat_range*180/np.pi,\
                                        compression='gzip',compression_opts=comp)
            Lon = openh5.create_dataset("Longitude",data=lon_range*180/np.pi,\
                                        compression='gzip',compression_opts=comp)
            # data
            for key in dest.keys():
                tmp_data = openh5.create_dataset(key,data=interm[key],\
                                                 compression='gzip',compression_opts=comp)
            openh5.close()

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

def maketable(x,y,z,xname,yname,zname,resultsf,fname):
    #print 2D data from plot to file
    if not os.path.exists(resultsf+'/tables'):
        os.mkdir(resultsf+'/tables')
    f = open(resultsf+'/tables/'+fname,'w')
    f.write('# %s %s %s\n'%(xname,yname,zname))
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            f.write('%#.6e %#.6e %#.6e\n'%(x[i],y[j],z[i,j]))
    f.close()

def vertical_lat(input,grid,output,rg,sigmaref,z,slice=[0,360],save=True,axis=False,csp=500,wind_vectors=False,use_p=True):
    # generic pressure/latitude plot function

    # Set the reference pressure
    if use_p:
        # sigmaref = normalize pressure units
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
        # loni, lati = np.meshgrid(lon[mask_ind],lat)
        # d_lon = np.shape(loni)

        ##################
        #    Averages    #
        ##################
        # Averaging in time and longitude
        if tsp > 1:
            Zonall = np.nanmean(z['value'][:,mask_ind[:,0],:,:],axis=1)
            # Vl = np.nanmean(rg.V[:,:,:,:],axis=1)
            # Wl = np.nanmean(rg.W[:,:,:,:],axis=1)
            Zonallt = np.nanmean(Zonall[:,:,:],axis=2)
            # Vlt = np.nanmean(Vl[:,:,:],axis=2)
            # Wlt = np.nanmean(Wl[:,:,:],axis=2)
            del Zonall
            if wind_vectors == True:
                Vl = np.nanmean(rg.V[:,mask_ind[:,0],:,:],axis=1)
                Wl = np.nanmean(rg.W[:,mask_ind[:,0],:,:],axis=1)
                Vlt = np.nanmean(Vl[:,:,:],axis=2)
                Wlt = np.nanmean(Wl[:,:,:],axis=2)
                del Vl, Wl
        else:
            Zonallt = np.nanmean(z['value'][:,:,:,0][:,mask_ind[:,0],:],axis=1)
            # Vlt = np.nanmean(rg.V[:,:,:,0],axis=1)
            # Wlt = np.nanmean(rg.W[:,:,:,0],axis=1)
            if wind_vectors == True:
                Vlt = np.nanmean(rg.V[:,:,:,0][:,mask_ind[:,0],:],axis=1)
                Wlt = np.nanmean(rg.W[:,:,:,0][:,mask_ind[:,0],:],axis=1)

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
            Zonallt = np.nanmean(Zonall[:,0,:,:],axis=2)
            del Zonall
            if wind_vectors == True:
                Vlt = np.nanmean(Vl[:,0,:,:],axis=2)
                Wlt = np.nanmean(Wl[:,0,:,:],axis=2)
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

    if use_p:
        # need to set desired pressure range (major PITA!)
        # prange = np.where(np.logical_and(rg.Pressure[:,0]>=np.min(Pref),rg.Pressure[:,0]<=np.max(Pref)))
        prange = np.where(rg.Pressure[:,0]>=np.min(Pref))
        ycoord = rg.Pressure[prange[0],0]/1e5
        zvals = Zonallt[:,prange[0]].T
    else:
        hrange = np.where(np.logical_and(rg.Altitude[:,0]>=np.min(sigmaref),rg.Altitude[:,0]<=np.max(sigmaref)))
        ycoord = rg.Altitude[hrange[0],0]/1000
        zvals = Zonallt[:,hrange[0]].T

    # Contour plot
    clevels = 40 # may want to make this adjustable
    # clevels = np.linspace(280,500,45)
    # print(np.max(zvals))
    if isinstance(axis,axes.SubplotBase):
        C = axis.contourf(latp*180/np.pi,ycoord,zvals,clevels,cmap=z['cmap'])
        ax = axis
    elif axis == False:
        plt.figure(figsize=(5,4))
        C = plt.contourf(latp*180/np.pi,ycoord,zvals,clevels,cmap=z['cmap'])
        ax = plt.gca()
    else:
        raise IOError("'axis = {}' but {} is not an axes.SubplotBase instance".format(axis,axis))

    if wind_vectors == True:
        vspacing = np.int(np.shape(rg.lat)[0]/10)
        if use_p:
            wspacing = np.int(np.shape(rg.Pressure)[0]/10)
            Vlt = Vlt[:,prange[0]]
            Wlt = Wlt[:,prange[0]]
            yqcoord = rg.Pressure[::wspacing,0][prange[0]]
        else:
            wspacing = np.int(np.shape(rg.Altitude)[0]/10)
            Vlt = Vlt[:,hrange[0]]
            Wlt = Wlt[:,hrange[0]]
            yqcoord = rg.Altitude[::wspacing,0][hrange[0]]

        Vq = Vlt[::vspacing,::wspacing].ravel()
        Wq = Wlt[::vspacing,::wspacing].ravel()
        #preq = rg.Pressure[:,0][::spacing,::spacing].ravel()
        #latq = lati[::spacing,::spacing].ravel()
        latq, preq = np.meshgrid(rg.lat[::vspacing,0],yqcoord)
        del Vlt, Wlt
        plt.quiver(latq.ravel(),preq.ravel()/1e5,Vq,Wq,color='0.5')

    clb = plt.colorbar(C,extend='both',ax=ax)
    clb.set_label(z['label'])
    if isinstance(csp,list) or isinstance(csp,tuple):
        levp = csp
    else:
        if csp == 'match':
            levp = 40
        else:
            if use_p:
                levp = np.arange(np.ceil(np.nanmin(Zonallt[:,prange[0]])/csp)*csp,np.floor(np.nanmax(Zonallt[:,prange[0]])/csp)*csp,csp)
            else:
                levp = np.arange(np.ceil(np.nanmin(Zonallt[:,hrange[0]])/csp)*csp,np.floor(np.nanmax(Zonallt[:,hrange[0]])/csp)*csp,csp)

    c2 = ax.contour(latp*180/np.pi,ycoord,zvals,levels=levp,colors='w',linewidths=1)
    plt.clabel(c2,inline=False,fontsize=6,fmt='%d',use_clabeltext=True)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    #ax.invert_yaxis()
    # plt.quiver(latq.ravel(),preq.ravel()/1e5,Vq/np.max(Vq),Wq/np.max(Wq),color='0.5')
    if z['plog'] == True:
        ax.set_yscale("log")
    ax.set_xlabel('Latitude (deg)')
    if use_p:
        ax.set_ylabel('Pressure (bar)')
        # ax.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
        ax.set_ylim(np.max(rg.Pressure[prange[0],0])/1e5,np.min(Pref)/1e5)

        if ax.get_ylim()[1] > ax.get_ylim()[0]:
            ax.invert_yaxis()
    else:
        ax.set_ylabel('Altitude (km)')

    if len(slice) == 2:
        ax.set_title('Time = %#.3f-%#.3f days, Lon = (%#.3f,%#.3f)'%(output.time[0],output.time[-1],slice[0],slice[1]),fontsize=10)
    else:
        ax.set_title('Time = %#.3f-%#.3f days, Lon = (%#.3f,)'%(output.time[0],output.time[-1],slice[0]),fontsize=10)

    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    if use_p:
        z['name'] += '_p'
    else:
        z['name'] += '_h'
    if save == True:
        # save the plot to file designated by z
        if len(slice)==2:
            plt.savefig(input.resultsf+'/figures/%s_ver_i%d_l%d_lon%#.2f-%#.2f.pdf'%(z['name'],output.ntsi,output.nts,slice[0],slice[1]))
        else:
            plt.savefig(input.resultsf+'/figures/%s_ver_i%d_l%d_lon%#.2f.pdf'%(z['name'],output.ntsi,output.nts,slice[0]))
        plt.close()
    if z['mt'] == True:
        if len(slice)==2:
            fname = '%s_ver_i%d_l%d_lon%#.2f-%#.2f.dat'%(z['name'],output.ntsi,output.nts,slice[0],slice[1])
        else:
            fname = '%s_ver_i%d_l%d_lon%#.2f.dat'%(z['name'],output.ntsi,output.nts,slice[0])
        if use_p:
            maketable(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,Zonallt[:,prange[0]],'Latitude(d)','Pressure(bar)',z['name'],input.resultsf,fname)
        else:
            maketable(latp*180/np.pi,rg.Altitude[hrange[0],0],Zonallt[:,hrange[0]],'Latitude(d)','Altitude(m)',z['name'],input.resultsf,fname)

def vertical_lon(input,grid,output,rg,sigmaref,z,slice=[0,360],save=True,axis=False,csp=500,wind_vectors=False,use_p=True):
    # generic pressure/longitude plot function

    # Set the reference pressure
    if use_p:
        # sigmaref = normalize pressure units
        Pref = input.P_Ref*sigmaref

    d_sig = np.size(sigmaref)
    tsp = output.nts-output.ntsi+1

    if not isinstance(slice,list):
        raise IOError("'slice' argument must be a list")

    lat = z['lat']
    lon = z['lon']
    if len(slice) == 2:
        # Set the latitude-longitude grid
        if slice[1]-slice[0] > 180:
            raise IOError("'slice' cannot exceed a range of 180 degrees")
        if slice[1] < slice[0]:
            raise IOError("'slice' values must be arranged (small, large) because I am too lazy to code it better")
        # if slice[0] < 0:
        #     mask_ind = np.logical_or((lon-360)>=slice[0],lon<=slice[1])
        # else:
        #     mask_ind = np.logical_and(lon>=slice[0],lon<=slice[1])
        mask_ind = np.logical_and(lat>=slice[0],lat<=slice[1])
        # loni, lati = np.meshgrid(lon,lat[mask_ind])
        # d_lat = np.shape(lati)

        ##################
        #    Averages    #
        ##################

        # Averaging in time and latitude (weighted by a cosine(lat))
        if tsp > 1:
            Meridl = np.nanmean(z['value'][mask_ind[:,0],:,:,:]*np.cos(lat[mask_ind[:,0],0]*np.pi/180)[:,None,None,None],axis=0)
            # Vl = np.nanmean(rg.V[:,:,:,:],axis=1)
            # Wl = np.nanmean(rg.W[:,:,:,:],axis=1)
            Meridlt = np.nanmean(Meridl[:,:,:],axis=2)
            # Vlt = np.nanmean(Vl[:,:,:],axis=2)
            # Wlt = np.nanmean(Wl[:,:,:],axis=2)
            del Meridl
            if wind_vectors == True:
                Ul = np.nanmean(rg.U[mask_ind[:,0],:,:,:]*np.cos(lat[mask_ind[:,0],0]*np.pi/180)[:,None,None,None],axis=0)
                Wl = np.nanmean(rg.W[mask_ind[:,0],:,:,:]*np.cos(lat[mask_ind[:,0],0]*np.pi/180)[:,None,None,None],axis=0)
                Ult = np.nanmean(Ul[:,:,:],axis=2)
                Wlt = np.nanmean(Wl[:,:,:],axis=2)
                del Ul, Wl
        else:
            Meridlt = np.nanmean(z['value'][:,:,:,0][mask_ind[:,0],:,:]*np.cos(lat[mask_ind[:,0],0]*np.pi/180)[:,None,None],axis=0)
            # Ult = np.nanmean(rg.V[:,:,:,0],axis=1)
            # Wlt = np.nanmean(rg.W[:,:,:,0],axis=1)
            if wind_vectors == True:
                Ult = np.nanmean(rg.U[:,:,:,0][mask_ind[:,0],:,:]*np.cos(lat[mask_ind[:,0],0]*np.pi/180)[:,None,None],axis=0)
                Wlt = np.nanmean(rg.W[:,:,:,0][mask_ind[:,0],:,:]*np.cos(lat[mask_ind[:,0],0]*np.pi/180)[:,None,None],axis=0)

    elif len(slice) == 1:
        if slice[0] in lat:
            Meridl = z['value'][lat[:,0]==slice[0],:,:,:]
            if wind_vectors == True:
                Ul = rg.U[lat[:,0]==slice[0],:,:,:]
                Wl = rg.W[lat[:,0]==slice[0],:,:,:]
        else:
            Meridl = np.zeros((1,len(lon),d_sig,tsp))
            if wind_vectors == True:
                Ul = np.zeros((1,len(lon),d_sig,tsp))
                Wl = np.zeros((1,len(lon),d_sig,tsp))
            # interpolate to slice given
            for t in tsp:
                for lev in np.arange(d_sig):
                    Meridl[0,:,lev,tsp] = interp.griddata(np.vstack([lon,lat]).T,z['value'][:,:,lev,tsp],(lon,slice[0]))
                    if wind_vectors == True:
                        Ul[0,:,lev,tsp] = interp.griddata(np.vstack([lon,lat]).T,rg.U[:,:,lev,tsp],(lon,slice[0]))
                        Wl[0,:,lev,tsp] = interp.griddata(np.vstack([lon,lat]).T,rg.W[:,:,lev,tsp],(lon,slice[0]))


        # Averaging in time
        if tsp > 1:
            Meridlt = np.nanmean(Meridl[0,:,:,:],axis=2)
            del Meridl
            if wind_vectors == True:
                Ult = np.nanmean(Ul[0,:,:,:],axis=2)
                Wlt = np.nanmean(Wl[0,:,:,:],axis=2)
                del Ul, Wl
        else:
            Meridlt = Meridl[0,:,:,0]
            if wind_vectors == True:
                Ult = Ul[0,:,:,0]
                Wlt = Wl[0,:,:,0]

    else:
        raise IOError("'slice' must have 1 or 2 values")

    #################
    # Create figure #
    #################
    # Longitude
    lonp = lon[:,0]*np.pi/180

    if use_p:
        # need to set desired pressure range (major PITA!)
        # prange = np.where(np.logical_and(rg.Pressure[:,0]>=np.min(Pref),rg.Pressure[:,0]<=np.max(Pref)))
        prange = np.where(rg.Pressure[:,0]>=np.min(Pref))
        ycoord = rg.Pressure[prange[0],0]/1e5
        zvals = Meridlt[:,prange[0]].T
    else:
        hrange = np.where(np.logical_and(rg.Altitude[:,0]>=np.min(sigmaref),rg.Altitude[:,0]<=np.max(sigmaref)))
        ycoord = rg.Altitude[hrange[0],0]/1000
        zvals = Meridlt[:,hrange[0]].T

    # Contour plot
    clevels = 40 # may want to make this adjustable
    # clevels = np.linspace(-100,6400,66)
    # print(np.min(zvals),np.max(zvals))
    if isinstance(axis,axes.SubplotBase):
        C = axis.contourf(lonp*180/np.pi,ycoord,zvals,clevels,cmap=z['cmap'])
        ax = axis
    elif axis == False:
        plt.figure(figsize=(5,4))
        C = plt.contourf(lonp*180/np.pi,ycoord,zvals,clevels,cmap=z['cmap'])
        ax = plt.gca()
    else:
        raise IOError("'axis = {}' but {} is not an axes.SubplotBase instance".format(axis,axis))

    if wind_vectors == True:
        vspacing = np.int(np.shape(rg.lon)[0]/10)
        if use_p:
            wspacing = np.int(np.shape(rg.Pressure)[0]/10)
            Ult = Ult[:,prange[0]]
            Wlt = Wlt[:,prange[0]]
            yqcoord = rg.Pressure[::wspacing,0][prange[0]]
        else:
            wspacing = np.int(np.shape(rg.Altitude)[0]/10)
            Ult = Ult[:,hrange[0]]
            Wlt = Wlt[:,hrange[0]]
            yqcoord = rg.Altitude[::wspacing,0][hrange[0]]

        Uq = Ult[::vspacing,::wspacing].ravel()
        Wq = Wlt[::vspacing,::wspacing].ravel()
        #preq = rg.Pressure[:,0][::spacing,::spacing].ravel()
        #latq = lati[::spacing,::spacing].ravel()
        lonq, preq = np.meshgrid(rg.lon[::vspacing,0],yqcoord)
        del Ult, Wlt
        plt.quiver(lonq.ravel(),preq.ravel()/1e5,Uq,Wq,color='0.5')

    clb = plt.colorbar(C,extend='both',ax=ax)
    clb.set_label(z['label'])
    if isinstance(csp,list) or isinstance(csp,tuple):
        levp = csp
    else:
        if csp == 'match':
            levp = 40
        else:
            if use_p:
                levp = np.arange(np.ceil(np.nanmin(Meridlt[:,prange[0]])/csp)*csp,np.floor(np.nanmax(Meridlt[:,prange[0]])/csp)*csp,csp)
            else:
                levp = np.arange(np.ceil(np.nanmin(Meridlt[:,hrange[0]])/csp)*csp,np.floor(np.nanmax(Meridlt[:,hrange[0]])/csp)*csp,csp)

    c2 = ax.contour(lonp*180/np.pi,ycoord,zvals,levels=levp,colors='w',linewidths=1)
    plt.clabel(c2,inline=False,fontsize=6,fmt='%d',use_clabeltext=True)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    #ax.invert_yaxis()
    # plt.quiver(latq.ravel(),preq.ravel()/1e5,Vq/np.max(Vq),Wq/np.max(Wq),color='0.5')
    if z['plog'] == True:
        ax.set_yscale("log")
    ax.set_xlabel('Longitude (deg)')
    if use_p:
        ax.set_ylabel('Pressure (bar)')
        ax.plot(lonp*180/np.pi,np.zeros_like(lonp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
        ax.set_ylim(np.max(rg.Pressure[prange[0],0])/1e5,np.min(Pref)/1e5)

        if ax.get_ylim()[1] > ax.get_ylim()[0]:
            ax.invert_yaxis()
    else:
        ax.set_ylabel('Altitude (km)')

    if len(slice) == 2:
        ax.set_title('Time = %#.3f-%#.3f days, Lat = (%#.3f,%#.3f)'%(output.time[0],output.time[-1],slice[0],slice[1]),fontsize=10)
    else:
        ax.set_title('Time = %#.3f-%#.3f days, Lat = (%#.3f,)'%(output.time[0],output.time[-1],slice[0]),fontsize=10)

    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    if use_p:
        z['name'] += '_p'
    else:
        z['name'] += '_h'
    if save == True:
        # save the plot to file designated by z
        if len(slice)==2:
            plt.savefig(input.resultsf+'/figures/%s_ver_i%d_l%d_lat%#.2f-%#.2f.pdf'%(z['name'],output.ntsi,output.nts,slice[0],slice[1]))
        else:
            plt.savefig(input.resultsf+'/figures/%s_ver_i%d_l%d_lat%#.2f.pdf'%(z['name'],output.ntsi,output.nts,slice[0]))
        plt.close()
    if z['mt'] == True:
        if len(slice)==2:
            fname = '%s_ver_i%d_l%d_lat%#.2f-%#.2f.dat'%(z['name'],output.ntsi,output.nts,slice[0],slice[1])
        else:
            fname = '%s_ver_i%d_l%d_lat%#.2f.dat'%(z['name'],output.ntsi,output.nts,slice[0])
        if use_p:
            maketable(lonp*180/np.pi,rg.Pressure[prange[0],0]/1e5,Meridlt[:,prange[0]],'Longitude(d)','Pressure(bar)',z['name'],input.resultsf,fname)
        else:
            maketable(lonp*180/np.pi,rg.Altitude[hrange[0],0],Meridlt[:,hrange[0]],'Longitude(d)','Altitude(m)',z['name'],input.resultsf,fname)

def horizontal_lev(input,grid,output,rg,Plev,z,save=True,axis=False,wind_vectors=False,use_p=True):
    # Set the latitude-longitude grid.
    loni, lati = np.meshgrid(rg.lon[:,0],rg.lat[:,0])

    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    zlev = np.zeros(np.shape(loni)+(tsp,))
    if wind_vectors == True:
         Uii = np.zeros(np.shape(loni)+(tsp,))
         Vii = np.zeros(np.shape(loni)+(tsp,))

    if use_p:
        vcoord = rg.Pressure
    else:
        vcoord = rg.Altitude

    for t in np.arange(tsp):
        # interpolate
        # above is index of higher pressure bound (lower index)
        # below is index of lower pressure bound (higher index)
        if np.size(vcoord[vcoord>Plev]) == 0:
            above = np.where(vcoord==np.max(vcoord))[0][0]
            if use_p:
                below = above + 1
            else:
                below = above - 1
        elif np.size(vcoord[vcoord<Plev]) == 0:
            below = np.where(vcoord==np.min(vcoord))[0][0]
            if use_p:
                above = below - 1
            else:
                above = below + 1
        else:
            above = np.where(vcoord==np.min(vcoord[vcoord>Plev]))[0][0]
            below = np.where(vcoord==np.max(vcoord[vcoord<Plev]))[0][0]

        if len(np.shape(z['value'])) == 4:
            # single level of a 3D field
            zlev[:,:,t] = (z['value'][:,:,below,t]*(vcoord[above,t] - Plev)\
                + z['value'][:,:,above,t]*(Plev - vcoord[below,t]))\
                / (vcoord[above,t]-vcoord[below,t])
            if use_p:
                title = 'Time = %#.3f-%#.3f days, Plev = %#.3f bar'%(output.time[0],output.time[-1],Plev/1e5)
                fname = '%s_lev%#.3fmbar_i%d_l%d'%(z['name'],Plev/100,output.ntsi,output.nts)
            else:
                title = 'Time = %#.3f-%#.3f days, lev = %#.3f km'%(output.time[0],output.time[-1],Plev/1e3)
                fname = '%s_lev%#.3fkm_i%d_l%d'%(z['name'],Plev/1000,output.ntsi,output.nts)
        elif len(np.shape(z['value'])) == 3:
            # 2D field, only one level available (e.g., insolation or surface temp)
            zlev[:,:,t] = z['value'][:,:,t]
            title = 'Time = %#.3f-%#.3f days'%(output.time[0],output.time[-1])
            fname = '%s_i%d_l%d'%(z['name'],output.ntsi,output.nts)
            wind_vectors = False

        if wind_vectors == True:
            Uii[:,:,t] = (rg.U[:,:,below,t]*(vcoord[above,t] - Plev)\
                + rg.U[:,:,above,t]*(Plev - vcoord[below,t]))\
                / (vcoord[above,t]-vcoord[below,t])
            Vii[:,:,t] = (rg.V[:,:,below,t]*(vcoord[above,t] - Plev)\
                + rg.V[:,:,above,t]*(Plev - vcoord[below,t]))\
                / (vcoord[above,t]-vcoord[below,t])

    # Averaging in time
    if tsp > 1:
        zlevt = np.nanmean(zlev,axis=2)
        del zlev
        if wind_vectors == True:
            Uiii = np.nanmean(Uii,axis=2)
            Viii = np.nanmean(Vii,axis=2)
            del Uii, Vii
    else:
        zlevt = zlev[:,:,0]
        del zlev
        if wind_vectors == True:
            Uiii = Uii[:,:,0]
            Viii = Vii[:,:,0]
            del Uii, Vii

    # smoothing
    zlevt = ndimage.gaussian_filter(zlevt,sigma=2,order=0)

    #################
    # Create Figure #
    #################

    lonp = rg.lon[:,0]
    latp = rg.lat[:,0]

    clevels = 50
    # clevels = np.linspace(900,1470,58)
    if isinstance(axis,axes.SubplotBase):
        if z['llswap']:
            C = axis.contourf(latp,lonp,zlevt.T,clevels,cmap=z['cmap'])
        else:
            C = axis.contourf(lonp,latp,zlevt,clevels,cmap=z['cmap'])
        ax = axis
    elif axis == False:
        plt.figure(figsize=(5,4))
        if z['llswap']:
            C = plt.contourf(latp,lonp,zlevt.T,clevels,cmap=z['cmap'])
        else:
            C = plt.contourf(lonp,latp,zlevt,clevels,cmap=z['cmap'])
        ax = plt.gca()
    else:
        raise IOError("'axis = {}' but {} is not an axes.SubplotBase instance".format(axis,axis))

    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0

    if z['llswap']:
        ax.set_xlabel('Latitude (deg)')
        ax.set_ylabel('Longitude (deg)')
    else:
        ax.set_ylabel('Latitude (deg)')
        ax.set_xlabel('Longitude (deg)')

    ax.set_title(title)

    if wind_vectors == True:
        d_z = np.shape(Uiii)
        spacing = np.int(np.shape(Uiii)[0]/10)
        U = Uiii[::spacing,::spacing].ravel()
        V = Viii[::spacing,::spacing].ravel()
        lonq = loni[::spacing,::spacing].ravel()
        latq = lati[::spacing,::spacing].ravel()
        del Uiii, Viii
        if z['llswap']:
            q = plt.quiver(latq,lonq,V,U,color='0.5')
        else:
            q = plt.quiver(lonq,latq,U,V,color='0.5')
        ax.quiverkey(q,X=0.85,Y=-0.1,U=np.max(np.sqrt(U**2+V**2)),label='%#.2f m/s'%np.max(np.sqrt(U**2+V**2)),labelpos='E')

    clb = plt.colorbar(C)
    clb.set_label(z['label'])

    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    if save == True:
        plt.savefig(input.resultsf+'/figures/'+fname+'.pdf')
        plt.close()
    if z['mt'] == True:
        dname = fname+'.dat'
        if z['name'] == 'temperature-uv':
            zname = 'temperature'
        else:
            zname = z['name']
        maketable(lonp,latp,zlevt.T,'Longitude(d)','Latitude(d)',zname,input.resultsf,dname)

def CurlF(fr,flat,flon,lat_range,lon_range,Altitude,A):
    curlFz = np.zeros_like(fr)
    curlFlat = np.zeros_like(fr)
    curlFlon = np.zeros_like(fr)

    curlFz = -1.0*(np.gradient(flon*np.cos(lat_range[:,None,None]),lat_range,axis=0)\
                  -np.gradient(flat,lon_range,axis=1))
    # curlFz = -1.0*(dFunc_dLat(flon*np.cos(lati),res_deg)-dFunc_dLon(flat,res_deg))
    curlFz /= (np.cos(lat_range[:,None,None])*(A+Altitude[None,None,:]))

    curlFlat = -1.0*(np.gradient(fr,lon_range,axis=1)\
                    -np.gradient((A+Altitude[None,None,:])*flon,Altitude,axis=2))
    #curlFlat = -1.0*(dFunc_dLon(fr,res_deg)/np.cos(lati)-dFunc_dZ((A+alti)*flon,alti,nv,dz))
    curlFlat /= (A+Altitude[None,None,:])

    curlFlon = -1.0*(np.gradient((A+Altitude[None,None,:])*flat,Altitude,axis=2)\
                    -np.gradient(fr,lat_range,axis=0))
    #curlFlon = -1.0*(dFunc_dZ((A+alti)*flat,alti,nv,dz)-dFunc_dLat(fr,res_deg))
    curlFlon /= (A+Altitude[None,None,:])

    return curlFz, curlFlat, curlFlon

def streamf_moc_plot(input,grid,output,rg,sigmaref,save=True,axis=False,wind_vectors=False,mt=False,plog=True):
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
    # prange = np.where(np.logical_and(rg.Pressure>=np.min(Pref),rg.Pressure<=np.max(Pref)))
    prange = np.where(rg.Pressure[:,0]>=np.min(Pref))


    # Contour plot
    if isinstance(axis,axes.SubplotBase):
        C = axis.contourf(rg.lat[:,0],rg.Pressure[prange[0],0]/1e5,sf[:,prange[0]].T,40,cmap = 'viridis')
        ax = axis
    elif axis == False:
        plt.figure(figsize=(5,4))
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
    c2 = ax.contour(rg.lat[:,0],rg.Pressure[:,0]/1e5,sf.T,levels=[0.0],colors='w',linewidths=1)
    clb = plt.colorbar(C)
    clb.set_label(r'Eulerian streamfunction (kg s$^{-1}$)')
    if plog == True:
        ax.set_yscale("log")
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Pressure (bar)')
    # ax.plot(rg.lat[:,0],np.zeros_like(rg.lat[:,0])+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
    # if np.min(rg.Pressure[prange[0],0]) < np.max(output.Pressure[:,grid.nv-1,:]):
    ax.set_ylim(np.max(rg.Pressure[prange[0],0])/1e5,np.min(Pref)/1e5)
    # else:
    #     ax.set_ylim(np.max(rg.Pressure[prange[0],0])/1e5,np.max(output.Pressure[:,grid.nv-1,:])/1e5)
    ax.set_title('Time = %#.1f-%#.1f days, Lon = (0,360)'%(output.time[0],output.time[-1]),fontsize=10)
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    if save == True:
        plt.savefig(input.resultsf+'/figures/streamf_ver_i%d_l%d.pdf'%(output.ntsi,output.nts))
        plt.close()
    if mt == True:
        fname = 'streamf_ver_i%d_l%d.dat'%(output.ntsi,output.nts)
        maketable(latp*180/np.pi,rg.Pressure[prange[0],0]/1e5,Zonallt[:,prange[0]],'Latitude(d)','Pressure(bar)','streamfunc',input.resultsf,fname)

def profile(input,grid,output,z,stride=50):
    # Pref = input.P_Ref*sigmaref
    # d_sig = np.size(sigmaref)

    tsp = output.nts-output.ntsi+1

    for column in np.arange(0,grid.point_num,stride):
        if tsp > 1:
            P = np.mean(output.Pressure[column,:,:],axis=1)
            x = np.mean(z['value'][column,:,:],axis=1)
        else:
            P = output.Pressure[column,:,0]
            x = z['value'][column,:,0]

        plt.semilogy(x,P/1e5,'k-',alpha=0.5,lw=1)
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
    Atot = input.A**2
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
    plt.plot(output.time, -1*S)
    plt.xlabel('Time (days)')
    plt.ylabel('Super-rotation index')
    plt.tight_layout()

    if not os.path.exists(input.resultsf+'/figures'):
            os.mkdir(input.resultsf+'/figures')
    plt.savefig(input.resultsf+'/figures/SRindex_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def Get_Prange(input,grid,rg,args,xtype='lat',use_p=True):
    # Sigma (normalized pressure) values for the plotting
    if not isinstance(args.slice,list):
        raise IOError("'slice' argument must be a list")

    if use_p:
        if (args.vertical_top[0]=='default'):
            args.vertical_top[0] = np.min(rg.Pressure)/100

        if np.max(input.P_Ref)/np.float(args.vertical_top[0]) > 1000:
            sigmaref = np.logspace(np.log10(np.max(rg.Pressure)),np.log10(np.float(args.vertical_top[0])*100),20)/input.P_Ref
        else:
            sigmaref = np.linspace(np.max(rg.Pressure),np.float(args.vertical_top[0])*100,20)/input.P_Ref

    else:
        if (args.vertical_top[0]=='default'):
            sigmaref = grid.Altitude
        else:
            sigmaref = grid.Altitude[np.where(grid.Altitude<=np.float(args.vertical_top[0])*input.Top_altitude)[0]]
    return sigmaref


def RTbalance(input,grid,output):
    #not finished!
    asr = fsw_dn[:,input.nv-1,:]*grid.areasT[:,None]
    olr = flw_up[:,input.nv-1,:]*grid.areasT[:,None]
