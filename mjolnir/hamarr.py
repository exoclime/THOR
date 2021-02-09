#---- code by Russell Deitrick and Urs Schroffenegger -------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.axes as axes
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import hsv_to_rgb
import matplotlib.patheffects as pe

import math

import scipy.interpolate as interp
import scipy.ndimage as ndimage
import os
import h5py
import time
import subprocess as spr
import pyshtools as chairs
import pdb

import pathlib
import psutil

try:
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import pycuda.gpuarray as gpuarray
    has_pycuda = True
except ImportError as e:
    print(e)
    has_pycuda = False

plt.rcParams['image.cmap'] = 'magma'
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'Helvetica-Normal'



class input_new:
    def __init__(self, resultsf, simID):
        fileh5 = resultsf + '/esp_output_planet_' + simID + '.h5'
        if os.path.exists(fileh5):
            openh5 = h5py.File(fileh5,'r')
        else:
            fileh5_old = resultsf + '/esp_output_' + simID + '.h5'
            if os.path.exists(fileh5_old):
                openh5 = h5py.File(fileh5_old,'r')
            else:
                raise IOError(fileh5 + ' or ' + fileh5_old + ' not found!')

        self.resultsf = resultsf
        self.simID = simID

        for key in openh5.keys():
            setattr(self,key,openh5[key][...])

        #special cases (things we test on a lot, etc)
        self.RT = "radiative_transfer" in openh5 and openh5["radiative_transfer"][0] == 1.0
        self.TSRT = "two_streams_radiative_transfer" in openh5 and openh5["two_streams_radiative_transfer"][0] == 1.0
        if self.TSRT:
            if 'alf_w0_g0_per_band' in openh5 and openh5['alf_w0_g0_per_band'][0] == 1.0:
                self.has_w0_g0 = True
            else:
                self.has_w0_g0 = False
        else:
            self.has_w0_g0 = False

        self.chemistry = "chemistry" in openh5 and openh5["chemistry" ][0] == 1

        if not hasattr(self,'surface'):
            self.surface = False
        #some bw compatibility things
        if hasattr(self,'vulcan'):
            self.chemistry = self.vulcan
        if hasattr(self,'diff_fac'):
            self.diff_ang = self.diff_fac
        openh5.close()

def LatLong2Cart(lat, lon):
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return x, y, z


def Cart2LatLong(x, y, z):
    lon = np.arctan2(y, x)
    lat = np.arctan2(z, np.sqrt(x**2 + y**2))
    return lat, lon


def Rot_y(theta, x, y, z):
    xnew = np.cos(theta) * x + np.sin(theta) * z
    znew = -np.sin(theta) * x + np.cos(theta) * z
    return xnew, y, znew


def Rot_z(phi, x, y, z):
    xnew = np.cos(phi) * x - np.sin(phi) * y
    ynew = np.sin(phi) * x + np.cos(phi) * y
    return xnew, ynew, z

class grid_new:
    def __init__(self, resultsf, simID, rotation=False, theta_y=0, theta_z=0):
        fileh5 = resultsf + '/esp_output_grid_' + simID + '.h5'
        if os.path.exists(fileh5):
            openh5 = h5py.File(fileh5,'r')
        else:
            raise IOError(fileh5 + ' not found!')

        self.point_num = np.int(openh5['point_num'][0])
        self.nv = np.int(openh5['nv'][0])
        self.nvi = self.nv + 1
        for key in openh5.keys():
            if key != 'nv' and key != 'point_num':
                if key == 'pntloc':
                    setattr(self,key,np.reshape(openh5[key][...],(self.point_num,6)))
                elif key == 'grad' or key == 'div' or key == 'curlz':
                    setattr(self,key,np.reshape(openh5[key][...],(self.point_num,7,3)))
                elif key == 'lonlat':
                    self.lon = (openh5[key][::2]) % (2*np.pi)
                    self.lat = openh5[key][1::2]
                    if rotation == True:
                        xtmp, ytmp, ztmp = LatLong2Cart(self.lat, self.lon)
                        xtmp1, ytmp1, ztmp1 = Rot_z(theta_z, xtmp, ytmp, ztmp)
                        xtmp2, ytmp2, ztmp2 = Rot_y(theta_y, xtmp1, ytmp1, ztmp1)
                        self.lat, self.lon = Cart2LatLong(xtmp2, ytmp2, ztmp2)
                        self.lon = self.lon % (2 * np.pi)
                else:
                    setattr(self,key,openh5[key][...])
        openh5.close()

class output_new:
    def __init__(self, resultsf, simID, ntsi, nts, input, grid, stride=1):
        #need to test best way to read/create virtual dataset
        #if I read in source files with "r", I need to close the virtual data file
        # and then reopen it with "r" to make data visible
        #if I read in source files with "r+", there is no need to close and reopen

        #some basic info
        self.ntsi = ntsi
        self.nts = nts
        # 1-D arrays
        self.time = np.zeros(nts - ntsi + 1)
        self.nstep = np.zeros(nts - ntsi + 1)
        self.ConvData = np.zeros(nts - ntsi + 1)
        self.GlobalE = np.zeros(nts - ntsi + 1)
        self.GlobalEnt = np.zeros(nts - ntsi + 1)
        self.GlobalMass = np.zeros(nts - ntsi + 1)
        self.GlobalAMx = np.zeros(nts - ntsi + 1)
        self.GlobalAMy = np.zeros(nts - ntsi + 1)
        self.GlobalAMz = np.zeros(nts - ntsi + 1)

        # dictionary of field variables (2-D or 3-D arrays)
        # key is the key in h5 file, value is desired attribute name in this class
        outputs = {'Rho': 'Rho', 'Pressure': 'Pressure', 'Mh': 'Mh', 'Wh': 'Wh'}

        if input.output_mean:
            outputs['Rho_mean'] = 'Rho_mean'
            outputs['Pressure_mean'] = 'Pressure_mean'
            outputs['Mh_mean'] = 'Mh_mean'
            outputs['Wh_mean'] = 'Wh_mean'

        #add things to outputs that can be checked for in input or grid
        # if input.RT or input.TSRT:
        #     outputs['Qheat'] = 'qheat'

        if input.RT:
            outputs['tau'] = 'tau'  #have to be careful about slicing this (sw = ::2, lw = 1::2)
            outputs['flw_up'] = 'flw_up'
            outputs['flw_dn'] = 'flw_dn'
            outputs['fsw_dn'] = 'fsw_dn'
            # outputs['DGQheat'] = 'DGqheat'

        if input.TSRT:
            outputs['F_up_tot'] = 'f_up_tot'
            outputs['F_down_tot'] = 'f_down_tot'
            outputs['F_dir_tot'] = 'f_dir_tot'
            outputs['F_net'] = 'f_net'
            outputs['Alf_Qheat'] = 'TSqheat'
            outputs['F_up_TOA_spectrum'] = 'spectrum'
            outputs['alf_spectrum'] = 'incoming_spectrum'
            outputs['lambda_wave'] = 'wavelength'

        # calc volume element
        Atot = input.A**2
        solid_ang = grid.areasT/Atot
        if input.DeepModel:
            rint = ((input.A + grid.Altitudeh[1:])**3 - (input.A + grid.Altitudeh[:-1])**3) / 3.0
        else:
            rint = Atot*(grid.Altitudeh[1:] - grid.Altitudeh[:-1])
        self.Vol0 = solid_ang[:, None] * rint[None, :]

        for t in np.arange(ntsi - 1, nts, stride):
            fileh5 = resultsf + '/esp_output_' + simID + '_' + np.str(t + 1) + '.h5'
            if os.path.exists(fileh5):
                openh5 = h5py.File(fileh5, 'r+')
            else:
                raise IOError(fileh5 + ' not found!')

            if t == ntsi - 1:
                # add things to outputs dictionary that require checking for existence in openh5
                if 'Qheat' in openh5.keys():
                    outputs['Qheat'] = 'qheat'
                if 'DGQheat' in openh5.keys():
                    outputs['DGQheat'] = 'DGqheat'
                if 'Etotal' in openh5.keys():
                    self.ConvData[t-ntsi+1] = True
                    outputs['Etotal'] = 'Etotal'
                    outputs['Mass'] = 'Mass'
                    outputs['AngMomx'] = 'AngMomx'
                    outputs['AngMomy'] = 'AngMomy'
                    outputs['AngMomz'] = 'AngMomz'
                    outputs['Entropy'] = 'Entropy'
                else:
                    print('Warning: conservation diagnostics not available in file %s' % fileh5)
                    self.ConvData[t-ntsi+1] = False

                if 'insol' in openh5.keys():
                    outputs['insol'] = 'Insol'
                if 'tracer' in openh5.keys():
                    outputs['tracer'] = 'tracer' #ch4 ::5, co 1::5, h2o 2::5, co2 3::5, nh3 4::5
                if 'Tsurface' in openh5.keys():
                    outputs['Tsurface'] = 'Tsurface'
                if 'Rd' in openh5.keys():
                    self.ConstRdCp = False
                    outputs['Rd'] = 'Rd'
                    outputs['Cp'] = 'Cp'
                else:
                    self.ConstRdCp = True

                # for rename transition of col_mu_star, can be removed later
                if 'cos_zenith_angles' in openh5.keys():
                    outputs['cos_zenith_angles'] = 'mustar'
                elif 'zenith_angles' in openh5.keys():
                    outputs['zenith_angles'] = 'mustar'
                elif 'col_mu_star' in openh5.keys():
                    outputs['col_mu_star'] = 'mustar'

                if input.TSRT:
                    if 'w0_band' in openh5.keys():
                        outputs['w0_band'] = 'w0_band'
                    if 'g0_band' in openh5.keys():
                        outputs['g0_band'] = 'g0_band'

                #create VDS layout shape
                for key in outputs.keys():
                    setattr(self, outputs[key]+'_lo', h5py.VirtualLayout(shape=
                                (np.shape(openh5[key])[0],nts-ntsi+1),dtype='f8'))

            #shove source reference into layout
            for key in outputs.keys():
                getattr(self, outputs[key]+'_lo')[:,t-ntsi+1] = h5py.VirtualSource(openh5[key])

            # 1-D arrays, unlikely to overflow memory
            self.time[t-ntsi+1] = openh5['simulation_time'][0] / 86400
            self.nstep[t-ntsi+1] = openh5['nstep'][0]
            if self.ConvData[t-ntsi+1]:
                self.GlobalE[t - ntsi + 1] = openh5['GlobalE'][0]
                self.GlobalMass[t - ntsi + 1] = openh5['GlobalMass'][0]
                self.GlobalAMx[t - ntsi + 1] = openh5['GlobalAMx'][0]
                self.GlobalAMy[t - ntsi + 1] = openh5['GlobalAMy'][0]
                self.GlobalAMz[t - ntsi + 1] = openh5['GlobalAMz'][0]
                self.GlobalEnt[t - ntsi + 1] = openh5['GlobalEnt'][0]

            openh5.close()

        fileVDS = resultsf+'/esp_output_'+simID+'_'+str(ntsi)+'_'+str(nts)+'_VDS.h5'
        self.openVDS = h5py.File(fileVDS,'w',libver='latest')

        #create data sets.
        #will need to reshape them (grid.point_num,grid.nv,nts-ntsi+1) in mjolnir
        #self .Pressure = self.openVDS.create_virtual_dataset('Pressure',Pressure_lo,fillvalue=0)
        for key in outputs.keys():
            setattr(self, outputs[key], self.openVDS.create_virtual_dataset(outputs[key],getattr(self, outputs[key]+'_lo'),fillvalue=0))

    # for use at end of plotting... closes h5 file with virtual data set
    def closeVDS(self):
        self.openVDS.close()

    #loads arrays from disk into memory and reshapes them (destroys reference to virtual data set)
    def load_reshape(self,grid,keys):
        tlen = self.nts-self.ntsi+1
        for key in keys:
            if key in self.openVDS.keys():
                data_ref = getattr(self, key)
                if isinstance(data_ref,h5py.Dataset):
                    if data_ref.size*8 > psutil.virtual_memory().available:
                        raise IOError("Dataset %s too large for memory; not loading"%key)
                    data = data_ref[...]
                    if np.shape(data) == (grid.point_num*grid.nv,tlen):
                        data = np.reshape(data,(grid.point_num,grid.nv,tlen))
                    elif np.shape(data) == (grid.point_num*grid.nvi,tlen):
                        data = np.reshape(data,(grid.point_num,grid.nvi,tlen))
                    # need to handle special formats (Mh, tracers...)
                    elif np.shape(data) == (grid.point_num*grid.nv*3,tlen):
                        data_new = np.zeros((3,grid.point_num,grid.nv,tlen))
                        data_new[0] = np.reshape(data[::3],(grid.point_num,grid.nv,tlen))
                        data_new[1] = np.reshape(data[1::3],(grid.point_num,grid.nv,tlen))
                        data_new[2] = np.reshape(data[2::3],(grid.point_num,grid.nv,tlen))
                        data = data_new
                    if key == 'tracer':  #special case for passive tracers
                        self.ch4 = np.reshape(data[::5],(grid.point_num,grid.nv,tlen))
                        self.co = np.reshape(data[1::5],(grid.point_num,grid.nv,tlen))
                        self.h2o = np.reshape(data[2::5],(grid.point_num,grid.nv,tlen))
                        self.co2 = np.reshape(data[3::5],(grid.point_num,grid.nv,tlen))
                        self.nh3 = np.reshape(data[4::5],(grid.point_num,grid.nv,tlen))
                    elif key == 'tau':
                        self.tau_sw = np.reshape(data[::2],(grid.point_num,grid.nv,tlen))
                        self.tau_lw = np.reshape(data[1::2],(grid.point_num,grid.nv,tlen))
                    elif key == 'spectrum':
                        self.spectrum = np.reshape(data, (grid.point_num,-1,tlen))
                    elif key == 'incoming_spectrum':
                        self.incoming_spectrum = np.reshape(data, (-1,tlen))
                    elif key == 'w0_band':
                        self.w0_band = np.reshape(data, (grid.point_num,grid.nv,-1,tlen))
                    elif key == 'g0_band':
                        self.g0_band = np.reshape(data, (grid.point_num,grid.nv,-1,tlen))
                    elif key == 'Etotal' or key == 'Entropy' or key == 'AngMomz':
                        data = np.reshape(data,(grid.point_num,grid.nv,tlen))/self.Vol0[:,:,None]
                        setattr(self, key, data)
                    else:
                        setattr(self, key, data)
                else:
                    print("Dataset %s has already been loaded into memory"%key)
            else:
                raise IOError("Invalid data set specified: %s"%key)

    def load_reshape_all(self,grid):
        keys = self.openVDS.keys()
        self.load_reshape(grid,keys)
        self.closeVDS()


class rg_out_new:
    def __init__(self, resultsf, simID, ntsi, nts, input, grid, pressure_vert=True, pgrid_ref='auto'):
        RT = False
        if "core_benchmark" in dir(input):
            # need to switch to 'radiative_tranfer' flag
            if input.core_benchmark[0] == 0 and input.RT:
                RT = True
                surf = 0
                if input.surface:
                    surf = 1
            else:
                surf = 0

        chem = 0
        if hasattr(input, "chemistry"):
            if input.chemistry == 1:
                chem = 1

        #some basic info
        self.ntsi = ntsi
        self.nts = nts
        # Read model results
        for t in np.arange(ntsi - 1, nts):
            if pressure_vert == True:
                if pgrid_ref == 'auto':
                    pgrid_folder = resultsf + '/pgrid_%d_%d_1'%(ntsi,nts)
                else:
                    pgrid_folder = resultsf + '/' + pgrid_ref[:-4]
                fileh5 = pgrid_folder + '/regrid_' + simID + '_' + np.str(t + 1)
            else:
                fileh5 = resultsf + '/regrid_height_' + simID + '_' + np.str(t + 1)
            fileh5 += '.h5'
            if os.path.exists(fileh5):
                openh5 = h5py.File(fileh5,'r+')
            else:
                print(fileh5 + ' not found, regridding now with default settings...')
                regrid(resultsf, simID, ntsi, nts, pgrid_ref=pgrid_ref)
                openh5 = h5py.File(fileh5,'r+')

            for key in openh5.keys():
                if t == ntsi - 1:  #set up VDS layout shape
                    if openh5[key].ndim == 1:
                        setattr(self,key,openh5[key][:])
                    else:
                        vars()[key+'_lo'] = h5py.VirtualLayout(shape=
                                            np.shape(openh5[key])+(nts-ntsi+1,),dtype='f8')
                if key+'_lo' in vars():
                    if openh5[key].ndim == 2:
                        vars()[key+'_lo'][:,:,t-ntsi+1] = h5py.VirtualSource(openh5[key])
                    elif openh5[key].ndim == 3:
                        vars()[key+'_lo'][:,:,:,t-ntsi+1] = h5py.VirtualSource(openh5[key])
                    else:
                        raise IOError('Property %s with unknown dimensions in regrid file'%key)

        fileVDS = resultsf+'/regrid_'+simID+'_'+str(ntsi)+'_'+str(nts)+'_VDS.h5'
        self.openVDS = h5py.File(fileVDS,'w',libver='latest')

        for key in openh5.keys():
            if key+'_lo' in vars():
                setattr(self,key,self.openVDS.create_virtual_dataset(key,vars()[key+'_lo'],fillvalue=0))

    def closeVDS(self):
        self.openVDS.close()

    def load(self,keys):
        tlen = self.nts-self.ntsi+1
        for key in keys:
            if key in self.openVDS.keys():
                data_ref = getattr(self, key)
                if isinstance(data_ref,h5py.Dataset):
                    if data_ref.size*8 > psutil.virtual_memory().available:
                        raise IOError("Dataset %s too large for memory; not loading"%key)
                    data = data_ref[...]
                    setattr(self,key,data)
                else:
                    print("Dataset %s has already been loaded into memory"%key)
            else:
                print("Invalid data set specified: %s; not loading"%key)

    def load_all(self):
        keys = self.openVDS.keys()
        self.load(keys)
        self.closeVDS()


class GetOutput:
    def __init__(self, resultsf, simID, ntsi, nts, stride=1, openrg=0, pressure_vert=True, rotation=False, theta_y=0, theta_z=0, pgrid_ref='auto'):
        self.input = input_new(resultsf, simID)
        self.grid = grid_new(resultsf, simID, rotation=rotation, theta_y=theta_y, theta_z=theta_z)
        if openrg == 1:  #checking for regrid before output prevents file open conflict
            self.rg = rg_out_new(resultsf, simID, ntsi, nts, self.input, self.grid,
                             pressure_vert=pressure_vert, pgrid_ref=pgrid_ref)
        self.output = output_new(resultsf, simID, ntsi, nts, self.input, self.grid, stride=stride)

def define_Pgrid(resultsf, simID, ntsi, nts, stride, overwrite=False):
    pfile = 'pgrid_%d_%d_%d.txt' % (ntsi, nts, stride)
    if not os.path.exists(resultsf + '/' + pfile) or overwrite == True:
        #check for existence of pgrid file

        # first we need the grid size
        fileh5 = resultsf + '/esp_output_grid_' + simID + '.h5'
        if os.path.exists(fileh5):
            openh5 = h5py.File(fileh5,'r')
        else:
            raise IOError(fileh5 + ' not found!')
        nv = np.int(openh5['nv'][0])
        point_num = np.int(openh5['point_num'][0])
        Altitude = openh5['Altitude'][...]
        Altitudeh = openh5['Altitudeh'][...]
        openh5.close()

        # we also need to know if we included a surface
        fileh5 = resultsf + '/esp_output_planet_' + simID + '.h5'
        if os.path.exists(fileh5):
            openh5 = h5py.File(fileh5,'r')
        else:
            raise IOError(fileh5 + ' not found!')
        if openh5['core_benchmark'][0] > 0:
            surf = 0
        else:
            if "radiative_transfer" in openh5:
                surf = openh5['surface'][0]
            else:
                surf = 0

        openh5.close()

        # now we'll loop over all the files to get the pressure_mean
        num_out = np.int((nts - ntsi) / stride) + 1
        pressure_mean = np.zeros((point_num, nv, num_out))
        for i in np.arange(num_out):
            t = ntsi + stride * i
            fileh5 = resultsf + '/esp_output_' + simID + '_' + np.str(t) + '.h5'
            if os.path.exists(fileh5):
                openh5 = h5py.File(fileh5,'r')
            else:
                raise IOError(fileh5 + ' not found!')

            pressure_mean[:, :, i] = np.reshape(openh5['Pressure_mean'][...], (point_num, nv))
            openh5.close()

        pgrid = np.mean(np.mean(pressure_mean, axis=0), axis=1)
        if surf == 1:
            extrap_low = (Altitudeh[0] - Altitude[1]) / (Altitude[0] - Altitude[1])
            Psurf = pressure_mean[:, 1, :] + extrap_low * (pressure_mean[:, 0, :] - pressure_mean[:, 1, :])
            Psurf_mean = np.mean(np.mean(Psurf, axis=0), axis=0)
            pgrid = np.concatenate((np.array([Psurf_mean]), pgrid))

        f = open(resultsf + '/' + pfile, 'w')
        for j in np.arange(len(pgrid)):
            f.write('%d %#.6e\n' % (j, pgrid[j]))
        f.close()
    else:
        print(pfile + ' already exists in this directory! Skipping...')
        #raise IOError(pfile + ' already exists in this directory!')
    return pfile


if has_pycuda:
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


def create_rg_map(resultsf, simID, idx1, idx2, rotation=False, theta_z=0, theta_y=0):
    outall = GetOutput(resultsf, simID, idx1, idx2, rotation=rotation, theta_z=theta_z, theta_y=theta_y)
    input = outall.input
    grid = outall.grid

    # vector positions of ico grid points
    v_ico = np.vstack([np.cos(grid.lon) * np.cos(grid.lat),
                       np.sin(grid.lon) * np.cos(grid.lat),
                       np.sin(grid.lat)]).T

    # set horizontal angular resolution of latlon grid roughly same as ico grid
    ang_res = (4.0 / 2**(input.glevel - 4)) * np.pi / 180

    # 1D lat and lon arrays
    lat_range = np.arange(-np.pi / 2 + ang_res / 2, np.pi / 2, ang_res)
    lon_range = np.arange(0, 2 * np.pi, ang_res)
    loni, lati = np.meshgrid(lon_range, lat_range)
    num_ll = len(lat_range) * len(lon_range)

    # vector positions of latitude-longitude points
    v_ll = np.vstack([(np.cos(loni) * np.cos(lati)).ravel(),
                      (np.sin(loni) * np.cos(lati)).ravel(),
                      (np.sin(lati)).ravel()]).T

    # arrays for nearest neighbors and weights
    near3 = np.zeros((num_ll, 3), dtype=np.int32)
    weight3 = np.zeros((num_ll, 3))
    near = np.zeros(num_ll, dtype=np.int32)

    # cuda functions and sizes
    find_nearest = regrid_tools.get_function("find_nearest")
    calc_weights = regrid_tools.get_function("calc_weights")
    gridGPU = (np.int(np.floor(num_ll / 256)) + 1, 1, 1)
    blockGPU = (256, 1, 1)

    # nearest neighbor search
    find_nearest(cuda.Out(near), cuda.In(v_ll.ravel()), cuda.In(v_ico.ravel()),
                 np.int32(num_ll), np.int32(grid.point_num),
                 block=blockGPU, grid=gridGPU)
    near3[:, 0] = near

    # next 2 neighbors and weights from Moller-Trumbore algorithm
    calc_weights(cuda.Out(weight3.ravel()), cuda.InOut(near3.ravel()),
                 cuda.In(v_ll.ravel()), cuda.In(v_ico.ravel()),
                 cuda.In(grid.pntloc.ravel()), np.int32(num_ll),
                 np.int32(grid.point_num), block=blockGPU, grid=gridGPU)

    # save into numpy zip file
    np.savez(resultsf + '/regrid_map.npz', lon=lon_range, lat=lat_range,
             near3=near3, weight3=weight3)


def vertical_regrid_field(source_array, nv, x, xnew):
    # handles set up and running of vertical interpolation in pycuda
    vert_lin_interp = regrid_tools.get_function("vert_lin_interp")
    x_non_mono_check = np.zeros(np.shape(source_array)[:-1], dtype=np.int32)
    xnew_non_mono_check = np.zeros(1, dtype=np.int32)
    gridgpu = (np.int(np.floor(len(x_non_mono_check.ravel()) / 256)) + 1, 1, 1)
    blockgpu = (256, 1, 1)

    y = np.ascontiguousarray(source_array.ravel())
    ynew = np.zeros(np.shape(source_array)[:-1] + (len(xnew),))
    vert_lin_interp(cuda.In(x), cuda.In(y), cuda.In(xnew),
                    cuda.Out(ynew.ravel()), np.int32(len(xnew)),
                    np.int32(len(x_non_mono_check.ravel())), np.int32(nv),
                    cuda.InOut(x_non_mono_check.ravel()),
                    cuda.InOut(xnew_non_mono_check), grid=gridgpu, block=blockgpu)
    if (x_non_mono_check == 1).any():
        raise ValueError('Pressure interpolation failed! Pressure not monotonically decreasing with height')
    if (xnew_non_mono_check[0] == 1):
        raise ValueError('Pressure interpolation failed! Destination pgrid not monotonically decreasing with height')
    return ynew


def regrid(resultsf, simID, ntsi, nts, pgrid_ref='auto', overwrite=False, comp=4,
           rotation=False, theta_z=0, theta_y=0, mask_surf=True):
    # New regridding process, smoother and faster than previous one
    # For each point on destination grid (lat-lon), find 3 points (near3) on the ico grid
    # that define a triangle containing the destination point. Then use the
    # Moller-Trumbore algorithm to calculate the weights (weight3) that interpolate the 3
    # points to the destination point. PyCuda is used to accelerate the process.
    # Interpolating each field is then a vectorized process using Numpy, taking
    #field_on_ll_grid = weight3 * field_on_ico_grid[near3]
    # Remapping from height to pressure is done on all columns via pycuda

    if not os.path.exists(resultsf + '/regrid_map.npz') or overwrite:
        # if numpy zip file containing weights d.n.e., then create it
        create_rg_map(resultsf, simID, ntsi, ntsi, rotation=rotation, theta_z=theta_z, theta_y=theta_y)

    # open numpy zip file containing weights
    rg_map = np.load(resultsf + '/regrid_map.npz')
    lon_range = rg_map['lon']
    lat_range = rg_map['lat']
    near3 = rg_map['near3']
    weight3 = rg_map['weight3']

    # open first file to get some basic properties
    outall = GetOutput(resultsf, simID, ntsi, ntsi, rotation=rotation, theta_z=theta_z, theta_y=theta_y)
    input = outall.input
    grid = outall.grid
    output = outall.output
    output.load_reshape_all(grid)

    print('Regrid data in folder ' + resultsf + '...\n')

    # handling of vertical coordinate
    if pgrid_ref == 'auto':
        pgrid_file = resultsf+'/'+define_Pgrid(resultsf,simID,ntsi,nts,1,overwrite=overwrite)
        # try:
        #     files_found = spr.check_output('ls ' + resultsf + '/pgrid_*.txt', shell=True).split()
        # except:
        #     raise IOError('No pgrid file found in "%s/", please run "pgrid -i <first> -l <last> %s"' % (resultsf, resultsf))
        # if len(files_found) > 1:
        #     raise IOError('Multiple pgrid files found in "%s/", please specify with -pgrid flag' % resultsf)
        # else:
        #     pgrid_file = files_found[0].decode()
    else:
        if os.path.exists(resultsf + '/' + pgrid_ref):
            pgrid_file = resultsf + '/' + pgrid_ref
        else:
            raise IOError('pgrid file %s not found!' % (resultsf + '/' + pgrid_ref))
    pgrid_folder = pgrid_file[:-4]  #define folder as pgrid without file suffix
    print('Vertical coordinate = pressure from file %s' % pgrid_file)
    Pref = np.loadtxt(pgrid_file, unpack=True, usecols=1)
    p_output_path = pathlib.Path(pgrid_folder)
    if not p_output_path.exists():
        p_output_path.mkdir()

    # destination grid and set some dimensions
    loni, lati = np.meshgrid(lon_range, lat_range)
    d_lon = np.shape(loni)
    tsp = nts - ntsi + 1
    d_sig = np.size(Pref)

    # handling of module properties
    RT = 0  # gray radiative transfer module
    if hasattr(input, "core_benchmark"):
        # need to switch to 'radiative_tranfer' flag
        if input.core_benchmark[0] == 0:
            if input.RT:
                RT = 1
                surf = 0
                if input.surface:
                    surf = 1
                    extrap_low = (grid.Altitudeh[0] - grid.Altitude[1]) / (grid.Altitude[0] - grid.Altitude[1])
                    Psurf = output.Pressure[:, 1, :] + extrap_low * (output.Pressure[:, 0, :] - output.Pressure[:, 1, :])
            else:
                RT = 0
                surf = 0

        else:
            surf = 0

    chem = 0  # chemistry module
    if hasattr(input, "chemistry"):
        if input.chemistry == 1:
            chem = 1

    # begin regrid loop over times
    for t in np.arange(ntsi, nts + 1):
        # check for existing h5 files
        proceed = 0
        skip_height_file = 0
        fileh5p = pgrid_folder + '/regrid_' + simID + '_' + np.str(t)

        fileh5h = resultsf + '/regrid_height_' + simID + '_' + np.str(t)
        if rotation == True:  # tell the user they asked for a rotated grid
            print('Applied rotation (theta_z,theta_y) = (%f,%f) to grid\n'
                  % (theta_z * 180 / np.pi, theta_y * 180 / np.pi))
        fileh5p += '.h5'
        fileh5h += '.h5'
        if os.path.exists(fileh5p) and os.path.exists(fileh5h):  # regrid files exist
            if overwrite == True:  # overwrite existing files
                proceed = 1
            else:  # skip existing files
                print(fileh5p + ' and ' + fileh5h + ' already present! Skipping time = %d' % t)
        else:  # >= one file missing
            proceed = 1
            if os.path.exists(fileh5h):  #height file already exists
                if overwrite == True:
                    skip_height_file = 0
                else:
                    print(fileh5h + ' already present! Skipping height file for time = %d' % t)
                    skip_height_file = 1

        # begin regridding
        if proceed == 1:
            print('Regridding time = %d...' % t)

            # read in one file at a time to prevent mem overflow
            outall = GetOutput(resultsf, simID, t, t, rotation=rotation, theta_z=theta_z, theta_y=theta_y)
            output = outall.output
            output.load_reshape_all(grid)

            # interpolate in z to middle of layers using this
            interpz = (grid.Altitude - grid.Altitudeh[:-1]) / (grid.Altitudeh[1:] - grid.Altitudeh[:-1])

            # set up dictionary with source arrays (add new quantities here!)
            source = {'Temperature': (output.Pressure / (output.Rd * output.Rho))[:, :, 0],
                      'W': (output.Wh[:, :-1, 0] + (output.Wh[:, 1:, 0] - output.Wh[:, :-1, 0]) * interpz[None, :]) / output.Rho[:, :, 0],
                      'Rho': output.Rho[:, :, 0],
                      'Mh': output.Mh[:, :, :, 0],
                      'Pressure': output.Pressure[:, :, 0],
                      'Rd': output.Rd[:, :, 0],
                      'Cp': output.Cp[:, :, 0],
                      'Temperature_mean': (output.Pressure_mean/(output.Rd*output.Rho_mean))[:,:,0],
                      'W_mean': (output.Wh_mean[:, :-1, 0] + (output.Wh_mean[:, 1:, 0]-output.Wh_mean[:, :-1, 0])*interpz[None, :])/output.Rho_mean[:, :, 0],
                      'Rho_mean': output.Rho_mean[:, :, 0],
                      'Mh_mean': output.Mh_mean[:, :, :, 0],
                      'Pressure_mean': output.Pressure_mean[:, :, 0]}

            if input.RT or input.TSRT:
                if hasattr(output,'qheat'):
                    source['qheat'] = output.qheat[:, :, 0]

            if input.RT == 1:
                source['flw_up'] = output.flw_up[:, :-1, 0] + (output.flw_up[:, 1:, 0] - output.flw_up[:, :-1, 0]) * interpz[None, :]
                source['flw_dn'] = output.flw_dn[:, :-1, 0] + (output.flw_dn[:, 1:, 0] - output.flw_dn[:, :-1, 0]) * interpz[None, :]
                source['fsw_dn'] = output.fsw_dn[:, :-1, 0] + (output.fsw_dn[:, 1:, 0] - output.fsw_dn[:, :-1, 0]) * interpz[None, :]
                fnet_tmp = output.flw_up[:,:,0] - (output.flw_dn[:,:,0]+output.fsw_dn[:,:,0])
                source['DGf_net'] = fnet_tmp[:,:-1] + (fnet_tmp[:, 1:] - fnet_tmp[:, :-1]) * interpz[None, :]
                source['tau_sw'] = output.tau_sw[:, :, 0]
                source['tau_lw'] = output.tau_lw[:, :, 0]
                if hasattr(output,'DGqheat'):
                    source['DGqheat'] = output.DGqheat[:, :, 0]
                source['insol'] = output.Insol[:, 0]
                if surf == 1:
                    source['Tsurface'] = output.Tsurface[:, 0]
                    source['Psurf'] = Psurf[:, 0]

            if hasattr(output,'Etotal'):
                source['Etotal'] = output.Etotal[:, :, 0]
                source['Entropy'] = output.Entropy[:, :, 0]
                source['AngMomz'] = output.AngMomz[:, :, 0]

            if input.TSRT:
                source['mustar'] = output.mustar[:, 0]
                source['TSf_net'] = output.f_net[:, :-1, 0] + (output.f_net[:, 1:, 0] - output.f_net[:, :-1, 0]) * interpz[None, :]
                source['TSqheat'] = output.TSqheat[:, :, 0]
                source['f_up_tot'] = output.f_up_tot[:, :-1, 0] + (output.f_up_tot[:, 1:, 0] - output.f_up_tot[:, :-1, 0]) * interpz[None, :]
                source['f_down_tot'] = output.f_down_tot[:, :-1, 0] + (output.f_down_tot[:, 1:, 0] - output.f_down_tot[:, :-1, 0]) * interpz[None, :]
            if chem == 1:
                source['ch4'] = output.ch4[:, :, 0] / output.Rho[:,:,0]
                source['co'] = output.co[:, :, 0]/ output.Rho[:,:,0]
                source['h2o'] = output.h2o[:, :, 0]/ output.Rho[:,:,0]
                source['co2'] = output.co2[:, :, 0]/ output.Rho[:,:,0]
                source['nh3'] = output.nh3[:, :, 0]/ output.Rho[:,:,0]

            # calculate zonal and meridional velocity (special step for Mh)
            source['U'] = (-source['Mh'][0]*np.sin(grid.lon[:,None])+\
                           source['Mh'][1]*np.cos(grid.lon[:, None]))/source['Rho']
            source['V'] = (-source['Mh'][0]*np.sin(grid.lat[:,None])*np.cos(grid.lon[:,None])\
                      -source['Mh'][1]*np.sin(grid.lat[:,None])*np.sin(grid.lon[:,None])\
                           + source['Mh'][2]*np.cos(grid.lat[:, None]))/source['Rho']
            source['U_mean'] = (-source['Mh_mean'][0]*np.sin(grid.lon[:,None])+\
                           source['Mh_mean'][1]*np.cos(grid.lon[:, None]))/source['Rho_mean']
            source['V_mean'] = (-source['Mh_mean'][0]*np.sin(grid.lat[:,None])*np.cos(grid.lon[:,None])\
                      -source['Mh_mean'][1]*np.sin(grid.lat[:,None])*np.sin(grid.lon[:,None])\
                           + source['Mh_mean'][2]*np.cos(grid.lat[:, None]))/source['Rho_mean']

            # set up intermediate arrays (icogrid and pressure)
            interm = {}
            for key in source.keys():
                if key == 'Mh' or key == 'Mh_mean':
                    pass
                elif np.shape(source[key]) == (grid.point_num,):
                    # 2D field (e.g., insolation) -> not needed
                    interm[key] = np.zeros((d_lon[0], d_lon[1]))
                else:
                    interm[key] = np.zeros((d_lon[0], d_lon[1], grid.nv))

            for key in interm.keys():
                if np.shape(interm[key]) == (d_lon[0], d_lon[1]):
                    # 2D field (e.g., insolation)
                    tmp = np.sum(weight3[:, :] * source[key][near3], axis=1)
                    interm[key][:, :] = tmp.reshape((d_lon[0], d_lon[1]))
                else:
                    tmp = np.sum(weight3[:, :, None] * source[key][near3], axis=1)
                    interm[key][:, :] = tmp.reshape((d_lon[0], d_lon[1], grid.nv))

            # dealing with RV and PV
            # potential temp (only for constant Rd and Cp, currently)
            pt = interm['Temperature'] * (interm['Pressure'] / input.P_Ref)**(-interm['Rd'] / interm['Cp'])
            dptdz = np.gradient(pt, grid.Altitude, axis=2)
            # these are not very accurate at the 'edges'
            dptdlat = np.gradient(pt, lat_range, axis=0) / (input.A + grid.Altitude[None, None, :])
            dptdlon = np.gradient(pt, lon_range, axis=1) / (input.A + grid.Altitude[None, None, :]) / np.cos(lat_range[:, None, None])

            # relative vorticity
            curlVz, curlVlat, curlVlon = CurlF(interm['W'], interm['V'], interm['U'],
                                               lat_range, lon_range, grid.Altitude, input.A)
            interm['RVw'] = curlVz  # vertical component
            interm['RVv'] = curlVlat  # horizontal components
            interm['RVu'] = curlVlon  # these are probably small in most cases

            curlVlon += 2 * input.Omega  # add solid body rotation
            interm['PV'] = (curlVz * dptdz + curlVlat * dptdlat + curlVlon * dptdlon) / interm['Rho']

            # set up destination arrays (lat-lon and pressure/height)
            dest = {}
            for key in interm.keys():
                if key == 'Mh' or key == 'Mh_mean' or key == 'Pressure':
                    pass  # don't need these any further
                elif np.shape(interm[key]) == (d_lon[0], d_lon[1]):
                    # 2D field (e.g., insolation)
                    dest[key] = np.zeros((d_lon[0], d_lon[1]))
                else:
                    dest[key] = np.zeros((d_lon[0], d_lon[1], d_sig))

            # regrid to pressure coordinate
            # need to make sure the arrays are contiguous
            x = np.ascontiguousarray(interm['Pressure'][:, :, ::-1].ravel())
            xnew = np.ascontiguousarray(Pref[::-1])

            for key in dest.keys():
                if np.shape(interm[key]) == (d_lon[0], d_lon[1]):
                    # 2D field (e.g., insolation)
                    dest[key] = interm[key]

                else:
                    dest[key][:, :, :] = vertical_regrid_field(interm[key][:, :, ::-1], grid.nv, x, xnew)[:, :, ::-1]
                    if surf and mask_surf:
                        dest[key][(Pref[None, None, :] >
                                   interm['Psurf'][:, :, None])] = np.nan

            # create h5 files (pressure grid)
            openh5 = h5py.File(fileh5p, "w")
            print('Writing file ' + fileh5p + '...')
            # coordinates
            Pre = openh5.create_dataset("Pressure", data=Pref,
                                        compression='gzip', compression_opts=comp)
            Lat = openh5.create_dataset("Latitude", data=lat_range * 180 / np.pi,
                                        compression='gzip', compression_opts=comp)
            Lon = openh5.create_dataset("Longitude", data=lon_range * 180 / np.pi,
                                        compression='gzip', compression_opts=comp)
            # data
            for key in dest.keys():
                tmp_data = openh5.create_dataset(key, data=dest[key],
                                                 compression='gzip', compression_opts=comp)
            openh5.close()

            # create h5 files (height grid)
            if not skip_height_file:
                openh5 = h5py.File(fileh5h, "w")
                print('Writing file ' + fileh5h + '...')
                Alt = openh5.create_dataset("Altitude", data=grid.Altitude,
                                            compression='gzip', compression_opts=comp)
                Lat = openh5.create_dataset("Latitude", data=lat_range * 180 / np.pi,
                                            compression='gzip', compression_opts=comp)
                Lon = openh5.create_dataset("Longitude", data=lon_range * 180 / np.pi,
                                            compression='gzip', compression_opts=comp)
                # data
                for key in dest.keys():
                    tmp_data = openh5.create_dataset(key, data=interm[key],
                                                     compression='gzip', compression_opts=comp)
                openh5.close()


def KE_spect(input, grid, output, sigmaref, coord='icoh', lmax_adjust=0):
    tsp = output.nts - output.ntsi + 1
    lmax_grid = np.int(np.floor(np.sqrt(grid.point_num)) / 2 - 1)

    if coord == 'icoh':
        W = 0.5 * (output.Wh[:, 1:, :] + output.Wh[:, :-1, :])
        Wx = W * np.cos(grid.lat[:, None, None]) * np.cos(grid.lon[:, None, None])
        Wy = W * np.cos(grid.lat[:, None, None]) * np.sin(grid.lon[:, None, None])
        Wz = W * np.sin(grid.lat[:, None, None])

        Vx = (output.Mh[0] + Wx)/output.Rho
        Vy = (output.Mh[1] + Wy)/output.Rho
        Vz = (output.Mh[2] + Wz)/output.Rho

        KE = 0.5 * (Vx**2 + Vy**2 + Vz**2)
        lmax = np.int(lmax_grid + lmax_adjust)  # sets lmax based on grid size

        x_coeffs = np.zeros((2, lmax + 1, lmax + 1, grid.nv, tsp), dtype=complex)
        y_coeffs = np.zeros((2, lmax + 1, lmax + 1, grid.nv, tsp), dtype=complex)
        z_coeffs = np.zeros((2, lmax + 1, lmax + 1, grid.nv, tsp), dtype=complex)
        KE_coeffs = np.zeros((2, lmax + 1, lmax + 1, grid.nv, tsp), dtype=complex)
        KE_power = np.zeros((lmax + 1, grid.nv, tsp))
        waven = np.arange(lmax + 1)  # total spherical wavenumber

        cmap = cm.get_cmap('cividis')
        fig, ax = plt.subplots(1, 1)

        if tsp == 1:
            for lev in np.arange(grid.nv):
                KE_coeffs[:, :, :, lev, 0], chiz = chairs.expand.SHExpandLSQ(KE[:, lev, 0], grid.lat * 180 / np.pi, grid.lon * 180 / np.pi, lmax)
                KE_power[:, lev, 0] = chairs.spectralanalysis.spectrum(KE_coeffs[:,:,:,lev,0], unit='per_lm')
                ax.plot(waven, KE_power[:, lev, 0], 'k-', c=cmap(lev / grid.nv), lw=1)
        else:
            for t in np.arange(tsp):
                for lev in np.arange(grid.nv):
                    KE_coeffs[:, :, :, lev, t], chiz = chairs.expand.SHExpandLSQ(KE[:, lev, t], grid.lat * 180 / np.pi, grid.lon * 180 / np.pi, lmax)
                    KE_power[:, lev, t] = chairs.spectralanalysis.spectrum(KE_coeffs[:,:,:,lev,t], unit='per_lm')
                    # ax.plot(waven, KE_power[:, lev, t], 'k-', c=cmap(lev / grid.nv), lw=1)

            KE_power_mean = np.mean(KE_power,axis=2)
            for lev in np.arange(grid.nv):
                ax.plot(waven, KE_power_mean[:, lev], 'k-', c=cmap(lev/grid.nv), lw=1)

    else:
        raise IOError("Invalid coord option! Valid options are 'icoh'")

    # ax.plot(waven,np.mean(KE_power[:,:,t],axis=1),'k-',lw=2)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.vlines(lmax_grid, ax.get_ylim()[0], ax.get_ylim()[1], zorder=1000, linestyle='--')
    ax.set(ylabel='KE (m$^2$ s$^{-2}$)', xlabel='Spherical wavenumber $n$')

    if not os.path.exists(input.resultsf + '/figures'):
        os.mkdir(input.resultsf + '/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf + '/figures/KEspectrum_%i_%i_%s.pdf' % (output.ntsi, output.nts, coord))
    plt.close()

    lev = 0
    norm = colors.Normalize(vmin=np.min(KE[:, lev, 0]), vmax=np.max(KE[:, lev, 0]))

    fig = plt.figure()
    fig.subplots_adjust(left=0.1, right=0.97)
    plt.subplot(2, 1, 1)
    # if coord == 'llp':  #not complete
    #     plt.imshow(KE[:, :, 0, 0], origin='lower', extent=(0, 360, -90, 90), norm=norm, aspect='auto')
    if coord == 'icoh':
        cont = plt.tricontourf(grid.lon * 180 / np.pi, grid.lat * 180 / np.pi, KE[:, lev, 0], levels=30)
        for cc in cont.collections:
            cc.set_edgecolor("face")  # fixes a stupid bug in matplotlib 2.0
    plt.xlabel('Longitude ($^{\circ}$)')
    plt.ylabel('Latitude ($^{\circ}$)')
    plt.title("Model output, lowest level")
    clb = plt.colorbar()
    clb.set_label('Kinetic energy (m$^2$ s$^{-2}$)')

    KEcomp = chairs.expand.MakeGridDH(KE_coeffs[:, :, :, lev, 0], sampling=2)
    lat = np.linspace(-90, 90, np.shape(KEcomp)[0])
    lon = np.linspace(0, 360, np.shape(KEcomp)[1])

    plt.subplot(2, 1, 2)
    plt.imshow(np.real(KEcomp), origin='upper', extent=(0, 360, -90, 90), norm=norm, aspect='auto')
    plt.xlabel('Latitude ($^{\circ}$)')
    plt.ylabel('Longitude ($^{\circ}$)')
    plt.title("Spherical harmonics reconstruction, lowest level")
    clb = plt.colorbar()
    clb.set_label('Kinetic energy (m$^{2}$ s$^{-2}$)')

    plt.tight_layout()
    plt.savefig(input.resultsf + '/figures/KEmap_lowest_%i_%s.pdf' % (output.ntsi, coord))
    plt.close()


def maketable(x, y, z, xname, yname, zname, resultsf, fname):
    # print 2D data from plot to file
    if not os.path.exists(resultsf + '/tables'):
        os.mkdir(resultsf + '/tables')
    f = open(resultsf + '/tables/' + fname, 'w')
    f.write('# %s %s %s\n' % (xname, yname, zname))
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            f.write('%#.6e %#.6e %#.6e\n' % (x[i], y[j], z[i, j]))
    f.close()


def vertical_lat(input, grid, output, rg, sigmaref, z, slice=['default'], save=True, axis=None, csp=500, wind_vectors=False, use_p=True, clevs=[40]):
    # generic pressure/latitude plot function

    # Set the reference pressure
    if use_p:
        # sigmaref = normalize pressure units
        Pref = input.P_Ref * sigmaref

    d_sig = np.size(sigmaref)
    tsp = output.nts - output.ntsi + 1

    if not isinstance(slice, list):
        raise IOError("'slice' argument must be a list")
    if slice[0] == 'default':
        slice = [0, 360]
    else:
        slice = [np.float(slice[0]),np.float(slice[1])]

    lat = z['lat']
    lon = z['lon']
    if len(slice) == 2:
        # Set the latitude-longitude grid
        if slice[1] - slice[0] > 360:
            raise IOError("'slice' cannot exceed a range of 360 degrees")
        if slice[1] < slice[0]:
            raise IOError("'slice' values must be arranged (small, large) because I am too lazy to code it better")
        if slice[0] < 0:
            mask_ind = np.logical_or((lon - 360) >= slice[0], lon <= slice[1])
        else:
            mask_ind = np.logical_and(lon >= slice[0], lon <= slice[1])
        # loni, lati = np.meshgrid(lon[mask_ind],lat)
        # d_lon = np.shape(loni)

        ##################
        #    Averages    #
        ##################
        # Averaging in time and longitude
        if tsp > 1:
            Zonall = np.nanmean(z['value'][:, mask_ind[:], :, :], axis=1)
            # Vl = np.nanmean(rg.V[:,:,:,:],axis=1)
            # Wl = np.nanmean(rg.W[:,:,:,:],axis=1)
            Zonallt = np.nanmean(Zonall[:, :, :], axis=2)
            # Vlt = np.nanmean(Vl[:,:,:],axis=2)
            # Wlt = np.nanmean(Wl[:,:,:],axis=2)
            del Zonall
            if wind_vectors == True:
                Vl = np.nanmean(rg.V[:, mask_ind[:], :, :], axis=1)
                Wl = np.nanmean(rg.W[:, mask_ind[:], :, :], axis=1)
                Vlt = np.nanmean(Vl[:, :, :], axis=2)
                Wlt = np.nanmean(Wl[:, :, :], axis=2)
                del Vl, Wl
        else:
            Zonallt = np.nanmean(z['value'][:, :, :, 0][:, mask_ind[:], :], axis=1)
            # Vlt = np.nanmean(rg.V[:,:,:,0],axis=1)
            # Wlt = np.nanmean(rg.W[:,:,:,0],axis=1)
            if wind_vectors == True:
                Vlt = np.nanmean(rg.V[:, :, :, 0][:, mask_ind[:], :], axis=1)
                Wlt = np.nanmean(rg.W[:, :, :, 0][:, mask_ind[:], :], axis=1)

    elif len(slice) == 1:
        if slice[0] in lon:
            Zonall = z['value'][:, lon[:] == slice[0], :, :]
            if wind_vectors == True:
                Vl = rg.V[:, lon[:] == slice[0], :, :]
                Wl = rg.W[:, lon[:] == slice[0], :, :]
        else:
            Zonall = np.zeros((len(lat), 1, d_sig, tsp))
            if wind_vectors == True:
                Vl = np.zeros((len(lat), 1, d_sig, tsp))
                Wl = np.zeros((len(lat), 1, d_sig, tsp))
            # interpolate to slice given
            for t in tsp:
                for lev in np.arange(d_sig):
                    Zonall[:, 0, lev, tsp] = interp.griddata(np.vstack([lon, lat]).T, z['value'][:, :, lev, tsp], (slice[0], lat))
                    if wind_vectors == True:
                        Vl[:, 0, lev, tsp] = interp.griddata(np.vstack([lon, lat]).T, rg.V[:, :, lev, tsp], (slice[0], lat))
                        Wl[:, 0, lev, tsp] = interp.griddata(np.vstack([lon, lat]).T, rg.W[:, :, lev, tsp], (slice[0], lat))

        # Averaging in time
        if tsp > 1:
            Zonallt = np.nanmean(Zonall[:, 0, :, :], axis=2)
            del Zonall
            if wind_vectors == True:
                Vlt = np.nanmean(Vl[:, 0, :, :], axis=2)
                Wlt = np.nanmean(Wl[:, 0, :, :], axis=2)
                del Vl, Wl
        else:
            Zonallt = Zonall[:, 0, :, 0]
            if wind_vectors == True:
                Vlt = Vl[:, 0, :, 0]
                Wlt = Wl[:, 0, :, 0]

    else:
        raise IOError("'slice' must have 1 or 2 values")

    #################
    # Create figure #
    #################
    # Latitude
    latp = lat[:] * np.pi / 180

    if use_p:
        # need to set desired pressure range (major PITA!)
        # prange = np.where(np.logical_and(rg.Pressure[:,0]>=np.min(Pref),rg.Pressure[:,0]<=np.max(Pref)))
        prange = np.where(rg.Pressure[:] >= np.min(Pref))
        ycoord = rg.Pressure[prange[0]] / 1e5
        zvals = Zonallt[:, prange[0]].T
    else:
        hrange = np.where(np.logical_and(rg.Altitude[:] >= np.min(sigmaref), rg.Altitude[:] <= np.max(sigmaref)))
        ycoord = rg.Altitude[hrange[0]] / 1000
        zvals = Zonallt[:, hrange[0]].T

    # Contour plot
    if len(clevs) == 1:
        clevels = np.int(clevs[0])
    elif len(clevs) == 3:
        if isinstance(clevs[2],str) and 'log' in clevs[2]:
            clevels = np.logspace(np.log10(np.float(clevs[0])),np.log10(np.float(clevs[1])),np.int(clevs[2][:-3]))
        else:
            clevels = np.linspace(np.int(clevs[0]),np.int(clevs[1]),np.int(clevs[2]))
    else:
        raise IOError("clevs not valid!")
    # print(np.max(zvals))

    if isinstance(axis, axes.SubplotBase):
        ax = axis
        fig = plt.gcf()
    elif (isinstance(axis, tuple)
          and isinstance(axis[0], plt.Figure)
          and isinstance(axis[1], axes.SubplotBase)):
        fig = axis[0]
        ax = axis[1]
    elif axis is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    else:
        raise IOError("'axis = {}' but {} is neither an axes.SubplotBase instance nor a (axes.SubplotBase, plt.Figure) instance".format(axis, axis))

    fig.set_tight_layout(True)

    C = ax.contourf(latp * 180 / np.pi, ycoord, zvals, clevels, cmap=z['cmap'])

    if wind_vectors == True:
        vspacing = np.int(np.shape(rg.Latitude)[0] / 10)
        if use_p:
            wspacing = np.int(np.shape(rg.Pressure)[0] / 10)
            Vlt = Vlt[:, prange[0]]
            Wlt = Wlt[:, prange[0]]
            yqcoord = rg.Pressure[::wspacing][prange[0]]
        else:
            wspacing = np.int(np.shape(rg.Altitude)[0] / 10)
            Vlt = Vlt[:, hrange[0]]
            Wlt = Wlt[:, hrange[0]]
            yqcoord = rg.Altitude[::wspacing][hrange[0]]

        Vq = Vlt[::vspacing, ::wspacing].ravel()
        Wq = Wlt[::vspacing, ::wspacing].ravel()
        #preq = rg.Pressure[:,0][::spacing,::spacing].ravel()
        #latq = lati[::spacing,::spacing].ravel()
        latq, preq = np.meshgrid(rg.Latitude[::vspacing], yqcoord)
        del Vlt, Wlt
        ax.quiver(latq.ravel(), preq.ravel() / 1e5, Vq, Wq, color='0.5')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    kwargs = {'format': '%#3.4e'}
    #clb = plt.colorbar(C, extend='both', ax=ax)
    clb = fig.colorbar(C, cax=cax, extend='both', **kwargs)
    clb.set_label(z['label'])

    if isinstance(csp, list) or isinstance(csp, tuple):
        levp = csp
    else:
        if csp == 'match':
            levp = 40
        else:
            if use_p:
                levp = np.arange(np.ceil(np.nanmin(Zonallt[:, prange[0]]) / csp) * csp, np.floor(np.nanmax(Zonallt[:, prange[0]]) / csp) * csp, csp)
            else:
                levp = np.arange(np.ceil(np.nanmin(Zonallt[:, hrange[0]]) / csp) * csp, np.floor(np.nanmax(Zonallt[:, hrange[0]]) / csp) * csp, csp)

    c2 = ax.contour(latp * 180 / np.pi, ycoord, zvals, levels=levp, colors='w', linewidths=1)
    ax.clabel(c2, inline=False, fontsize=6, fmt='%d', use_clabeltext=True)
    for cc in C.collections:
        cc.set_edgecolor("face")  # fixes a stupid bug in matplotlib 2.0
    # ax.invert_yaxis()
    # plt.quiver(latq.ravel(),preq.ravel()/1e5,Vq/np.max(Vq),Wq/np.max(Wq),color='0.5')
    if z['plog'] == True:
        ax.set_yscale("log")
    ax.set_xlabel('Latitude (deg)')
    if use_p:
        ax.set_ylabel('Pressure (bar)')
        # ax.plot(latp*180/np.pi,np.zeros_like(latp)+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
        #ax.set_ylim(np.max(rg.Pressure[prange[0], 0])/1e5, np.min(Pref)/1e5)
        ax.set_ylim(np.max(rg.Pressure[prange[0]]) / 1e5, np.min(rg.Pressure[prange[0]]) / 1e5)

        if ax.get_ylim()[1] > ax.get_ylim()[0]:
            ax.invert_yaxis()
    else:
        ax.set_ylabel('Altitude (km)')

    if len(slice) == 2:
        ax.set_title('Time = %#.3f-%#.3f days, Lon = (%#.3f,%#.3f)' % (output.time[0], output.time[-1], slice[0], slice[1]), fontsize=10)
    else:
        ax.set_title('Time = %#.3f-%#.3f days, Lon = (%#.3f,)' % (output.time[0], output.time[-1], slice[0]), fontsize=10)

    if use_p:
        z['name'] += '_p'
    else:
        z['name'] += '_h'

    pfile = None
    if save == True:
        output_path = pathlib.Path(input.resultsf) / 'figures'
        if not output_path.exists():
            output_path.mkdir()

        # save the plot to file designated by z
        if len(slice) == 2:
            fname = '%s_ver_i%d_l%d_lon%#.2f-%#.2f' % (z['name'], output.ntsi, output.nts, slice[0], slice[1])
        else:
            fname = '%s_ver_i%d_l%d_lon%#.2f' % (z['name'], output.ntsi, output.nts, slice[0])

        pfile = output_path / (fname.replace(".", "+") + '.pdf')
        plt.savefig(pfile)

        plt.close()
        pfile = str(pfile)

    if z['mt'] == True:
        if len(slice) == 2:
            fname = '%s_ver_i%d_l%d_lon%#.2f-%#.2f.dat' % (z['name'], output.ntsi, output.nts, slice[0], slice[1])
        else:
            fname = '%s_ver_i%d_l%d_lon%#.2f.dat' % (z['name'], output.ntsi, output.nts, slice[0])
        if use_p:
            maketable(latp * 180 / np.pi, rg.Pressure[prange[0]] / 1e5, Zonallt[:, prange[0]], 'Latitude(d)', 'Pressure(bar)', z['name'], input.resultsf, fname)
        else:
            maketable(latp * 180 / np.pi, rg.Altitude[hrange[0]], Zonallt[:, hrange[0]], 'Latitude(d)', 'Altitude(m)', z['name'], input.resultsf, fname)

    return pfile


def vertical_lon(input, grid, output, rg, sigmaref, z, slice='default', save=True, axis=None, csp=500, wind_vectors=False, use_p=True, clevs=[40]):
    # generic pressure/longitude plot function

    # Set the reference pressure
    if use_p:
        # sigmaref = normalize pressure units
        Pref = input.P_Ref * sigmaref

    d_sig = np.size(sigmaref)
    tsp = output.nts - output.ntsi + 1

    if not isinstance(slice, list):
        raise IOError("'slice' argument must be a list")
    if slice[0] == 'default':
        slice = [-90, 90]
    else:
        slice = [np.float(slice[0]),np.float(slice[1])]

    lat = z['lat']
    lon = z['lon']
    if len(slice) == 2:
        # Set the latitude-longitude grid
        if slice[1] - slice[0] > 180:
            raise IOError("'slice' cannot exceed a range of 180 degrees")
        if slice[1] < slice[0]:
            raise IOError("'slice' values must be arranged (small, large) because I am too lazy to code it better")
        # if slice[0] < 0:
        #     mask_ind = np.logical_or((lon-360)>=slice[0],lon<=slice[1])
        # else:
        #     mask_ind = np.logical_and(lon>=slice[0],lon<=slice[1])
        mask_ind = np.logical_and(lat >= slice[0], lat <= slice[1])
        # loni, lati = np.meshgrid(lon,lat[mask_ind])
        # d_lat = np.shape(lati)

        ##################
        #    Averages    #
        ##################

        # Averaging in time and latitude (weighted by a cosine(lat))
        if tsp > 1:
            Meridl = np.nanmean(z['value'][mask_ind[:], :, :, :] * np.cos(lat[mask_ind[:]] * np.pi / 180)[:, None, None, None], axis=0)
            # Vl = np.nanmean(rg.V[:,:,:,:],axis=1)
            # Wl = np.nanmean(rg.W[:,:,:,:],axis=1)
            Meridlt = np.nanmean(Meridl[:, :, :], axis=2)
            # Vlt = np.nanmean(Vl[:,:,:],axis=2)
            # Wlt = np.nanmean(Wl[:,:,:],axis=2)
            del Meridl
            if wind_vectors == True:
                Ul = np.nanmean(rg.U[mask_ind[:], :, :, :] * np.cos(lat[mask_ind[:], 0] * np.pi / 180)[:, None, None, None], axis=0)
                Wl = np.nanmean(rg.W[mask_ind[:], :, :, :] * np.cos(lat[mask_ind[:], 0] * np.pi / 180)[:, None, None, None], axis=0)
                Ult = np.nanmean(Ul[:, :, :], axis=2)
                Wlt = np.nanmean(Wl[:, :, :], axis=2)
                del Ul, Wl
        else:
            Meridlt = np.nanmean(z['value'][:, :, :, 0][mask_ind[:], :, :] * np.cos(lat[mask_ind[:]] * np.pi / 180)[:, None, None], axis=0)
            # Ult = np.nanmean(rg.V[:,:,:,0],axis=1)
            # Wlt = np.nanmean(rg.W[:,:,:,0],axis=1)
            if wind_vectors == True:
                Ult = np.nanmean(rg.U[:, :, :, 0][mask_ind[:], :, :] * np.cos(lat[mask_ind[:]] * np.pi / 180)[:, None, None], axis=0)
                Wlt = np.nanmean(rg.W[:, :, :, 0][mask_ind[:], :, :] * np.cos(lat[mask_ind[:]] * np.pi / 180)[:, None, None], axis=0)

    elif len(slice) == 1:
        if slice[0] in lat:
            Meridl = z['value'][lat[:] == slice[0], :, :, :]
            if wind_vectors == True:
                Ul = rg.U[lat[:] == slice[0], :, :, :]
                Wl = rg.W[lat[:] == slice[0], :, :, :]
        else:
            Meridl = np.zeros((1, len(lon), d_sig, tsp))
            if wind_vectors == True:
                Ul = np.zeros((1, len(lon), d_sig, tsp))
                Wl = np.zeros((1, len(lon), d_sig, tsp))
            # interpolate to slice given
            for t in tsp:
                for lev in np.arange(d_sig):
                    Meridl[0, :, lev, tsp] = interp.griddata(np.vstack([lon, lat]).T, z['value'][:, :, lev, tsp], (lon, slice[0]))
                    if wind_vectors == True:
                        Ul[0, :, lev, tsp] = interp.griddata(np.vstack([lon, lat]).T, rg.U[:, :, lev, tsp], (lon, slice[0]))
                        Wl[0, :, lev, tsp] = interp.griddata(np.vstack([lon, lat]).T, rg.W[:, :, lev, tsp], (lon, slice[0]))

        # Averaging in time
        if tsp > 1:
            Meridlt = np.nanmean(Meridl[0, :, :, :], axis=2)
            del Meridl
            if wind_vectors == True:
                Ult = np.nanmean(Ul[0, :, :, :], axis=2)
                Wlt = np.nanmean(Wl[0, :, :, :], axis=2)
                del Ul, Wl
        else:
            Meridlt = Meridl[0, :, :, 0]
            if wind_vectors == True:
                Ult = Ul[0, :, :, 0]
                Wlt = Wl[0, :, :, 0]

    else:
        raise IOError("'slice' must have 1 or 2 values")

    #################
    # Create figure #
    #################
    # Longitude
    lonp = lon[:] * np.pi / 180

    if use_p:
        # need to set desired pressure range (major PITA!)
        # prange = np.where(np.logical_and(rg.Pressure[:,0]>=np.min(Pref),rg.Pressure[:,0]<=np.max(Pref)))
        prange = np.where(rg.Pressure[:] >= np.min(Pref))
        ycoord = rg.Pressure[prange[0]] / 1e5
        zvals = Meridlt[:, prange[0]].T
    else:
        hrange = np.where(np.logical_and(rg.Altitude[:] >= np.min(sigmaref), rg.Altitude[:] <= np.max(sigmaref)))
        ycoord = rg.Altitude[hrange[0]] / 1000
        zvals = Meridlt[:, hrange[0]].T

    # Contour plot
    if len(clevs) == 1:
        clevels = np.int(clevs[0])
    elif len(clevs) == 3:
        clevels = np.linspace(np.int(clevs[0]), np.int(clevs[1]), np.int(clevs[2]))
    else:
        raise IOError("clevs not valid!")
    # print(np.min(zvals),np.max(zvals))
    if isinstance(axis, axes.SubplotBase):
        ax = axis
        fig = plt.gcf()
    elif (isinstance(axis, tuple)
          and isinstance(axis[0], plt.Figure)
          and isinstance(axis[1], axes.SubplotBase)):
        fig = axis[0]
        ax = axis[1]
    elif axis is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    else:
        raise IOError("'axis = {}' but {} is neither an axes.SubplotBase instance nor a (axes.SubplotBase, plt.Figure) instance".format(axis, axis))

    fig.set_tight_layout(True)

    C = ax.contourf(lonp * 180 / np.pi, ycoord, zvals, clevels, cmap=z['cmap'])

    if wind_vectors == True:
        vspacing = np.int(np.shape(rg.Longitude)[0] / 10)
        if use_p:
            wspacing = np.int(np.shape(rg.Pressure)[0] / 10)
            Ult = Ult[:, prange[0]]
            Wlt = Wlt[:, prange[0]]
            yqcoord = rg.Pressure[::wspacing, 0][prange[0]]
        else:
            wspacing = np.int(np.shape(rg.Altitude)[0] / 10)
            Ult = Ult[:, hrange[0]]
            Wlt = Wlt[:, hrange[0]]
            yqcoord = rg.Altitude[::wspacing][hrange[0]]

        Uq = Ult[::vspacing, ::wspacing].ravel()
        Wq = Wlt[::vspacing, ::wspacing].ravel()
        #preq = rg.Pressure[:,0][::spacing,::spacing].ravel()
        #latq = lati[::spacing,::spacing].ravel()
        lonq, preq = np.meshgrid(rg.Longitude[::vspacing], yqcoord)
        del Ult, Wlt
        ax.quiver(lonq.ravel(), preq.ravel() / 1e5, Uq, Wq, color='0.5')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    kwargs = {'format': '%#3.4e'}
    #clb = plt.colorbar(C, extend='both', ax=ax)
    clb = fig.colorbar(C, cax=cax, extend='both', **kwargs)
    clb.set_label(z['label'])

    if isinstance(csp, list) or isinstance(csp, tuple):
        levp = csp
    else:
        if csp == 'match':
            levp = 40
        else:
            if use_p:
                levp = np.arange(np.ceil(np.nanmin(Meridlt[:, prange[0]]) / csp) * csp, np.floor(np.nanmax(Meridlt[:, prange[0]]) / csp) * csp, csp)
            else:
                levp = np.arange(np.ceil(np.nanmin(Meridlt[:, hrange[0]]) / csp) * csp, np.floor(np.nanmax(Meridlt[:, hrange[0]]) / csp) * csp, csp)

    c2 = ax.contour(lonp * 180 / np.pi, ycoord, zvals, levels=levp, colors='w', linewidths=1)
    ax.clabel(c2, inline=False, fontsize=6, fmt='%d', use_clabeltext=True)

    for cc in C.collections:
        cc.set_edgecolor("face")  # fixes a stupid bug in matplotlib 2.0
    # ax.invert_yaxis()
    # plt.quiver(latq.ravel(),preq.ravel()/1e5,Vq/np.max(Vq),Wq/np.max(Wq),color='0.5')
    if z['plog'] == True:
        ax.set_yscale("log")
    ax.set_xlabel('Longitude (deg)')
    if use_p:
        ax.set_ylabel('Pressure (bar)')
        #ax.plot(lonp * 180 / np.pi, np.zeros_like(lonp) + np.max(output.Pressure[:, grid.nv - 1, :]) / 1e5, 'r--')
        #ax.set_ylim(np.max(rg.Pressure[prange[0]]) / 1e5, np.min(Pref) / 1e5)
        ax.set_ylim(np.max(rg.Pressure[prange[0]]) / 1e5, np.min(rg.Pressure[prange[0]]) / 1e5)

        if ax.get_ylim()[1] > ax.get_ylim()[0]:
            ax.invert_yaxis()
    else:
        ax.set_ylabel('Altitude (km)')

    if len(slice) == 2:
        ax.set_title('Time = %#.3f-%#.3f days, Lat = (%#.3f,%#.3f)' % (output.time[0], output.time[-1], slice[0], slice[1]), fontsize=10)
    else:
        ax.set_title('Time = %#.3f-%#.3f days, Lat = (%#.3f,)' % (output.time[0], output.time[-1], slice[0]), fontsize=10)

    if not os.path.exists(input.resultsf + '/figures'):
        os.mkdir(input.resultsf + '/figures')

    if use_p:
        z['name'] += '_p'
    else:
        z['name'] += '_h'
    pfile = None
    if save == True:
        # save the plot to file designated by z
        if len(slice) == 2:
            fname = '%s_ver_i%d_l%d_lat%#.2f-%#.2f' % (z['name'], output.ntsi, output.nts, slice[0], slice[1])
            pfile = input.resultsf + '/figures/' + fname.replace(".", "+") + '.pdf'
        else:
            fname = '%s_ver_i%d_l%d_lat%#.2f' % (z['name'], output.ntsi, output.nts, slice[0])
            pfile = input.resultsf + '/figures/' + fname.replace(".", "+") + '.pdf'

        plt.savefig(pfile)

        plt.close()

        pfile = str(pfile)

    if z['mt'] == True:
        if len(slice) == 2:
            fname = '%s_ver_i%d_l%d_lat%#.2f-%#.2f.dat' % (z['name'], output.ntsi, output.nts, slice[0], slice[1])
        else:
            fname = '%s_ver_i%d_l%d_lat%#.2f.dat' % (z['name'], output.ntsi, output.nts, slice[0])
        if use_p:
            maketable(lonp * 180 / np.pi, rg.Pressure[prange[0]] / 1e5, Meridlt[:, prange[0]], 'Longitude(d)', 'Pressure(bar)', z['name'], input.resultsf, fname)
        else:
            maketable(lonp * 180 / np.pi, rg.Altitude[hrange[0]], Meridlt[:, hrange[0]],
                      'Longitude(d)', 'Altitude(m)', z['name'], input.resultsf, fname)
    return pfile


def horizontal_lev(input, grid, output, rg, Plev, z, save=True, axis=False, wind_vectors=False, use_p=True, clevs=[40]):
    # Set the latitude-longitude grid.
    loni, lati = np.meshgrid(rg.Longitude[:], rg.Latitude[:])

    d_lon = np.shape(loni)
    tsp = output.nts - output.ntsi + 1

    zlev = np.zeros(np.shape(loni) + (tsp,))
    if wind_vectors == True:
        Uii = np.zeros(np.shape(loni) + (tsp,))
        Vii = np.zeros(np.shape(loni) + (tsp,))

    if use_p:
        vcoord = rg.Pressure
    else:
        vcoord = rg.Altitude

    for t in np.arange(tsp):
        # interpolate
        # above is index of higher pressure bound (lower index)
        # below is index of lower pressure bound (higher index)
        if np.size(vcoord[vcoord > Plev]) == 0:
            above = np.where(vcoord == np.max(vcoord))[0][0]
            if use_p:
                below = above + 1
            else:
                below = above - 1
        elif np.size(vcoord[vcoord < Plev]) == 0:
            below = np.where(vcoord == np.min(vcoord))[0][0]
            if use_p:
                above = below - 1
            else:
                above = below + 1
        else:
            above = np.where(vcoord == np.min(vcoord[vcoord > Plev]))[0][0]
            below = np.where(vcoord == np.max(vcoord[vcoord < Plev]))[0][0]

        if len(np.shape(z['value'])) == 4:
            # single level of a 3D field
            zlev[:, :, t] = (z['value'][:, :, below, t] * (vcoord[above] - Plev)
                             + z['value'][:, :, above, t] * (Plev - vcoord[below]))\
                             / (vcoord[above] - vcoord[below])
            if use_p:
                title = 'Time = %#.3f-%#.3f days, Plev = %#.3f bar' % (output.time[0], output.time[-1], Plev / 1e5)
                fname = '%s_lev%#.3fmbar_i%d_l%d' % (z['name'], Plev / 100, output.ntsi, output.nts)
            else:
                title = 'Time = %#.3f-%#.3f days, lev = %#.3f km' % (output.time[0], output.time[-1], Plev / 1e3)
                fname = '%s_lev%#.3fkm_i%d_l%d' % (z['name'], Plev / 1000, output.ntsi, output.nts)
        elif len(np.shape(z['value'])) == 3:
            # 2D field, only one level available (e.g., insolation or surface temp)
            zlev[:, :, t] = z['value'][:, :, t]
            title = 'Time = %#.3f-%#.3f days' % (output.time[0], output.time[-1])
            fname = '%s_i%d_l%d' % (z['name'], output.ntsi, output.nts)
            wind_vectors = False

        if wind_vectors == True:
            Uii[:, :, t] = (rg.U[:, :, below, t] * (vcoord[above] - Plev)
                            + rg.U[:, :, above, t] * (Plev - vcoord[below]))\
                / (vcoord[above] - vcoord[below])
            Vii[:, :, t] = (rg.V[:, :, below, t] * (vcoord[above] - Plev)
                            + rg.V[:, :, above, t] * (Plev - vcoord[below]))\
                / (vcoord[above] - vcoord[below])

    # Averaging in time
    if tsp > 1:
        zlevt = np.nanmean(zlev, axis=2)
        del zlev
        if wind_vectors == True:
            Uiii = np.nanmean(Uii, axis=2)
            Viii = np.nanmean(Vii, axis=2)
            del Uii, Vii
    else:
        zlevt = zlev[:, :, 0]
        del zlev
        if wind_vectors == True:
            Uiii = Uii[:, :, 0]
            Viii = Vii[:, :, 0]
            del Uii, Vii

    # smoothing
    #zlevt = ndimage.gaussian_filter(zlevt, sigma=2, order=0)

    #################
    # Create Figure #
    #################

    lonp = rg.Longitude[:]
    latp = rg.Latitude[:]

    if len(clevs) == 1:
        clevels = np.int(clevs[0])
    elif len(clevs) == 3:
        clevels = np.linspace(np.int(clevs[0]), np.int(clevs[1]), np.int(clevs[2]))
    else:
        raise IOError("clevs not valid!")
    # clevels = np.linspace(900,1470,58)
    if isinstance(axis, axes.SubplotBase):
        ax = axis
        fig = plt.gcf()
    elif (isinstance(axis, tuple)
          and isinstance(axis[0], plt.Figure)
          and isinstance(axis[1], axes.SubplotBase)):
        fig = axis[0]
        ax = axis[1]
    elif axis is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    else:
        raise IOError("'axis = {}' but {} is neither an axes.SubplotBase instance nor a (axes.SubplotBase, plt.Figure) instance".format(axis, axis))

    fig.set_tight_layout(True)

    if z['llswap']:
        C = ax.contourf(latp, lonp, zlevt.T, clevels, cmap=z['cmap'])
    else:
        C = ax.contourf(lonp, latp, zlevt, clevels, cmap=z['cmap'])

    for cc in C.collections:
        cc.set_edgecolor("face")  # fixes a stupid bug in matplotlib 2.0

    if z['llswap']:
        ax.set_xlabel('Latitude (deg)')
        ax.set_ylabel('Longitude (deg)')
    else:
        ax.set_ylabel('Latitude (deg)')
        ax.set_xlabel('Longitude (deg)')

    ax.set_title(title)

    if wind_vectors == True:
        d_z = np.shape(Uiii)
        spacing = np.int(np.shape(Uiii)[0] / 10)
        U = Uiii[::spacing, ::spacing].ravel()
        V = Viii[::spacing, ::spacing].ravel()
        lonq = loni[::spacing, ::spacing].ravel()
        latq = lati[::spacing, ::spacing].ravel()
        del Uiii, Viii
        if z['llswap']:
            q = ax.quiver(latq, lonq, V, U, color='0.5')
        else:
            q = ax.quiver(lonq, latq, U, V, color='0.5')
        ax.quiverkey(q, X=0.85, Y=-0.1, U=np.max(np.sqrt(U**2 + V**2)), label='%#.2f m/s' % np.max(np.sqrt(U**2 + V**2)), labelpos='E')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    kwargs = {'format': '%#3.4e'}
    clb = fig.colorbar(C, cax=cax, orientation='vertical', **kwargs)
    clb.set_label(z['label'])


    pfile = None
    if save == True:
        output_path = pathlib.Path(input.resultsf) / 'figures'
        if not output_path.exists():
            output_path.mkdir()

        pfile = output_path / (fname.replace(".", "+") + '.pdf')
        plt.savefig(pfile)

        plt.close()

        pfile = str(pfile)

    if z['mt'] == True:
        dname = fname + '.dat'
        if z['name'] == 'temperature-uv':
            zname = 'temperature'
        else:
            zname = z['name']
        maketable(lonp, latp, zlevt.T, 'Longitude(d)', 'Latitude(d)', zname, input.resultsf, dname)
    return pfile


def CurlF(fr, flat, flon, lat_range, lon_range, Altitude, A):
    curlFz = np.zeros_like(fr)
    curlFlat = np.zeros_like(fr)
    curlFlon = np.zeros_like(fr)

    curlFz = -1.0 * (np.gradient(flon * np.cos(lat_range[:, None, None]), lat_range, axis=0)
                     - np.gradient(flat, lon_range, axis=1))
    # curlFz = -1.0*(dFunc_dLat(flon*np.cos(lati),res_deg)-dFunc_dLon(flat,res_deg))
    curlFz /= (np.cos(lat_range[:, None, None]) * (A + Altitude[None, None, :]))

    curlFlat = -1.0 * (np.gradient(fr, lon_range, axis=1)
                       - np.gradient((A + Altitude[None, None, :]) * flon, Altitude, axis=2))
    #curlFlat = -1.0*(dFunc_dLon(fr,res_deg)/np.cos(lati)-dFunc_dZ((A+alti)*flon,alti,nv,dz))
    curlFlat /= (A + Altitude[None, None, :])

    curlFlon = -1.0 * (np.gradient((A + Altitude[None, None, :]) * flat, Altitude, axis=2)
                       - np.gradient(fr, lat_range, axis=0))
    #curlFlon = -1.0*(dFunc_dZ((A+alti)*flat,alti,nv,dz)-dFunc_dLat(fr,res_deg))
    curlFlon /= (A + Altitude[None, None, :])

    return curlFz, curlFlat, curlFlon


def streamf_moc_plot(input, grid, output, rg, sigmaref, save=True, axis=False, wind_vectors=False, mt=False, plog=True, clevs=[40]):
    # special plotting function for the mass streamfunction

    # Set the reference pressure
    Pref = input.P_Ref * sigmaref
    d_sig = np.size(sigmaref)
    tsp = output.nts - output.ntsi + 1

    if tsp > 1:
        Vavgl = np.nanmean(rg.V, axis=1)
        Vavglt = np.nanmean(Vavgl, axis=2)
        if wind_vectors == True:
            Wavgl = np.nanmean(rg.W, axis=1)
            Wavglt = np.nanmean(Wavgl, axis=2)
    else:
        Vavglt = np.nanmean(rg.V[:, :, :, 0], axis=1)
        if wind_vectors == True:
            Wavglt = np.nanmean(rg.V[:, :, :, 0], axis=1)

    sf = np.zeros_like(Vavglt)
    arg = 2 * np.pi * input.A * np.cos(rg.Latitude[:, None] * np.pi / 180) / input.Gravit * Vavglt
    for ilat in np.arange(np.shape(Vavglt)[0]):
        for ilev in np.arange(np.shape(Vavglt)[1]):
            if ilev == 0:
                sf[ilat, -1] = arg[ilat, -1] * rg.Pressure[-1]
            else:
                sf[ilat, -(ilev + 1)] = np.trapz(arg[ilat, -1:-(ilev + 2):-1], x=rg.Pressure[-1:-(ilev + 2):-1][:])

    # need to set desired pressure range (major PITA!)
    # prange = np.where(np.logical_and(rg.Pressure>=np.min(Pref),rg.Pressure<=np.max(Pref)))
    prange = np.where(rg.Pressure[:] >= np.min(Pref))

    # Contour plot
    if len(clevs) == 1:
        clevels = np.int(clevs[0])
    elif len(clevs) == 3:
        clevels = np.linspace(np.int(clevs[0]), np.int(clevs[1]), np.int(clevs[2]))
    else:
        raise IOError("clevs not valid!")

    if isinstance(axis, axes.SubplotBase):
        C = axis.contourf(rg.Latitude[:], rg.Pressure[prange[0]] / 1e5, sf[:, prange[0]].T, clevels, cmap='viridis')
        ax = axis
    elif axis == False:
        plt.figure(figsize=(5, 4))
        C = plt.contourf(rg.Latitude[:], rg.Pressure[prange[0]] / 1e5, sf[:, prange[0]].T, clevels, cmap='viridis')

        ax = plt.gca()
    else:
        raise IOError("'axis = {}' but {} is not an axes.SubplotBase instance".format(axis, axis))

    for cc in C.collections:
        cc.set_edgecolor("face")  # fixes a stupid bug in matplotlib 2.0

    if wind_vectors == True:
        vspacing = np.int(np.shape(rg.Latitude)[0] / 10)
        wspacing = np.int(np.shape(rg.Pressure)[0] / 10)
        Vlt = Vavglt[:, prange[0]]
        Wlt = Wavglt[:, prange[0]]
        Vq = Vlt[::vspacing, ::wspacing].ravel()
        Wq = Wlt[::vspacing, ::wspacing].ravel()
        #preq = rg.Pressure[:,0][::spacing,::spacing].ravel()
        #latq = lati[::spacing,::spacing].ravel()
        latq, preq = np.meshgrid(rg.Latitude[::vspacing], rg.Pressure[::wspacing][prange[0]])
        del Vlt, Wlt
        plt.quiver(latq.ravel(), preq.ravel() / 1e5, Vq, Wq, color='0.5')

    ax.invert_yaxis()
    c2 = ax.contour(rg.Latitude[:], rg.Pressure[:] / 1e5, sf.T, levels=[0.0], colors='w', linewidths=1)
    clb = plt.colorbar(C)
    clb.set_label(r'Eulerian streamfunction (kg s$^{-1}$)')
    if plog == True:
        ax.set_yscale("log")
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Pressure (bar)')
    # ax.plot(rg.Latitude[:,0],np.zeros_like(rg.Latitude[:,0])+np.max(output.Pressure[:,grid.nv-1,:])/1e5,'r--')
    # if np.min(rg.Pressure[prange[0],0]) < np.max(output.Pressure[:,grid.nv-1,:]):
    ax.set_ylim(np.max(rg.Pressure[prange[0]]) / 1e5, np.min(rg.Pressure[prange[0]]) / 1e5)

    # else:
    #     ax.set_ylim(np.max(rg.Pressure[prange[0],0])/1e5,np.max(output.Pressure[:,grid.nv-1,:])/1e5)
    ax.set_title('Time = %#.1f-%#.1f days, Lon = (0,360)' % (output.time[0], output.time[-1]), fontsize=10)
    if not os.path.exists(input.resultsf + '/figures'):
        os.mkdir(input.resultsf + '/figures')
    plt.tight_layout()
    pfile = None
    if save == True:
        pfile = input.resultsf + '/figures/streamf_ver_i%d_l%d.pdf' % (output.ntsi, output.nts)
        plt.savefig(pfile)
        plt.close()

        pfile = str(pfile)

    if mt == True:
        fname = 'streamf_ver_i%d_l%d.dat' % (output.ntsi, output.nts)
        maketable(latp * 180 / np.pi, rg.Pressure[prange[0]] / 1e5, Zonallt[:, prange[0]],
                  'Latitude(d)', 'Pressure(bar)', 'streamfunc', input.resultsf, fname)
    return pfile


def profile(input, grid, output, z, stride=50, axis=None, save=True):
    # Pref = input.P_Ref*sigmaref
    # d_sig = np.size(sigmaref)

    # prepare pressure array
    output.load_reshape(grid, ['Pressure'])

    # get figure and axis
    if isinstance(axis, axes.SubplotBase):
        ax = axis
        fig = plt.gcf()
    elif (isinstance(axis, tuple)
          and isinstance(axis[0], plt.Figure)
          and isinstance(axis[1], axes.SubplotBase)):
        fig = axis[0]
        ax = axis[1]
    elif axis is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    else:
        raise IOError("'axis = {}' but {} is neither an axes.SubplotBase instance nor a (axes.SubplotBase, plt.Figure) instance".format(axis, axis))

    fig.set_tight_layout(True)

    tsp = output.nts - output.ntsi + 1

    col_lon = []
    col_lat = []
    col_lor = []
    for column in np.arange(0, grid.point_num, stride):
        if tsp > 1:
            P = np.mean(output.Pressure[column, :, :], axis=1)
            x = np.mean(z['value'][column, :, :], axis=1)
        else:
            P = output.Pressure[column, :, 0]
            x = z['value'][column, :, 0]

        #color = hsv_to_rgb([column / grid.point_num, 1.0, 1.0])
        lon = grid.lon[column]
        lat = grid.lat[column]
        lon_norm = lon / (2.0 * math.pi)
        lat_norm = lat / (math.pi / 2.0)
        if lat > 0.0:
            color = hsv_to_rgb([lon_norm, 1.0 - 0.7 * lat_norm, 1.0])
        else:
            color = hsv_to_rgb([lon_norm, 1.0, 1.0 + 0.7 * lat_norm])

        col_lon.append(lon)
        col_lat.append(lat)
        col_lor.append(color)
        ax.semilogy(x, P / 1e5, 'k-', alpha=0.5, lw=1.0,
                    path_effects=[pe.Stroke(linewidth=1.5, foreground=color), pe.Normal()])

        rp, = ax.plot(x[np.int(np.floor(grid.nv / 2))], P[np.int(np.floor(grid.nv / 2))] / 100000, 'k+', ms=5, alpha=0.5)
        gp, = ax.plot(x[np.int(np.floor(grid.nv * 0.75))], P[np.int(np.floor(grid.nv * 0.75))] / 100000, 'k*', ms=5, alpha=0.5)

    # add an insert showing the position of
    inset_pos = [0.8, 0.8, 0.18, 0.18]
    ax_inset = ax.inset_axes(inset_pos)
    ax_inset.scatter(col_lon, col_lat, c=col_lor, s=1.0)
    ax_inset.tick_params(axis='both',
                         which='both',
                         top=False,
                         bottom=False,
                         left=False,
                         right=False,
                         labelright=False,
                         labelleft=False,
                         labeltop=False,
                         labelbottom=False)

    # plt.plot(Tad,P/100,'r--')
    ax.invert_yaxis()
    ax.set_ylabel('Pressure (bar)')
    ax.set_xlabel(z['label'])
    ax.legend([rp, gp], ['z=0.5*ztop', 'z=0.75*ztop'], loc="lower right", fontsize='xx-small')
    ax.set_title('Time = %#.3f - %#.3f days' % (output.time[0], output.time[-1]))
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.tick_params(axis='x', labelrotation=45 )

    ax.grid(True)
    pfile = None
    if save == True:
        output_path = pathlib.Path(input.resultsf) / 'figures'
        if not output_path.exists():
            output_path.mkdir()

        fname = f"{z['name']}Pprofile_i{output.ntsi}_l{output.nts}"
        pfile = output_path / (fname.replace(".", "+") + '.pdf')
        plt.savefig(pfile)

        plt.close()
        pfile = str(pfile)

    return pfile


def CalcE_M_AM(input, grid, output, split):
    temperature = output.Pressure / (input.Rd * output.Rho)
    # dz = grid.Altitude[1]-grid.Altitude[0]
    # Vol = grid.areasT*dz

    # testing @@
    tsp = output.nts - output.ntsi + 1

    Atot = input.A**2
    solid_ang = grid.areasT / Atot

    rint = ((input.A + grid.Altitudeh[1:])**3 - (input.A + grid.Altitudeh[:-1])**3) / 3.0
    Vol0 = solid_ang[:, None] * rint[None, :]

    Eint = output.Rho * (input.Cp - input.Rd) * temperature
    Eg = output.Rho * input.Gravit * grid.Altitude[None, :, None]
    W = 0.5 * (output.Wh[:, 1:, :] + output.Wh[:, :-1, :])
    Wx = W * np.cos(grid.lat[:, None, None]) * np.cos(grid.lon[:, None, None])
    Wy = W * np.cos(grid.lat[:, None, None]) * np.sin(grid.lon[:, None, None])
    Wz = W * np.sin(grid.lat[:, None, None])

    Mtotx = output.Mh[0] + Wx
    Mtoty = output.Mh[1] + Wy
    Mtotz = output.Mh[2] + Wz
    Mtot2 = (Mtotx)**2 + (Mtoty)**2 + (Mtotz)**2

    Ek = 0.5 * Mtot2 / output.Rho

    output.Etotal = (Eint + Eg + Ek) * Vol0[:, :, None]

    output.Mass = output.Rho * Vol0[:, :, None]

    r = input.A + grid.Altitude

    rx = r[None, :, None] * np.cos(grid.lat[:, None, None]) * np.cos(grid.lon[:, None, None])
    ry = r[None, :, None] * np.cos(grid.lat[:, None, None]) * np.sin(grid.lon[:, None, None])
    rz = r[None, :, None] * np.sin(grid.lat[:, None, None])

    output.AngMomx = (ry * Mtotz - rz * Mtoty - output.Rho * input.Omega[0] *
                      rz * r[None, :, None] * np.cos(grid.lat[:, None, None]) *
                      np.cos(grid.lon[:, None, None])) * Vol0[:, :, None]
    output.AngMomy = (-rx * Mtotz + rz * Mtotx - output.Rho * input.Omega[0] *
                      rz * r[None, :, None] * np.cos(grid.lat[:, None, None]) *
                      np.sin(grid.lon[:, None, None])) * Vol0[:, :, None]
    output.AngMomz = (rx * Mtoty - ry * Mtotx + output.Rho * input.Omega[0] *
                      r[None, :, None]**2 * np.cos(grid.lat[:, None, None]) *
                      np.cos(grid.lat[:, None, None])) * Vol0[:, :, None]

    if split == False:
        output.GlobalE = np.sum(np.sum(output.Etotal, 0), 0)
        output.GlobalMass = np.sum(np.sum(output.Mass, 0), 0)
        output.GlobalAMx = np.sum(np.sum(output.AngMomx, 0), 0)
        output.GlobalAMy = np.sum(np.sum(output.AngMomy, 0), 0)
        output.GlobalAMz = np.sum(np.sum(output.AngMomz, 0), 0)
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
            output.WeatherE[ii] = np.sum(output.Etotal[output.Pressure[:, :, ii] <= split][:, ii])
            output.DeepE[ii] = np.sum(output.Etotal[output.Pressure[:, :, ii] > split][:, ii])
            output.WeatherMass[ii] = np.sum(output.Mass[output.Pressure[:, :, ii] <= split][:, ii])
            output.DeepMass[ii] = np.sum(output.Mass[output.Pressure[:, :, ii] > split][:, ii])
            output.WeatherAMx[ii] = np.sum(output.AngMomx[output.Pressure[:, :, ii] <= split][:, ii])
            output.DeepAMx[ii] = np.sum(output.AngMomx[output.Pressure[:, :, ii] > split][:, ii])
            output.WeatherAMy[ii] = np.sum(output.AngMomy[output.Pressure[:, :, ii] <= split][:, ii])
            output.DeepAMy[ii] = np.sum(output.AngMomy[output.Pressure[:, :, ii] > split][:, ii])
            output.WeatherAMz[ii] = np.sum(output.AngMomz[output.Pressure[:, :, ii] <= split][:, ii])
            output.DeepAMz[ii] = np.sum(output.AngMomz[output.Pressure[:, :, ii] > split][:, ii])


def CalcEntropy(input, grid, output, split):
    temperature = output.Pressure / (input.Rd * output.Rho)
    # dz = grid.Altitude[1]-grid.Altitude[0]
    # Vol = grid.areasT*dz

    Atot = input.A**2
    solid_ang = grid.areasT / Atot

    rint = ((input.A + grid.Altitudeh[1:])**3 - (input.A + grid.Altitudeh[:-1])**3) / 3.0
    Vol0 = solid_ang[:, None] * rint[None, :]

    kappa = input.Rd / input.Cp

    potT = temperature * (input.P_Ref / output.Pressure)**kappa
    S = input.Cp * np.log(potT)
    output.Entropy = S * Vol0[:, :, None]
    if split == False:
        output.GlobalEnt = np.sum(np.sum(output.Entropy, 0), 0)
    else:
        output.WeatherEnt = np.zeros(np.shape(output.Pressure)[2])
        output.DeepEnt = np.zeros(np.shape(output.Pressure)[2])
        for ii in np.arange(np.shape(output.Pressure)[2]):
            output.WeatherEnt[ii] = np.sum(output.Entropy[output.Pressure[:, :, ii] <= split][:, ii])
            output.DeepEnt[ii] = np.sum(output.Entropy[output.Pressure[:, :, ii] > split][:, ii])


def conservation(input, grid, output, split):
    # plot quantities that are interesting for conservation
    if split == False:
        if (output.ConvData == False).any():
            print('Calculating energy, mass, angular momentum...')
            CalcE_M_AM(input, grid, output, split)
        if (output.EntData == False).any():
            print('Calculating entropy...')
            CalcEntropy(input, grid, output, split)
        plots = ['global']

    else:
        CalcE_M_AM(input, grid, output, split)
        CalcEntropy(input, grid, output, split)
        plots = ['weather', 'deep']

    for ii in np.arange(len(plots)):
        fig = plt.figure(figsize=(12, 8))
        fig.suptitle('Time = %#.3f - %#.3f days, split = %e bar' % (output.time[0], output.time[-1], split / 1e5))
        fig.subplots_adjust(wspace=0.25, left=0.07, right=0.98, top=0.94, bottom=0.07)
        plt.subplot(2, 3, 1)
        if split == False:
            plt.plot(output.time, output.GlobalE, 'ko', linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time, output.WeatherE, 'bo', linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time, output.DeepE, 'ro', linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel('Total energy of atmosphere (J)')

        plt.subplot(2, 3, 2)
        if split == False:
            plt.plot(output.time, output.GlobalMass, 'ko', linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time, output.WeatherMass, 'bo', linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time, output.DeepMass, 'ro', linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel('Total mass of atmosphere (kg)')

        plt.subplot(2, 3, 3)
        if split == False:
            plt.plot(output.time, output.GlobalAMz, 'ko', linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time, output.WeatherAMz, 'bo', linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time, output.DeepAMz, 'ro', linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel(r'Z angular momentum of atmosphere (kg m$^2$ s$^{-1}$)')

        plt.subplot(2, 3, 4)
        if split == False:
            plt.plot(output.time, output.GlobalEnt, 'ko', linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time, output.WeatherEnt, 'bo', linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time, output.DeepEnt, 'ro', linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel('Total entropy of atmosphere (J K$^{-1}$)')

        plt.subplot(2, 3, 5)
        if split == False:
            plt.plot(output.time, output.GlobalAMx, 'ko', linestyle='--')
            plt.plot(output.time, output.GlobalAMy, 'o', color='0.5', linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time, output.WeatherAMx, 'bo', linestyle='--')
                plt.plot(output.time, output.WeatherAMy, 'co', linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time, output.DeepAMx, 'ro', linestyle='--')
                plt.plot(output.time, output.DeepAMy, 'o', color='orange', linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel(r'X, Y angular momentum of atmosphere (kg m$^2$ s$^{-1}$)')

        plt.subplot(2, 3, 6)
        if split == False:
            plt.plot(output.time, np.sqrt(output.GlobalAMx**2 + output.GlobalAMy**2), 'ko', linestyle='--')
        else:
            if plots[ii] == 'weather':
                plt.plot(output.time, np.sqrt(output.WeatherAMx**2 + output.WeatherAMy**2), 'bo', linestyle='--')
            elif plots[ii] == 'deep':
                plt.plot(output.time, np.sqrt(output.DeepAMx**2 + output.DeepAMy**2), 'ro', linestyle='--')
        plt.xlabel('Time (days)')
        plt.ylabel(r'Horizontal angular momentum of atmosphere (kg m$^2$ s$^{-1}$)')

        if not os.path.exists(input.resultsf + '/figures'):
            os.mkdir(input.resultsf + '/figures')
        plt.savefig(input.resultsf + '/figures/conservation_s%e_i%d_l%d_%s.pdf' % (split / 100, output.ntsi, output.nts, plots[ii]))
        plt.close()


def SRindex(input, grid, output):
    tsp = output.nts - output.ntsi + 1
    # total surface area
    Atot = input.A**2
    solid_ang = grid.areasT / Atot

    U = (-output.Mh[0] * np.sin(grid.lon[:, None, None]) + output.Mh[1] * np.cos(grid.lon[:, None, None])) / output.Rho

    # should I adjust radius for height?? yes!!
    mom = (input.A + grid.Altitude[None, :, None]) * np.cos(grid.lat[:, None, None]) *\
          ((input.A + grid.Altitude[None, :, None]) * input.Omega * np.cos(grid.lat[:, None, None]) + U)
    mom0 = (input.A + grid.Altitude[None, :, None]) * np.cos(grid.lat[:, None, None]) *\
        ((input.A + grid.Altitude[None, :, None]) * input.Omega * np.cos(grid.lat[:, None, None]) + np.zeros(np.shape(U[:, :, :])))

    rint = ((input.A + grid.Altitudeh[1:])**3 - (input.A + grid.Altitudeh[:-1])**3) / 3.0
    mom_grid = mom * rint[None, :, None] * output.Rho
    mom_grid0 = mom0 * rint[None, :, None] * output.Rho

    Mtot = np.sum(np.sum(mom_grid * solid_ang[:, None, None], axis=0), axis=0)
    Mtot0 = np.sum(np.sum(mom_grid0 * solid_ang[:, None, None], axis=0), axis=0)

    S = Mtot / Mtot0 - 1

    plt.semilogy(output.time, S)
    plt.plot(output.time, -1 * S)
    plt.xlabel('Time (days)')
    plt.ylabel('Super-rotation index')
    plt.tight_layout()

    if not os.path.exists(input.resultsf + '/figures'):
        os.mkdir(input.resultsf + '/figures')
    plt.savefig(input.resultsf + '/figures/SRindex_i%d_l%d.pdf' % (output.ntsi, output.nts))
    plt.close()

def phase_curve(input, grid, output, axis=None, save=True):
    #plot phase curve, averaged over output times

    # get figure and axis
    if isinstance(axis, axes.SubplotBase):
        ax = axis
        fig = plt.gcf()
    elif (isinstance(axis, tuple)
          and isinstance(axis[0], plt.Figure)
          and isinstance(axis[1], axes.SubplotBase)):
        fig = axis[0]
        ax = axis[1]
    elif axis is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    else:
        raise IOError("'axis = {}' but {} is neither an axes.SubplotBase instance nor a (axes.SubplotBase, plt.Figure) instance".format(axis, axis))

    angle = np.linspace(0,2*np.pi,100)

    pc = np.zeros_like(angle)
    for j in np.arange(len(output.time)):
        olr = output.flw_up[:,-1,j]
        for i in np.arange(len(angle)):
            mu = np.cos(grid.lat)*np.cos(grid.lon-angle[i])
            mu[mu<0] = 0
            los_flux = mu*olr
            pc_tmp = np.sum(los_flux*grid.areasT)/np.sum(grid.areasT)
            pc[i] = pc_tmp/(j+1)+pc[i]*(j)/(j+1)

    xx = 360-((np.pi+angle)%(2*np.pi))*180/np.pi
    ax.plot(np.sort(xx),pc[np.argsort(xx)])

    ax.set_title('Time = %#.3f - %#.3f days' % (output.time[0], output.time[-1]))

    ax.axvline(180,ls='--',c = 'k')
    ax.set_xlabel('Angle relative to transit (deg)')
    ax.set_ylabel(r'Flux (W m$^{-2}$)')
    fig.tight_layout()
    if save == True:
        output_path = pathlib.Path(input.resultsf) / 'figures'
        if not output_path.exists():
            output_path.mkdir()

        fname = f"phase_curve_i{output.ntsi}_l{output.nts}.pdf"
        pfile = output_path / (fname.replace(".", "+") + '.pdf')
        plt.savefig(pfile)

        plt.close()

def Get_Prange(input, grid, rg, args, xtype='lat', use_p=True):
    # Sigma (normalized pressure) values for the plotting
    if not isinstance(args.slice, list):
        raise IOError("'slice' argument must be a list")

    if use_p:
        if (args.vertical_top[0] == 'default'):
            args.vertical_top[0] = np.min(rg.Pressure) / 100

        if np.max(input.P_Ref) / np.float(args.vertical_top[0]) > 1000:
            sigmaref = np.logspace(np.log10(np.max(rg.Pressure)), np.log10(np.float(args.vertical_top[0]) * 100), 20) / input.P_Ref
        else:
            sigmaref = np.linspace(np.max(rg.Pressure), np.float(args.vertical_top[0]) * 100, 20) / input.P_Ref

    else:
        if (args.vertical_top[0] == 'default'):
            sigmaref = grid.Altitude
        else:
            sigmaref = grid.Altitude[np.where(grid.Altitude <= np.float(args.vertical_top[0]) * input.Top_altitude)[0]]
    return sigmaref


def RTbalance(input, grid, output):
    # not finished!
    asr = fsw_dn[:, input.nv - 1, :] * grid.areasT[:, None]
    olr = flw_up[:, input.nv - 1, :] * grid.areasT[:, None]

def spectrum(input, grid, output, z, stride=20, axis=None, save=True):
    output.load_reshape(grid, ['spectrum', 'incoming_spectrum'])
    # get figure and axis
    if isinstance(axis, axes.SubplotBase):
        ax = axis
        fig = plt.gcf()
    elif (isinstance(axis, tuple)
          and isinstance(axis[0], plt.Figure)
          and isinstance(axis[1], axes.SubplotBase)):
        fig = axis[0]
        ax = axis[1]
    elif axis is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    else:
        raise IOError("'axis = {}' but {} is neither an axes.SubplotBase instance nor a (axes.SubplotBase, plt.Figure) instance".format(axis, axis))

    fig.set_tight_layout(True)

    tsp = output.nts - output.ntsi + 1


    col_lon = []
    col_lat = []
    col_lor = []

    for column in np.arange(0, grid.point_num, stride):
        lamda = output.wavelength[:]*1e6
        if tsp > 1:
            spectrum= np.mean(output.spectrum[column, :, :], axis=1)
        else:
            spectrum= output.spectrum[column, :, 0]

        #color = hsv_to_rgb([column / grid.point_num, 1.0, 1.0])
        lon = grid.lon[column]
        lat = grid.lat[column]
        lon_norm = lon / (2.0 * math.pi)
        lat_norm = lat / (math.pi / 2.0)
        if lat > 0.0:
            color = hsv_to_rgb([lon_norm, 1.0 - 0.7 * lat_norm, 1.0])
        else:
            color = hsv_to_rgb([lon_norm, 1.0, 1.0 + 0.7 * lat_norm])

        col_lon.append(lon)
        col_lat.append(lat)
        col_lor.append(color)
        ax.plot(lamda, spectrum, 'k-', alpha=0.5, lw=1.0,
                    path_effects=[pe.Stroke(linewidth=1.5, foreground=color), pe.Normal()])

    lamda = output.wavelength[:]*1e6
    if tsp > 1:
        incoming_spectrum= np.mean(output.incoming_spectrum[:, :], axis=1)
    else:
        incoming_spectrum= output.incoming_spectrum[:, 0]
    R_SUN   = 695700000.0
    AU       = 149597870700.0
    flx_cst = pow((input.radius_star*R_SUN)/(input.planet_star_dist*AU), 2.0)*math.pi
    ax_twin = ax.twinx()
    ax_twin.plot(lamda, incoming_spectrum*flx_cst, 'k-', alpha=0.5, lw=1.0,
            path_effects=[pe.Stroke(linewidth=1.5, foreground='black'), pe.Normal()])
    ax_twin.set_ylabel('Incoming stellar flux [$W m^2$]')
    # add an insert showing the position of columns
    inset_pos = [0.8, 0.1, 0.18, 0.18]
    ax_inset = ax.inset_axes(inset_pos)
    ax_inset.scatter(col_lon, col_lat, c=col_lor, s=1.0)
    ax_inset.tick_params(axis='both',
                         which='both',
                         top=False,
                         bottom=False,
                         left=False,
                         right=False,
                         labelright=False,
                         labelleft=False,
                         labeltop=False,
                         labelbottom=False)

    # plt.plot(Tad,P/100,'r--')
    #ax.invert_yaxis()
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.set_ylabel('Flux [$W m^2$]')
    ax.set_xlabel('Wavelength [um]')
    ax.set_title('Time = %#.3f - %#.3f days' % (output.time[0], output.time[-1]))
    ax.grid(True)

    pfile = None
    if save == True:
        output_path = pathlib.Path(input.resultsf) / 'figures'
        if not output_path.exists():
            output_path.mkdir()

        fname = f"spectrum_i{output.ntsi}_l{output.nts}"
        pfile = output_path / (fname.replace(".", "+") + '.pdf')
        plt.savefig(pfile)

        plt.close()
        pfile = str(pfile)

    return pfile
