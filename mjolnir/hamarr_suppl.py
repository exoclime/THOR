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
from hamarr import GetOutput, curlF

#This is mostly old code that might still be useful in some situations
#I don't plan on maintaining this so you're on your own if you want to use it


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


def calc_RV_PV(grid,output,input,lons,lats,sigma,t_ind,fileh5,comp=4,pressure_vert=True,type='gd',lmax_set='grid'):
    #Calculates relative and potential vorticity on height levels, then interpolates to pressure.
    #It is easiest to calculate these on a lat-lon-altitude grid since conversion
    #to pressure requires transformation of vertical velocity and computing 3D curls
    #at the grid level would require heavy interpolation and invocation of stokes thm
    ang_res = 4.0/2**(input.glevel-4)
    # first interpolate to near thor grid resolution
    if type == 'sh' or type == 'SH':
        if lmax_set == 'grid':
            lmax = np.int(np.sqrt(grid.point_num)*0.5)
        elif isinstance(lmax_set,int):
            lmax = lmax_set
        else:
            raise ValueError("Invalid value of lmax in regrid")
        n = 2*lmax+2
        lat_range = np.arange(90,-90,-180/n)
        lon_range = np.arange(0,360,180/n)
    else:
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
        if type == 'gd' or type == 'GD':
            pt_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,pt[:,lev],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
            rho_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,output.Rho[:,lev,t_ind],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
            p_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,output.Pressure[:,lev,t_ind],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
            U_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,U[:,lev],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
            V_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,V[:,lev],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
            W_llz[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,W[:,lev],(loni[:,:,lev],lati[:,:,lev]),method='nearest')
        elif type == 'sh' or type == 'SH':
            pt_coeff, pt_chi2 = chairs.expand.SHExpandLSQ(pt[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
            pt_llz[:,:,lev] = np.real(chairs.expand.MakeGridDH(pt_coeff,sampling=2))
            rho_coeff, rho_chi2 = chairs.expand.SHExpandLSQ(output.Rho[:,lev,t_ind],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
            rho_llz[:,:,lev] = np.real(chairs.expand.MakeGridDH(rho_coeff,sampling=2))
            p_coeff, p_chi2 = chairs.expand.SHExpandLSQ(output.Pressure[:,lev,t_ind],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
            p_llz[:,:,lev] = np.real(chairs.expand.MakeGridDH(p_coeff,sampling=2))
            U_coeff, U_chi2 = chairs.expand.SHExpandLSQ(U[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
            U_llz[:,:,lev] = np.real(chairs.expand.MakeGridDH(U_coeff,sampling=2))
            V_coeff, V_chi2 = chairs.expand.SHExpandLSQ(V[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
            V_llz[:,:,lev] = np.real(chairs.expand.MakeGridDH(V_coeff,sampling=2))
            W_coeff, W_chi2 = chairs.expand.SHExpandLSQ(W[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
            W_llz[:,:,lev] = np.real(chairs.expand.MakeGridDH(W_coeff,sampling=2))
        elif type == 'spl' or type == 'SPL':
            theta = np.pi/2 - grid.lat
            f = interp.SmoothSphereBivariateSpline(theta,grid.lon,pt[:,lev],s=3.5)
            pt_llz[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
            f = interp.SmoothSphereBivariateSpline(theta,grid.lon,output.Rho[:,lev,t_ind],s=3.5)
            rho_llz[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
            f = interp.SmoothSphereBivariateSpline(theta,grid.lon,output.Pressure[:,lev,t_ind],s=3.5)
            p_llz[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
            f = interp.SmoothSphereBivariateSpline(theta,grid.lon,U[:,lev],s=3.5)
            U_llz[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
            f = interp.SmoothSphereBivariateSpline(theta,grid.lon,V[:,lev],s=3.5)
            V_llz[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
            f = interp.SmoothSphereBivariateSpline(theta,grid.lon,W[:,lev],s=3.5)
            W_llz[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)

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

def regrid_old(resultsf,simID,ntsi,nts,nlev=40,pgrid_ref='auto',overwrite=False,comp=4,
            pressure_vert=True,type='gd',vertical_top='default',rotation=False,theta_z=0,theta_y=0,
            lmax_set='grid',mask_surf=True):
    # runs over files and converts ico-height grid to lat-lon-pr grid
    outall = GetOutput(resultsf,simID,ntsi,ntsi,rotation=rotation,theta_z=theta_z,theta_y=theta_y)
    input = outall.input
    grid = outall.grid
    output = outall.output
    res_deg = 4.0/2**(input.glevel-4)

    print('Regrid data in folder '+resultsf+'...\n')
    if pressure_vert == True:
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
        # sigmaref = np.mean(output.Pressure,axis=0)/input.P_Ref
        # if pgrid_ref == 'mean':
        #     Pref = input.P_Ref*np.mean(sigmaref,axis=1)
        # elif pgrid_ref == 'first':
        #     Pref = input.P_Ref*sigmaref[:,0]
        # elif pgrid_ref == 'last':
        #     Pref = input.P_Ref*sigmaref[:,-1]
        # else:
        #     raise ValueError("Invalid value for pgrid_ref (valid = mean, first, or last)")
    else:
        print('Vertical coordinate = height')
        nlev = len(grid.Altitude)
        Pref = np.zeros(nlev)

    if type == 'sh' or type == 'SH':
        if lmax_set == 'grid':
            lmax = np.int(np.sqrt(grid.point_num)*0.5)
        elif isinstance(lmax_set,int):
            lmax = lmax_set
        else:
            raise ValueError("Invalid value of lmax in regrid")
        n = 2*lmax+2
        lat_range = np.arange(90,-90,-180/n)
        lon_range = np.arange(0,360,180/n)
    else:
        lat_range_tmp = np.arange(-90,90+res_deg,res_deg)
        # recenter so that there are an even number of latitude points
        lat_range = (lat_range_tmp[:-1]+lat_range_tmp[1:])/2
        lon_range = np.arange(0,360,res_deg)
    loni, lati = np.meshgrid(lon_range,lat_range)
    d_lon = np.shape(loni)
    tsp = nts-ntsi+1

    RT = 0
    if hasattr(input,"core_benchmark"):
        if input.core_benchmark[0] == 0: # need to switch to 'radiative_tranfer' flag
            RT = 1
            surf = 0
            if input.surface:
                surf = 1
                extrap_low = (grid.Altitudeh[0]-grid.Altitude[1])/(grid.Altitude[0]-grid.Altitude[1])
                Psurf = output.Pressure[:,1,:]+extrap_low*(output.Pressure[:,0,:]-output.Pressure[:,1,:])

                # if pgrid_ref == 'mean':
                #     Pref = np.concatenate((np.array([np.mean(Psurf)]),Pref))
                # elif pgrid_ref == 'first':
                #     Pref = np.concatenate((np.array([np.mean(Psurf[:,0])]),Pref))
                # elif pgrid_ref == 'last':
                #     Pref = np.concatenate((np.array([np.mean(Psurf[:,-1])]),Pref))
        else:
            surf = 0

    d_sig = np.size(Pref)

    chem = 0
    if hasattr(input,"chemistry"):
        if input.chemistry == 1:
            chem = 1

    for t in np.arange(ntsi,nts+1):
        #check for exising h5 files
        proceed = 0
        if pressure_vert == True:
            fileh5 = resultsf+'/regrid_'+simID+'_'+np.str(t)
        else:
            fileh5 = resultsf+'/regrid_height_'+simID+'_'+np.str(t)
        if rotation == True:
            print('Applied rotation (theta_z,theta_y) = (%f,%f) to grid\n'%(theta_z*180/np.pi,theta_y*180/np.pi))
            # fileh5 += '_rot'
        fileh5 += '.h5'
        if os.path.exists(fileh5):
            if overwrite == True:
                proceed = 1
            else:
                print(fileh5+' already present! Skipping time = %d'%t)
        else:
            proceed = 1

        if proceed == 1:
            outall = GetOutput(resultsf,simID,t,t,rotation=rotation,theta_z=theta_z,theta_y=theta_y)
            output = outall.output

            Temp_icoh = output.Pressure/(input.Rd*output.Rho)
            interpx = (grid.Altitude-grid.Altitudeh[:-1])/(grid.Altitudeh[1:]-grid.Altitudeh[:-1])
            # on even height grid, interpolation is excessive, but wth?
            W_icoh = output.Wh[:,:-1,:] + (output.Wh[:,1:,:]-output.Wh[:,:-1,:])*interpx[None,:,None]
            # del_hseq = np.gradient(output.Pressure,grid.Altitude,axis=1) + output.Rho*input.Gravit
            if RT == 1:
                flw_up_icoh = output.flw_up[:,:-1,:]+(output.flw_up[:,1:,:]-output.flw_up[:,:-1,:])*interpx[None,:,None]
                flw_dn_icoh = output.flw_dn[:,:-1,:]+(output.flw_dn[:,1:,:]-output.flw_dn[:,:-1,:])*interpx[None,:,None]

            print('Regridding time = %d...'%t)
            Temp_icop = np.zeros((grid.point_num,d_sig))
            Temp_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            Rho_icop = np.zeros((grid.point_num,d_sig))
            Rho_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            Mh_icop = np.zeros((3,grid.point_num,d_sig))

            W_icop = np.zeros((grid.point_num,d_sig))

            U_llp = np.zeros((d_lon[0],d_lon[1],d_sig))
            V_llp = np.zeros((d_lon[0],d_lon[1],d_sig))
            W_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            # del_hseq_icop = np.zeros((grid.point_num,d_sig))
            # del_hseq_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

            if type == 'sh' or type == 'SH':
                Temp_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                Temp_chi2 = np.zeros(d_sig)

                Rho_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                Rho_chi2 = np.zeros(d_sig)

                U_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                U_chi2 = np.zeros(d_sig)

                V_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                V_chi2 = np.zeros(d_sig)

                W_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                W_chi2 = np.zeros(d_sig)

            if RT == 1:
                tau_sw_icop = np.zeros((grid.point_num,d_sig))
                tau_sw_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                tau_lw_icop = np.zeros((grid.point_num,d_sig))
                tau_lw_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                flw_up_icop = np.zeros((grid.point_num,d_sig))
                flw_up_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                flw_dn_icop = np.zeros((grid.point_num,d_sig))
                flw_dn_llp = np.zeros((d_lon[0],d_lon[1],d_sig))

                insol_ll = np.zeros((d_lon[0],d_lon[1]))
                if surf == 1:
                    Tsurf_ll = np.zeros((d_lon[0],d_lon[1]))

                if type == 'sh' or type == 'SH':
                    tau_sw_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                    tau_sw_chi2 = np.zeros(d_sig)

                    tau_lw_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                    tau_lw_chi2 = np.zeros(d_sig)

                    flw_up_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                    flw_up_chi2 = np.zeros(d_sig)

                    flw_dn_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                    flw_dn_chi2 = np.zeros(d_sig)

                    insol_coeff = np.zeros((1,2,lmax+1,lmax+1))
                    insol_chi2 = np.zeros(1)

                    if surf == 1:
                        Tsurf_coeff = np.zeros((1,2,lmax+1,lmax+1))
                        Tsurf_chi2 = np.zeros(1)

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

                if type == 'sh' or type == 'SH':
                    ch4_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                    ch4_chi2 = np.zeros(d_sig)

                    co_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                    co_chi2 = np.zeros(d_sig)

                    h2o_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                    h2o_chi2 = np.zeros(d_sig)

                    co2_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                    co2_chi2 = np.zeros(d_sig)

                    nh3_coeff = np.zeros((d_sig,2,lmax+1,lmax+1))
                    nh3_chi2 = np.zeros(d_sig)

            if pressure_vert == True: # use pressure as vertical coordinate
                for i in np.arange(grid.point_num):
                    #interp to pressure grid
                    sigma = output.Pressure[i,:,0]
                    Temp_icop[i,:] = interp.pchip_interpolate(sigma[::-1],Temp_icoh[i,::-1,0],Pref[::-1])[::-1]
                    Rho_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.Rho[i,::-1,0],Pref[::-1])[::-1]
                    Mh_icop[0,i,:] = interp.pchip_interpolate(sigma[::-1],output.Mh[0,i,::-1,0],Pref[::-1])[::-1]
                    Mh_icop[1,i,:] = interp.pchip_interpolate(sigma[::-1],output.Mh[1,i,::-1,0],Pref[::-1])[::-1]
                    Mh_icop[2,i,:] = interp.pchip_interpolate(sigma[::-1],output.Mh[2,i,::-1,0],Pref[::-1])[::-1]
                    W_icop[i,:] = interp.pchip_interpolate(sigma[::-1],W_icoh[i,::-1,0],Pref[::-1])[::-1]
                    # del_hseq_icop[i,:] = interp.pchip_interpolate(sigma[::-1],del_hseq[i,::-1,0],Pref[::-1])[::-1]
                    if surf and mask_surf:
                        #when isobars intersect the surface, make values below surface undefined
                        Temp_icop[i,Pref>Psurf[i,0]] = np.nan
                        Rho_icop[i,Pref>Psurf[i,0]] = np.nan
                        Mh_icop[:,i,Pref>Psurf[i,0]] = np.nan
                        W_icop[i,Pref>Psurf[i,0]] = np.nan

                    if RT == 1:
                        tau_sw_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.tau_sw[i,::-1,0],Pref[::-1])[::-1]
                        tau_lw_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.tau_lw[i,::-1,0],Pref[::-1])[::-1]
                        flw_up_icop[i,:] = interp.pchip_interpolate(sigma[::-1],flw_up_icoh[i,::-1,0],Pref[::-1])[::-1]
                        flw_dn_icop[i,:] = interp.pchip_interpolate(sigma[::-1],flw_dn_icoh[i,::-1,0],Pref[::-1])[::-1]
                        if surf and mask_surf:
                            #when isobars intersect the surface, make values below surface undefined
                            tau_sw_icop[i,Pref>Psurf[i,0]] = np.nan
                            tau_lw_icop[i,Pref>Psurf[i,0]] = np.nan
                            flw_up_icop[i,Pref>Psurf[i,0]] = np.nan
                            flw_dn_icop[i,Pref>Psurf[i,0]] = np.nan

                    if chem == 1:
                        ch4_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.ch4[i,::-1,0],Pref[::-1])[::-1]
                        co_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.co[i,::-1,0],Pref[::-1])[::-1]
                        h2o_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.h2o[i,::-1,0],Pref[::-1])[::-1]
                        co2_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.co2[i,::-1,0],Pref[::-1])[::-1]
                        nh3_icop[i,:] = interp.pchip_interpolate(sigma[::-1],output.nh3[i,::-1,0],Pref[::-1])[::-1]
                        if surf and mask_surf:
                            #when isobars intersect the surface, make values below surface undefined
                            ch4_icop[i,Pref>Psurf[i,0]] = np.nan
                            co_icop[i,Pref>Psurf[i,0]] = np.nan
                            h2o_icop[i,Pref>Psurf[i,0]] = np.nan
                            co2_icop[i,Pref>Psurf[i,0]] = np.nan
                            nh3_icop[i,Pref>Psurf[i,0]] = np.nan

            else:  # keep height at vertical coordinate
                Temp_icop[:,:] = Temp_icoh[:,:,0]
                Rho_icop[:,:] = output.Rho[:,:,0]
                Mh_icop[:,:,:] = output.Mh[:,:,:,0]
                W_icop[:,:] = W_icoh[:,:,0]
                # del_hseq_icop[:,:] = del_hseq[:,:,0]
                if RT == 1:
                    tau_sw_icop[:,:] = output.tau_sw[:,:,0]
                    tau_lw_icop[:,:] = output.tau_lw[:,:,0]
                    flw_up_icop[:,:] = flw_up_icoh[:,:,0]
                    flw_dn_icop[:,:] = flw_dn_icoh[:,:,0]
                if chem == 1:
                    ch4_icop[:,:] = output.ch4[:,:,0]
                    co_icop[:,:] = output.co[:,:,0]
                    h2o_icop[:,:] = output.h2o[:,:,0]
                    co2_icop[:,:] = output.co2[:,:,0]
                    nh3_icop[:,:] = output.nh3[:,:,0]

            # Convert icosahedral grid into lon-lat grid
            U_icop = (-Mh_icop[0]*np.sin(grid.lon[:,None])+Mh_icop[1]*np.cos(grid.lon[:,None]))/Rho_icop
            V_icop = (-Mh_icop[0]*np.sin(grid.lat[:,None])*np.cos(grid.lon[:,None])\
                      -Mh_icop[1]*np.sin(grid.lat[:,None])*np.sin(grid.lon[:,None])\
                      +Mh_icop[2]*np.cos(grid.lat[:,None]))/Rho_icop

            # interp_fun = interp.NearestNDInterpolator(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,Temp_icop[:,0])
            for lev in np.arange(d_sig):
                if type == 'gd' or type == 'GD':
                    if lev == 0:
                        print('Using "griddata" for horizontal interpolation')
                    # interp_fun.values = Temp_icop[:,lev]
                    # Temp_llp[:,:,lev] = interp_fun((loni,lati))
                    Temp_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,Temp_icop[:,lev],(loni,lati),method='nearest')
                    # interp_fun.values = Rho_icop[:,lev]
                    # Rho_llp[:,:,lev] = interp_fun((loni,lati))
                    Rho_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,Rho_icop[:,lev],(loni,lati),method='nearest')
                    # interp_fun.values = U_icop[:,lev]
                    # U_llp[:,:,lev] = interp_fun((loni,lati))
                    U_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,U_icop[:,lev],(loni,lati),method='nearest')
                    # interp_fun.values = V_icop[:,lev]
                    # V_llp[:,:,lev] = interp_fun((loni,lati))
                    V_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,V_icop[:,lev],(loni,lati),method='nearest')
                    # interp_fun.values = W_icop[:,lev]
                    # W_llp[:,:,lev] = interp_fun((loni,lati))
                    W_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,W_icop[:,lev]/Rho_icop[:,lev],(loni,lati),method='nearest')
                    # del_hseq_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,del_hseq_icop[:,lev],(loni,lati),method='nearest')
                    if RT == 1:
                        # interp_fun.values = tau_sw_icop[:,lev]
                        # tau_sw_llp[:,:,lev] = interp_fun((loni,lati))
                        tau_sw_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,tau_sw_icop[:,lev],(loni,lati),method='nearest')
                        # interp_fun.values = tau_lw_icop[:,lev]
                        # tau_lw_llp[:,:,lev] = interp_fun((loni,lati))
                        tau_lw_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,tau_lw_icop[:,lev],(loni,lati),method='nearest')
                        # interp_fun.values = flw_up_icop[:,lev]
                        # flw_up_llp[:,:,lev] = interp_fun((loni,lati))
                        flw_up_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,flw_up_icop[:,lev],(loni,lati),method='nearest')
                        # interp_fun.values = flw_dn_icop[:,lev]
                        # flw_dn_llp[:,:,lev] = interp_fun((loni,lati))
                        flw_dn_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,flw_dn_icop[:,lev],(loni,lati),method='nearest')
                        if lev == 0:
                            # interp_fun.values = output.Insol[:,0]
                            # insol_ll[:,:] = interp_fun((loni,lati))
                            insol_ll[:,:] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,output.Insol[:,0],(loni,lati),method='nearest')
                            if surf == 1:
                                # interp_fun.values = output.Tsurface[:,0]
                                # Tsurf_ll[:,:] = interp_fun((loni,lati))
                                Tsurf_ll[:,:] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,output.Tsurface[:,0],(loni,lati),method='nearest')
                    if chem == 1:
                        ch4_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,ch4_icop[:,lev],(loni,lati),method='nearest')
                        co_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,co_icop[:,lev],(loni,lati),method='nearest')
                        h2o_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,h2o_icop[:,lev],(loni,lati),method='nearest')
                        co2_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,co2_icop[:,lev],(loni,lati),method='nearest')
                        nh3_llp[:,:,lev] = interp.griddata(np.vstack([grid.lon*180/np.pi,grid.lat*180/np.pi]).T,nh3_icop[:,lev],(loni,lati),method='nearest')
                elif type == 'sh' or type == 'SH':
                    if lev == 0:
                        print('Using "SHExpandLSQ" for horizontal interpolation')
                    Temp_coeff[lev,:,:,:], Temp_chi2[lev] = chairs.expand.SHExpandLSQ(Temp_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                    Temp_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(Temp_coeff[lev,:,:,:],sampling=2))
                    Rho_coeff[lev,:,:,:], Rho_chi2[lev] = chairs.expand.SHExpandLSQ(Rho_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                    Rho_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(Rho_coeff[lev,:,:,:],sampling=2))
                    U_coeff[lev,:,:,:], U_chi2[lev] = chairs.expand.SHExpandLSQ(U_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                    U_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(U_coeff[lev,:,:,:],sampling=2))
                    V_coeff[lev,:,:,:], V_chi2[lev] = chairs.expand.SHExpandLSQ(V_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                    V_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(V_coeff[lev,:,:,:],sampling=2))
                    W_coeff[lev,:,:,:], W_chi2[lev] = chairs.expand.SHExpandLSQ(W_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                    W_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(W_coeff[lev,:,:,:],sampling=2))
                    if RT == 1:
                        tau_sw_coeff[lev,:,:,:], tau_sw_chi2[lev] = chairs.expand.SHExpandLSQ(tau_sw_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                        tau_sw_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(tau_sw_coeff[lev,:,:,:],sampling=2))
                        tau_lw_coeff[lev,:,:,:], tau_lw_chi2[lev] = chairs.expand.SHExpandLSQ(tau_lw_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                        tau_lw_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(tau_lw_coeff[lev,:,:,:],sampling=2))
                        flw_up_coeff[lev,:,:,:], flw_up_chi2[lev] = chairs.expand.SHExpandLSQ(flw_up_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                        flw_up_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(flw_up_coeff[lev,:,:,:],sampling=2))
                        flw_dn_coeff[lev,:,:,:], flw_dn_chi2[lev] = chairs.expand.SHExpandLSQ(flw_dn_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                        flw_dn_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(flw_dn_coeff[lev,:,:,:],sampling=2))
                        if lev == 0:
                            insol_coeff[lev,:,:,:], insol_chi2[lev] = chairs.expand.SHExpandLSQ(output.Insol[:,0],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                            insol_ll[:,:] = np.real(chairs.expand.MakeGridDH(insol_coeff[lev,:,:,:],sampling=2))
                            if surf == 1:
                                Tsurf_coeff[lev,:,:,:], Tsurf_chi2[lev] = chairs.expand.SHExpandLSQ(output.Tsurface[:,0],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                                Tsurf_ll[:,:] = np.real(chairs.expand.MakeGridDH(Tsurf_coeff[lev,:,:,:],sampling=2))

                    if chem == 1:
                        ch4_coeff[lev,:,:,:], ch4_chi2[lev] = chairs.expand.SHExpandLSQ(ch4_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                        ch4_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(ch4_coeff[lev,:,:,:],sampling=2))
                        co_coeff[lev,:,:,:], co_chi2[lev] = chairs.expand.SHExpandLSQ(co_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                        co_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(co_coeff[lev,:,:,:],sampling=2))
                        co2_coeff[lev,:,:,:], co2_chi2[lev] = chairs.expand.SHExpandLSQ(co2_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                        co2_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(co2_coeff[lev,:,:,:],sampling=2))
                        h2o_coeff[lev,:,:,:], h2o_chi2[lev] = chairs.expand.SHExpandLSQ(h2o_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                        h2o_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(h2o_coeff[lev,:,:,:],sampling=2))
                        nh3_coeff[lev,:,:,:], nh3_chi2[lev] = chairs.expand.SHExpandLSQ(nh3_icop[:,lev],grid.lat*180/np.pi,grid.lon*180/np.pi,lmax)
                        nh3_llp[:,:,lev] = np.real(chairs.expand.MakeGridDH(nh3_coeff[lev,:,:,:],sampling=2))
                elif type == 'spl' or type == 'SPL':
                    if lev == 0:
                        print('Using "SmoothSphereBivariateSpline" for horizontal interpolation')
                    theta = np.pi/2 - grid.lat
                    f = interp.SmoothSphereBivariateSpline(theta,grid.lon,Temp_icop[:,lev],s=3.5)
                    Temp_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                    f = interp.SmoothSphereBivariateSpline(theta,grid.lon,Rho_icop[:,lev],s=3.5)
                    Rho_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                    f = interp.SmoothSphereBivariateSpline(theta,grid.lon,U_icop[:,lev],s=5)
                    U_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                    f = interp.SmoothSphereBivariateSpline(theta,grid.lon,V_icop[:,lev],s=5)
                    V_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                    f = interp.SmoothSphereBivariateSpline(theta,grid.lon,W_icop[:,lev],s=5)
                    W_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                    if RT == 1:
                        f = interp.SmoothSphereBivariateSpline(theta,grid.lon,tau_sw_icop[:,lev],s=3.5)
                        tau_sw_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                        f = interp.SmoothSphereBivariateSpline(theta,grid.lon,tau_lw_icop[:,lev],s=3.5)
                        tau_lw_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                        f = interp.SmoothSphereBivariateSpline(theta,grid.lon,flw_up_icop[:,lev],s=3.5)
                        flw_up_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                        f = interp.SmoothSphereBivariateSpline(theta,grid.lon,flw_dn_icop[:,lev],s=3.5)
                        flw_dn_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                        if lev == 0:
                            f = interp.SmoothSphereBivariateSpline(theta,grid.lon,output.Insol[:,0],s=3.5)
                            insol_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                            if surf == 1:
                                f = interp.SmoothSphereBivariateSpline(theta,grid.lon,output.Tsurface[:,0],s=3.5)
                                Tsurf_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                    if chem == 1:
                        f = interp.SmoothSphereBivariateSpline(theta,grid.lon,ch4_icop[:,lev],s=3.5)
                        ch4_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                        f = interp.SmoothSphereBivariateSpline(theta,grid.lon,co_icop[:,lev],s=3.5)
                        co_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                        f = interp.SmoothSphereBivariateSpline(theta,grid.lon,co2_icop[:,lev],s=3.5)
                        co2_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                        f = interp.SmoothSphereBivariateSpline(theta,grid.lon,h2o_icop[:,lev],s=3.5)
                        h2o_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                        f = interp.SmoothSphereBivariateSpline(theta,grid.lon,nh3_icop[:,lev],s=3.5)
                        nh3_llp[:,:,lev] = f(np.pi/2-lat_range[::-1]*np.pi/180,lon_range*np.pi/180)
                else:
                    raise ValueError('Invalid type for horizontal interpolation in regrid.\n Valid types are "gd"/"GD" (griddata),"sh"/"SH" (spherical harmonics), "spl"/"SPL" (spherical spline).')

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
            # del_hseqf = openh5.create_dataset("del_hseq",data=del_hseq_llp,compression='gzip',compression_opts=comp)
            if type == 'sh' or type == 'SH':
                Tcoeff = openh5.create_dataset("Temp_coeff",data=Temp_coeff,compression='gzip',compression_opts=comp)
                Tchi2 = openh5.create_dataset("Temp_chi2",data=Temp_chi2,compression='gzip',compression_opts=comp)
                Rhocoeff = openh5.create_dataset("Rho_coeff",data=Rho_coeff,compression='gzip',compression_opts=comp)
                Rhochi2 = openh5.create_dataset("Rho_chi2",data=Rho_chi2,compression='gzip',compression_opts=comp)
                Ucoeff = openh5.create_dataset("U_coeff",data=U_coeff,compression='gzip',compression_opts=comp)
                Uchi2 = openh5.create_dataset("U_chi2",data=U_chi2,compression='gzip',compression_opts=comp)
                Vcoeff = openh5.create_dataset("V_coeff",data=V_coeff,compression='gzip',compression_opts=comp)
                Vchi2 = openh5.create_dataset("V_chi2",data=V_chi2,compression='gzip',compression_opts=comp)
                Wcoeff = openh5.create_dataset("W_coeff",data=W_coeff,compression='gzip',compression_opts=comp)
                Wchi2 = openh5.create_dataset("W_chi2",data=W_chi2,compression='gzip',compression_opts=comp)

            if RT == 1:
                tau_sw = openh5.create_dataset("tau_sw",data=tau_sw_llp,compression='gzip',compression_opts=comp)
                tau_lw = openh5.create_dataset("tau_lw",data=tau_lw_llp,compression='gzip',compression_opts=comp)
                flw_up = openh5.create_dataset("flw_up",data=flw_up_llp,compression='gzip',compression_opts=comp)
                flw_dn = openh5.create_dataset("flw_dn",data=flw_dn_llp,compression='gzip',compression_opts=comp)
                insol = openh5.create_dataset("insol",data=insol_ll,compression='gzip',compression_opts=comp)
                if surf == 1:
                    Tsurf = openh5.create_dataset("Tsurface",data=Tsurf_ll,compression='gzip',compression_opts=comp)

                if type == 'sh' or type == 'SH':
                    tauswcoeff = openh5.create_dataset("tau_sw_coeff",data=tau_sw_coeff,compression='gzip',compression_opts=comp)
                    tauswchi2 = openh5.create_dataset("tau_sw_chi2",data=tau_sw_chi2,compression='gzip',compression_opts=comp)
                    taulwcoeff = openh5.create_dataset("tau_lw_coeff",data=tau_lw_coeff,compression='gzip',compression_opts=comp)
                    taulwchi2 = openh5.create_dataset("tau_lw_chi2",data=tau_lw_chi2,compression='gzip',compression_opts=comp)
                    fnetupcoeff = openh5.create_dataset("flw_up_coeff",data=flw_up_coeff,compression='gzip',compression_opts=comp)
                    fnetupchi2 = openh5.create_dataset("flw_up_chi2",data=flw_up_chi2,compression='gzip',compression_opts=comp)
                    fnetdncoeff = openh5.create_dataset("flw_dn_coeff",data=flw_dn_coeff,compression='gzip',compression_opts=comp)
                    fnetdnchi2 = openh5.create_dataset("flw_dn_chi2",data=flw_dn_chi2,compression='gzip',compression_opts=comp)
                    insolcoeff = openh5.create_dataset("insol_coeff",data=insol_coeff,compression='gzip',compression_opts=comp)
                    insolchi2 = openh5.create_dataset("insol_chi2",data=insol_chi2,compression='gzip',compression_opts=comp)
                    if surf == 1:
                        Tsurfcoeff = openh5.create_dataset("Tsurf_coeff",data=Tsurf_coeff,compression='gzip',compression_opts=comp)
                        Tsurfchi2 = openh5.create_dataset("Tsurf_chi2",data=Tsurf_chi2,compression='gzip',compression_opts=comp)

            if chem == 1:
                ch4 = openh5.create_dataset("ch4",data=ch4_llp,compression='gzip',compression_opts=comp)
                co = openh5.create_dataset("co",data=co_llp,compression='gzip',compression_opts=comp)
                h2o = openh5.create_dataset("h2o",data=h2o_llp,compression='gzip',compression_opts=comp)
                co2 = openh5.create_dataset("co2",data=co2_llp,compression='gzip',compression_opts=comp)
                nh3 = openh5.create_dataset("nh3",data=nh3_llp,compression='gzip',compression_opts=comp)
                if type == 'sh' or type == 'SH':
                    ch4coeff = openh5.create_dataset("ch4_coeff",data=ch4_coeff,compression='gzip',compression_opts=comp)
                    ch4chi2 = openh5.create_dataset("ch4_chi2",data=ch4_chi2,compression='gzip',compression_opts=comp)
                    cocoeff = openh5.create_dataset("co_coeff",data=co_coeff,compression='gzip',compression_opts=comp)
                    cochi2 = openh5.create_dataset("co_chi2",data=co_chi2,compression='gzip',compression_opts=comp)
                    co2coeff = openh5.create_dataset("co2_coeff",data=co2_coeff,compression='gzip',compression_opts=comp)
                    co2chi2 = openh5.create_dataset("co2_chi2",data=co2_chi2,compression='gzip',compression_opts=comp)
                    h2ocoeff = openh5.create_dataset("h2o_coeff",data=h2o_coeff,compression='gzip',compression_opts=comp)
                    h2ochi2 = openh5.create_dataset("h2o_chi2",data=h2o_chi2,compression='gzip',compression_opts=comp)
                    nh3coeff = openh5.create_dataset("nh3_coeff",data=nh3_coeff,compression='gzip',compression_opts=comp)
                    nh3chi2 = openh5.create_dataset("nh3_chi2",data=nh3_chi2,compression='gzip',compression_opts=comp)

            openh5.close()

            ## calculate relative and potential vorticity and add to regrid file
            calc_RV_PV(grid,output,input,lon_range,lat_range,Pref,0,fileh5,pressure_vert,type=type,lmax_set=lmax_set)
            # calc_moc_streamf(grid,output,input,lon_range,lat_range,Pref,0,fileh5,pressure_vert)



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
