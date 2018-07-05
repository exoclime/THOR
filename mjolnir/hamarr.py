import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import os
import h5py
import time

plt.rcParams['image.cmap'] = 'magma'

class input:
    def __init__(self,resultsf,simID):
        fileh5 = resultsf+'/esp_output_'+simID+'.h5'
        if os.path.exists(fileh5):
            openh5 = h5py.File(fileh5)
        else:
            raise IOError(fileh5+' not found!')
        self.A = openh5['A'][...]  #'A' is the key for this dataset
        self.Rd = openh5['Rd'][...]  #[...] is syntax for "gimme all the data under this key"
        self.Omega = openh5['Omega'][...]
        self.Mmol = openh5['Mmol'][...]
        self.P_Ref = openh5['P_Ref'][...]
        self.Top_altitude = openh5['Top_altitude'][...]
        self.Cp = openh5['Cp'][...]
        self.resultsf = resultsf
        self.simID = simID
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
        self.nv = np.int(openh5['nv'][0])
        openh5.close()
        self.nvi = self.nv+1
        self.lon = (self.lonlat[::2])%(2*np.pi)
        self.lat = self.lonlat[1::2]

class output:
    def __init__(self,resultsf,simID,ntsi,nts,grid):
        # Initialize arrays
        self.Rho = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.Pressure = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.Mh = np.zeros((3,grid.point_num,grid.nv,nts-ntsi+1))
        self.Wh = np.zeros((grid.point_num,grid.nvi,nts-ntsi+1))
        self.ntsi = ntsi
        self.nts = nts

        # Read model results
        for t in np.arange(ntsi-1,nts):
            fileh5 = resultsf+'/esp_output_'+simID+'_'+np.str(t+1)+'.h5'
            if os.path.exists(fileh5):
                openh5 = h5py.File(fileh5)
            else:
                raise IOError(fileh5+' not found!')
            Rhoi = openh5['Rho'][...]
            Pressurei = openh5['Pressure'][...]
            Mhi = openh5['Mh'][...]
            Whi = openh5['Wh'][...]
            openh5.close()

            self.Rho[:,:,t-ntsi+1] = np.reshape(Rhoi,(grid.point_num,grid.nv))
            self.Pressure[:,:,t-ntsi+1] = np.reshape(Pressurei,(grid.point_num,grid.nv))
            self.Mh[0,:,:,t-ntsi+1] = np.reshape(Mhi[::3],(grid.point_num,grid.nv))
            self.Mh[1,:,:,t-ntsi+1] = np.reshape(Mhi[1::3],(grid.point_num,grid.nv))
            self.Mh[2,:,:,t-ntsi+1] = np.reshape(Mhi[2::3],(grid.point_num,grid.nv))

            self.Wh[:,:,t-ntsi+1] = np.reshape(Whi,(grid.point_num,grid.nvi))

def temperature(input,grid,output,sigmaref):
    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)

    # Set the latitude-longitude grid.
    res_deg = 0.005
    loni, lati = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2,np.pi/2,res_deg))
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1
    ###############
    # Temperature #
    ###############

    # Initialize arrays
    Temperatureii = np.zeros(grid.nv)
    Temperaturei = np.zeros((grid.point_num,d_sig))
    Temperature = np.zeros((d_lon[0],d_lon[1],d_sig,tsp))

    # Compute temperatures
    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            sigma = output.Pressure[i,:,t]
            for lev in np.arange(grid.nv):
                Temperatureii[lev] = output.Pressure[i,lev,t]/(input.Rd*output.Rho[i,lev,t])
            # Interpolate atmospheric column ot the reference pressure.
            # Need to reverse index all arrays because pchip only works if x is increasing
            Temperaturei[i,:] = interp.pchip_interpolate(sigma[::-1],Temperatureii[::-1],Pref[::-1])[::-1]
        # Convert icosahedral grid into lon-lat grid
        for lev in np.arange(d_sig):
            Temperature[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Temperaturei[:,lev],(loni,lati),method='nearest')

    del Temperatureii, Temperaturei

    # Averaging in time and longitude.
    if tsp > 1:
        Temperaturel = np.mean(Temperature[:,:,:,:],axis=1)
        del Temperature
        Temperaturelt = np.mean(Temperaturel[:,:,:],axis=2)
        del Temperaturel
    else:
        Temperaturelt = np.mean(Temperature[:,:,:,0],axis=1)
        del Temperature

    #################
    # Create figure #
    #################

    # Latitude
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)

    # Contour plot
    C = plt.contourf(latp*180/np.pi,Pref/100,Temperaturelt.T,40)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.gca().invert_yaxis()
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (mba)')
    clb = plt.colorbar(C)
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/temperature_ver_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def u(input,grid,output,sigmaref):
    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)

    # Set the latitude-longitude grid
    res_deg = 0.005
    loni, lati = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2,np.pi/2,res_deg))
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    ##################
    # Zonal momentum #
    ##################

    # Initialize arrays
    ZonalMii = np.zeros(grid.nv)
    ZonalMi = np.zeros((grid.point_num,d_sig))
    ZonalM = np.zeros((d_lon[0],d_lon[1],d_sig,tsp))

    ZonalMtemp = (output.Mh[0]*(-np.sin(grid.lon[:,None,None])) + \
                output.Mh[1]*np.cos(grid.lon[:,None,None]))/output.Rho
    # ZonalMtemp1 = np.vstack((ZonalMtemp[:,:,0],ZonalMtemp[:,:,1]))
    #
    # Ptemp = np.vstack((output.Pressure[:,:,0],output.Pressure[:,:,1]))
    #
    # y = np.arange(grid.point_num*tsp)
    # x = np.arange(grid.nv)
    # xx, yy = np.meshgrid(x,y)
    # fint2d = interp.interp2d(Ptemp,yy,ZonalMtemp1)
    # #ZonalMtemp2 = np.reshape(fint2d(Pref,Ptemp[:,0]),(grid.point_num,len(Pref),tsp))
    # ZonalMtemp2 = fint2d(Pref,Ptemp[:,0])

    # Compute zonal Winds
    beg = time.time()
    # v = np.hstack((output.Pressure,ZonalMtemp))
    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            sigma = output.Pressure[i,:,t]
            # for lev in np.arange(grid.nv):
            #     ZonalMii[lev] = (output.Mh[0,i,lev,t]*(-np.sin(grid.lon[i])) + \
            #                      output.Mh[1,i,lev,t]*np.cos(grid.lon[i]) + \
            #                      output.Mh[2,i,lev,t]*(0))/output.Rho[i,lev,t]
            # Interpolate atmospheric column to the reference pressure
            import pdb; pdb.set_trace()

            ZonalMi[i,:] = interp.pchip_interpolate(sigma[::-1],ZonalMtemp[i,::-1,t],Pref[::-1])[::-1]
        # Convert icosahedral grid into lon-lat grid
        for lev in np.arange(d_sig):
            ZonalM[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,ZonalMi[:,lev],(loni,lati),method='nearest')

    # plt.plot(Pref,ZonalMtemp2[0,:],'b-')
    # # plt.plot(Pref,ZonalMi[0,:],'r-')
    # plt.plot(output.Pressure[0,:,0],ZonalMtemp[0,:,0],'k.')
    # plt.plot(Ptemp[0,:],ZonalMtemp1[0,:],'r-')
    # plt.show()
    del ZonalMii, ZonalMi
    end = time.time()
    print(end-beg)

    # Averaging in time and longitude
    if tsp > 1:
        ZonalMl = np.mean(ZonalM[:,:,:,:],axis=1)
        del ZonalM
        ZonalMlt = np.mean(ZonalMl[:,:,:],axis=2)
        del ZonalMl
    else:
        ZonalMlt = np.mean(ZonalM[:,:,:,0],axis=1)
        del ZonalM

    #################
    # Create figure #
    #################

    # Latitude
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)

    # Contour plot
    C = plt.contourf(latp*180/np.pi,Pref/100,ZonalMlt.T,40)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.gca().invert_yaxis()
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (mba)')
    clb = plt.colorbar(C)
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/u_ver_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def uv_lev(input,grid,output,Plev):
    # Set the latitude-longitude grid.
    res_deg = 0.005
    loni, lati = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2,np.pi/2,res_deg))
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    #######################
    # Winds & Temperature #
    #######################

    # Initialize arrays
    Pr = np.zeros((grid.nv,1))
    Mx = np.zeros((grid.nv,1))
    My = np.zeros((grid.nv,1))
    Mz = np.zeros((grid.nv,1))
    Rhot = np.zeros((grid.point_num,1))
    Mxf = np.zeros((grid.point_num,1))
    Myf = np.zeros((grid.point_num,1))
    Mzf = np.zeros((grid.point_num,1))
    Ui = np.zeros((grid.point_num,1))
    Vi = np.zeros((grid.point_num,1))
    Uii = np.zeros((d_lon[0],d_lon[1],tsp))
    Vii = np.zeros((d_lon[0],d_lon[1],tsp))

    # Compute winds and temperatures
    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            for lev in np.arange(grid.nv):
                Pr[lev] = output.Pressure[i,lev,t]
                Mx[lev] = output.Mh[0,i,lev,t]
                My[lev] = output.Mh[1,i,lev,t]
                Mz[lev] = output.Mh[2,i,lev,t]
            # Interpolate in pressure
            Rhot[i] = interp.interp1d(Pr.T[0],output.Rho[i,:,t],kind='linear',fill_value='extrapolate')(Plev)
            Mxf[i] = interp.interp1d(Pr.T[0],Mx.T[0],kind='linear',fill_value='extrapolate')(Plev)
            Myf[i] = interp.interp1d(Pr.T[0],My.T[0],kind='linear',fill_value='extrapolate')(Plev)
            Mzf[i] = interp.interp1d(Pr.T[0],Mz.T[0],kind='linear',fill_value='extrapolate')(Plev)
        for i in np.arange(grid.point_num):
            Ui[i] = (Mxf[i]*(-np.sin(grid.lon[i])) + \
                     Myf[i]*np.cos(grid.lon[i]) + \
                     Mzf[i]*(0))/Rhot[i]
            Vi[i] = (Mxf[i]*(-np.sin(grid.lat[i])*np.cos(grid.lon[i])) + \
                     Myf[i]*(-np.sin(grid.lat[i])*np.sin(grid.lon[i])) + \
                     Mzf[i]*np.cos(grid.lat[i]))/Rhot[i]
        # Convert icosahedral grid into lon-lat grid
        Uii[:,:,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Ui.T[0],(loni,lati),method='cubic')
        Vii[:,:,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Vi.T[0],(loni,lati),method='cubic')

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
    lonp = np.arange(0,2*np.pi,res_deg)
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)
    C = plt.contourf(lonp*180/np.pi,latp*180/np.pi,Uiii,50)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')

    clb = plt.colorbar(C)
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/u_lev%#.3fmbar_i%d_l%d.pdf'%((Plev/100),output.ntsi,output.nts))
    plt.close()

    plt.figure()
    lonp = np.arange(0,2*np.pi,res_deg)
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)
    C = plt.contourf(lonp*180/np.pi,latp*180/np.pi,Viii,50)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')

    clb = plt.colorbar(C)
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/grid.figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/v_lev%#.3fmbar_i%d_l%d.pdf'%((Plev/100),output.ntsi,output.nts))
    plt.close()

def temperature_u_lev(input,grid,output,Plev):
    # Set the latitude-longitude grid.
    res_deg = 0.005
    loni, lati = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2,np.pi/2,res_deg))
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    #######################
    # Winds & Temperature #
    #######################

    # Initialize arrays
    Temp = np.zeros((grid.nv,1))
    Pr = np.zeros((grid.nv,1))
    Mx = np.zeros((grid.nv,1))
    My = np.zeros((grid.nv,1))
    Mz = np.zeros((grid.nv,1))
    Tempi = np.zeros((grid.point_num,1))
    Rhot = np.zeros((grid.point_num,1))
    Mxf = np.zeros((grid.point_num,1))
    Myf = np.zeros((grid.point_num,1))
    Mzf = np.zeros((grid.point_num,1))
    Ui = np.zeros((grid.point_num,1))
    Vi = np.zeros((grid.point_num,1))
    Temperatureii = np.zeros((d_lon[0],d_lon[1],tsp))
    Uii = np.zeros((d_lon[0],d_lon[1],tsp))
    Vii = np.zeros((d_lon[0],d_lon[1],tsp))

    # Compute winds and temperatures
    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            for lev in np.arange(grid.nv):
                Pr[lev] = output.Pressure[i,lev,t]
                Temp[lev] = output.Pressure[i,lev,t]/(input.Rd*output.Rho[i,lev,t])
                Mx[lev] = output.Mh[0,i,lev,t]
                My[lev] = output.Mh[1,i,lev,t]
                Mz[lev] = output.Mh[2,i,lev,t]
            # Interpolate in pressure
            Rhot[i] = interp.interp1d(Pr.T[0],output.Rho[i,:,t],kind='linear',fill_value='extrapolate')(Plev)
            Mxf[i] = interp.interp1d(Pr.T[0],Mx.T[0],kind='linear',fill_value='extrapolate')(Plev)
            Myf[i] = interp.interp1d(Pr.T[0],My.T[0],kind='linear',fill_value='extrapolate')(Plev)
            Mzf[i] = interp.interp1d(Pr.T[0],Mz.T[0],kind='linear',fill_value='extrapolate')(Plev)
            Tempi[i] = interp.interp1d(Pr.T[0],Temp.T[0],kind='linear',fill_value='extrapolate')(Plev)
        for i in np.arange(grid.point_num):
            Ui[i] = (Mxf[i]*(-np.sin(grid.lon[i])) + \
                     Myf[i]*np.cos(grid.lon[i]) + \
                     Mzf[i]*(0))/Rhot[i]
            Vi[i] = (Mxf[i]*(-np.sin(grid.lat[i])*np.cos(grid.lon[i])) + \
                     Myf[i]*(-np.sin(grid.lat[i])*np.sin(grid.lon[i])) + \
                     Mzf[i]*np.cos(grid.lat[i]))/Rhot[i]
        # Convert icosahedral grid into lon-lat grid
        # import pdb; pdb.set_trace()
        Temperatureii[:,:,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Tempi.T[0],(loni,lati),method='cubic')
        Uii[:,:,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Ui.T[0],(loni,lati),method='cubic')
        Vii[:,:,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Vi.T[0],(loni,lati),method='cubic')

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
    C = plt.contourf(lonp*180/np.pi,latp*180/np.pi,Temperature,50)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')
    plt.quiver(lonq*180/np.pi,latq*180/np.pi,U,V,color='0.5')

    clb = plt.colorbar(C)
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/temperature-uv_lev%#.3fmbar_i%d_l%d.pdf'%(Plev/100,output.ntsi,output.nts))
    plt.close()

def potential_temp(input,grid,output,sigmaref):
    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)
    kappa_ad = input.Rd/input.Cp  # adiabatic coefficient

    # Set the latitude-longitude grid.
    res_deg = 0.005
    loni, lati = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2,np.pi/2,res_deg))
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    #########################
    # Potential temperature #
    #########################

    # Initialize arrays
    Thetaii = np.zeros(grid.nv)
    Thetai = np.zeros((grid.point_num,d_sig))
    Theta = np.zeros((d_lon[0],d_lon[1],d_sig,tsp))

    # Compute temperatures
    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            sigma = output.Pressure[i,:,t]
            for lev in np.arange(grid.nv):
                Thetaii[lev] = output.Pressure[i,lev,t]/(input.Rd*output.Rho[i,lev,t]) * \
                                    (output.Pressure[i,lev,t]/input.P_Ref)**(-kappa_ad)
            # Interpolate atmospheric column ot the reference pressure.
            # Need to reverse index all arrays because pchip only works if x is increasing
            Thetai[i,:] = interp.pchip_interpolate(sigma[::-1],Thetaii[::-1],Pref[::-1])[::-1]
        # Convert icosahedral grid into lon-lat grid
        for lev in np.arange(d_sig):
            Theta[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,Thetai[:,lev],(loni,lati),method='nearest')

    del Thetaii, Thetai

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
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)

    # Contour plot
    C = plt.contourf(latp*180/np.pi,Pref/100,Thetalt.T,40,linewidths=None)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    plt.gca().invert_yaxis()
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (mba)')
    clb = plt.colorbar(C)
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/potential_temp_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()

def potential_temp_alt(Rho, Pressure, lon, lat, num, nts, ntsi, Rd, ps0, nv, sigmaref, Cp, Alt, resultsf):
    # Set the reference pressure
    Pref = ps0*sigmaref
    d_sig = np.size(sigmaref)
    kappa_ad = Rd/Cp  # adiabatic coefficient

    # Set the latitude-longitude grid.
    res_deg = 0.005
    loni, lati = np.meshgrid(np.arange(0,2*np.pi,res_deg),\
        np.arange(-np.pi/2,np.pi/2,res_deg))
    d_lon = np.shape(loni)
    tsp = nts-ntsi+1

    #########################
    # Potential temperature #
    #########################

    # Initialize arrays
    Thetaii = np.zeros(nv)
    Thetai = np.zeros((num,d_sig))
    Theta = np.zeros((d_lon[0],d_lon[1],d_sig,tsp))
    Altitudei = np.zeros((num,d_sig))
    Altitude = np.zeros((d_lon[0],d_lon[1],d_sig,tsp))

    # Compute temperatures
    for t in np.arange(tsp):
        for i in np.arange(num):
            sigma = Pressure[i,:,t]
            for lev in np.arange(nv):
                Thetaii[lev] = Pressure[i,lev,t]/(Rd*Rho[i,lev,t]) * \
                                    (Pressure[i,lev,t]/ps0)**(-kappa_ad)
            # Interpolate atmospheric column ot the reference pressure.
            # Need to reverse index all arrays because pchip only works if x is increasing
            Thetai[i,:] = interp.pchip_interpolate(sigma[::-1],Thetaii[::-1],Pref[::-1])[::-1]
            Altitudei[i,:] = interp.pchip_interpolate(sigma[::-1],Alt[::-1],Pref[::-1])[::-1]
        # Convert icosahedral grid into lon-lat grid
        for lev in np.arange(d_sig):
            Theta[:,:,lev,t] = interp.griddata(np.vstack([lon,lat]).T,Thetai[:,lev],(loni,lati),method='nearest')
            Altitude[:,:,lev,t] = interp.griddata(np.vstack([lon,lat]).T,Altitudei[:,lev],(loni,lati),method='nearest')

    del Thetaii, Thetai

    # Averaging in time and longitude.
    if tsp > 1:
        Thetal = np.mean(Theta[:,:,:,:],axis=1)
        Altitudel = np.mean(Altitude[:,:,:,:],axis=1)
        del Theta, Altitude
        Thetalt = np.mean(Thetal[:,:,:],axis=2)
        Altitudelt = np.mean(Altitudel[:,:,:],axis=2)
        del Thetal, Altitudel

    else:
        Thetalt = np.mean(Theta[:,:,:,0],axis=1)
        Altitudelt = np.mean(Altitude[:,:,:,0],axis=1)
        del Theta, Altitude

    #################
    # Create figure #
    #################

    # Latitude
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)

    # Contour plot
    C = plt.contourf(latp*180/np.pi,np.mean(Altitudelt,0),Thetalt.T,40,linewidths=None)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    # plt.gca().invert_yaxis()
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (mba)')
    clb = plt.colorbar(C)
    if not os.path.exists(resultsf+'/figures'):
        os.mkdir(resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(resultsf+'/figures/potential_temp_alt_i%d_l%d.pdf'%(ntsi,nts))
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
    C = plt.contourf(latp*180/np.pi,Pref/100,Philt.T,40,linewidths=None)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0

    plt.gca().invert_yaxis()
    plt.xlabel('Latitude (deg)')
    plt.ylabel('Pressure (mba)')
    #plt.ylabel('Altitude (m)')
    clb = plt.colorbar(C)
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
    C = plt.contourf(lonp*180/np.pi,latp*180/np.pi,Phi[:,:,0,0],40,linewidths=None)
    for cc in C.collections:
        cc.set_edgecolor("face") #fixes a stupid bug in matplotlib 2.0
    # plt.gca().invert_yaxis()
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Longitude (deg)')
    clb = plt.colorbar(C)
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.tight_layout()
    plt.savefig(input.resultsf+'/figures/pot_vort_lev%d_i%d_l%d.pdf'%(sigmaref/input.P_Ref*1000,output.ntsi,output.nts))
    plt.close()

def vring(input,grid,output,sigmaref):
    # Set the reference pressure
    Pref = input.P_Ref*sigmaref
    d_sig = np.size(sigmaref)

    # Set the latitude-longitude grid
    res_deg = 0.005
    loni, lati = np.meshgrid([np.pi/2,3*np.pi/2],\
        np.arange(-np.pi/2,np.pi/2,res_deg))
    d_lon = np.shape(loni)
    tsp = output.nts-output.ntsi+1

    # Get x wind at lon = 90 (terminator in transit)
    # Initialize arrays
    XMii = output.Mh[0,:,:,:]/output.Rho[:,:,:]
    XMi = np.zeros((grid.point_num,d_sig))
    XM = np.zeros((d_lon[0],d_lon[1],d_sig,tsp))

    # Compute zonal Winds
    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            sigma = output.Pressure[i,:,t]
            XMi[i,:] = interp.pchip_interpolate(sigma[::-1],XMii[i,::-1,t],Pref[::-1])[::-1]
        # Convert icosahedral grid into lon-lat grid
        for lev in np.arange(d_sig):
            XM[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,XMi[:,lev],(loni,lati),method='cubic')

    vrotX = input.Omega*input.A*np.cos(lati[:,0])
    import pdb; pdb.set_trace()

    plt.plot(lati,XM[:,0,0,0]-vrotX,'r.')
    plt.plot(lati,XM[:,1,0,0]+vrotX,'b.')
    plt.plot(lati,vrotX,'g.')
    plt.show()

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

        plt.semilogy(T,P/100,'k-',alpha= 0.5,lw=1)

    Tad = T[15]*(P/P[15])**kappa

    plt.plot(Tad,P/100,'r--')
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure (mbar)')
    plt.xlabel('Temperature [K]')

    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.savefig(input.resultsf+'/figures/TPprofile_i%d_l%d.pdf'%(output.ntsi,output.nts))
    plt.close()
