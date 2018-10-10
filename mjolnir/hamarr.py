import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import os
import h5py
import time

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
        self.Mmol = openh5['Mmol'][...]
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
            self.hstest = openh5['hstest'][...]
            if self.hstest[0] == 0:
                self.Tstar = openh5['Tstar'][...]
                # self.planet_star_dist = openh5['planet_star_dist'][...]
                # self.radius_star = openh5['radius_star'][...]
                # self.diff_fac = openh5['diff_fac'][...]
                # self.Tlow = openh5['Tlow'][...]
                # self.albedo = openh5['albedo'][...]
                # self.tausw = openh5['tausw'][...]
                # self.taulw = openh5['taulw'][...]
        if 'vulcan' in openh5.keys():
            self.vulcan = openh5['vulcan'][...]
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
    def __init__(self,resultsf,simID,ntsi,nts,grid):
        # Initialize arrays
        self.Rho = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.Pressure = np.zeros((grid.point_num,grid.nv,nts-ntsi+1))
        self.Mh = np.zeros((3,grid.point_num,grid.nv,nts-ntsi+1))
        self.Wh = np.zeros((grid.point_num,grid.nvi,nts-ntsi+1))
        self.ntsi = ntsi
        self.nts = nts
        self.time = np.zeros(nts-ntsi+1)
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

        # Read model results
        for t in np.arange(ntsi-1,nts):
            fileh5 = resultsf+'/esp_output_'+simID+'_'+np.str(t+1)+'.h5'
            if os.path.exists(fileh5):
                openh5 = h5py.File(fileh5)
            else:
                raise IOError(fileh5+' not found!')

            import pdb; pdb.set_trace()
            Rhoi = openh5['Rho'][...]
            Pressurei = openh5['Pressure'][...]
            Mhi = openh5['Mh'][...]
            Whi = openh5['Wh'][...]
            time = openh5['simulation_time'][0]/86400
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
            openh5.close()

            self.Rho[:,:,t-ntsi+1] = np.reshape(Rhoi,(grid.point_num,grid.nv))
            self.Pressure[:,:,t-ntsi+1] = np.reshape(Pressurei,(grid.point_num,grid.nv))
            self.Mh[0,:,:,t-ntsi+1] = np.reshape(Mhi[::3],(grid.point_num,grid.nv))
            self.Mh[1,:,:,t-ntsi+1] = np.reshape(Mhi[1::3],(grid.point_num,grid.nv))
            self.Mh[2,:,:,t-ntsi+1] = np.reshape(Mhi[2::3],(grid.point_num,grid.nv))
            self.Wh[:,:,t-ntsi+1] = np.reshape(Whi,(grid.point_num,grid.nvi))
            self.time[t-ntsi+1] = time
            if 'Etotali' in locals():
                self.Etotal[:,:,t-ntsi+1] = np.reshape(Etotali,(grid.point_num,grid.nv))
                self.Mass[:,:,t-ntsi+1] = np.reshape(Massi,(grid.point_num,grid.nv))
                self.AngMomx[:,:,t-ntsi+1] = np.reshape(AngMomxi,(grid.point_num,grid.nv))
                self.AngMomy[:,:,t-ntsi+1] = np.reshape(AngMomyi,(grid.point_num,grid.nv))
                self.AngMomz[:,:,t-ntsi+1] = np.reshape(AngMomzi,(grid.point_num,grid.nv))

class GetOutput:
    def __init__(self,resultsf,simID,ntsi,nts):
        self.input = input(resultsf,simID)
        self.grid = grid(resultsf,simID)
        self.output = output(resultsf,simID,ntsi,nts,self.grid)

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
    C = plt.contourf(latp*180/np.pi,Pref/1e5,Temperaturelt.T,40)
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

def u(input,grid,output,sigmaref):
    # contour spacing
    csp = 500

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
    # v = np.hstack((output.Pressure,ZonalMtemp))
    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            sigma = output.Pressure[i,:,t]
            # for lev in np.arange(grid.nv):
            #     ZonalMii[lev] = (output.Mh[0,i,lev,t]*(-np.sin(grid.lon[i])) + \
            #                      output.Mh[1,i,lev,t]*np.cos(grid.lon[i]) + \
            #                      output.Mh[2,i,lev,t]*(0))/output.Rho[i,lev,t]
            # Interpolate atmospheric column to the reference pressure
            # import pdb; pdb.set_trace()

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
    C = plt.contourf(latp*180/np.pi,Pref/1e5,ZonalMlt.T,40,cmap='viridis')

    levp = np.arange(np.ceil(np.min(ZonalMlt)/csp)*csp,np.floor(np.max(ZonalMlt)/csp)*csp,csp)
    c2 = plt.contour(latp*180/np.pi,Pref/1e5,ZonalMlt.T,levels=levp,colors='w',linewidths=1)
    plt.clabel(c2,inline=1,fontsize=10)
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

def w_ver(input,grid,output,sigmaref):
    # contour spacing
    csp = 200

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
    # Vert. momentum #
    ##################

    # Initialize arrays
    VertMii = np.zeros(grid.nv)
    VertMi = np.zeros((grid.point_num,d_sig))
    VertM = np.zeros((d_lon[0],d_lon[1],d_sig,tsp))

    # Wh is output at the vertical interfaces, so we must average
    VertMtemp = 0.5*(output.Wh[:,1:,:]+output.Wh[:,:-1,:])/output.Rho
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
    # v = np.hstack((output.Pressure,ZonalMtemp))
    for t in np.arange(tsp):
        for i in np.arange(grid.point_num):
            sigma = output.Pressure[i,:,t]
            # for lev in np.arange(grid.nv):
            #     ZonalMii[lev] = (output.Mh[0,i,lev,t]*(-np.sin(grid.lon[i])) + \
            #                      output.Mh[1,i,lev,t]*np.cos(grid.lon[i]) + \
            #                      output.Mh[2,i,lev,t]*(0))/output.Rho[i,lev,t]
            # Interpolate atmospheric column to the reference pressure
            # import pdb; pdb.set_trace()

            VertMi[i,:] = interp.pchip_interpolate(sigma[::-1],VertMtemp[i,::-1,t],Pref[::-1])[::-1]
        # Convert icosahedral grid into lon-lat grid
        for lev in np.arange(d_sig):
            VertM[:,:,lev,t] = interp.griddata(np.vstack([grid.lon,grid.lat]).T,VertMi[:,lev],(loni,lati),method='nearest')

    # plt.plot(Pref,ZonalMtemp2[0,:],'b-')
    # # plt.plot(Pref,ZonalMi[0,:],'r-')
    # plt.plot(output.Pressure[0,:,0],ZonalMtemp[0,:,0],'k.')
    # plt.plot(Ptemp[0,:],ZonalMtemp1[0,:],'r-')
    # plt.show()
    del VertMii, VertMi

    # Averaging in time and longitude
    if tsp > 1:
        VertMl = np.mean(VertM[:,:,:,:],axis=1)
        del ZonalM
        VertMlt = np.mean(VertMl[:,:,:],axis=2)
        del VertMl
    else:
        VertMlt = np.mean(VertM[:,:,:,0],axis=1)
        del VertM

    #################
    # Create figure #
    #################

    # Latitude
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)

    # Contour plot
    C = plt.contourf(latp*180/np.pi,Pref/1e5,VertMlt.T,40,cmap='viridis')

    levp = np.arange(np.ceil(np.min(VertMlt)/csp)*csp,np.floor(np.max(VertMlt)/csp)*csp,csp)
    c2 = plt.contour(latp*180/np.pi,Pref/1e5,VertMlt.T,levels=levp,colors='w',linewidths=1)
    plt.clabel(c2,inline=1,fontsize=10)
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
    C = plt.contourf(lonp*180/np.pi,latp*180/np.pi,Uiii,50,cmap='viridis')
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
    lonp = np.arange(0,2*np.pi,res_deg)
    latp = np.arange(-np.pi/2,np.pi/2,res_deg)
    C = plt.contourf(lonp*180/np.pi,latp*180/np.pi,Viii,50,cmap='viridis')
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
    #plt.plot(grid.lon*180/np.pi,grid.lat*180/np.pi,'w.')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Temperature (K)')

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
    C = plt.contourf(latp*180/np.pi,Pref/1e5,Thetalt.T,40,linewidths=None)
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
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    clb = plt.colorbar(C)
    clb.set_label(r'Potential temperature (K)')
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

        plt.semilogy(T,P/1e5,'k-',alpha= 0.5,lw=1)
        plt.plot(T[np.int(np.floor(grid.nv/2))],P[np.int(np.floor(grid.nv/2))]/100000,'r+',ms =5,alpha=0.5)
        plt.plot(T[np.int(np.floor(grid.nv*0.75))],P[np.int(np.floor(grid.nv*0.75))]/100000,'g+',ms =5,alpha=0.5)


    Tad = T[15]*(P/P[15])**kappa

    # plt.plot(Tad,P/100,'r--')
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure (bar)')
    plt.xlabel('Temperature [K]')
    plt.title('Time = %#.3f - %#.3f days'%(output.time[0],output.time[-1]))
    if not os.path.exists(input.resultsf+'/figures'):
        os.mkdir(input.resultsf+'/figures')
    plt.savefig(input.resultsf+'/figures/TPprofile_i%d_l%d.pdf'%(output.ntsi,output.nts))
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
    dz = grid.Altitude[1]-grid.Altitude[0]
    Vol = grid.areasT*dz

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

    output.Etotal = (Eint+Eg+Ek)*Vol[:,None,None]

    output.Mass = output.Rho*Vol[:,None,None]

    r = input.A+grid.Altitude

    rx = r[None,:,None]*np.cos(grid.lat[:,None,None])*np.cos(grid.lon[:,None,None])
    ry = r[None,:,None]*np.cos(grid.lat[:,None,None])*np.sin(grid.lon[:,None,None])
    rz = r[None,:,None]*np.sin(grid.lat[:,None,None])

    output.AngMomx = (ry*Mtotz-rz*Mtoty-output.Rho*input.Omega[0]*\
                     rz*r[None,:,None]*np.cos(grid.lat[:,None,None])*\
                     np.cos(grid.lon[:,None,None]))*Vol[:,None,None]
    output.AngMomy = (-rx*Mtotz+rz*Mtotx-output.Rho*input.Omega[0]*\
                     rz*r[None,:,None]*np.cos(grid.lat[:,None,None])*\
                     np.sin(grid.lon[:,None,None]))*Vol[:,None,None]
    output.AngMomz = (rx*Mtoty-ry*Mtotx + output.Rho*input.Omega[0]*\
                     r[None,:,None]**2*np.cos(grid.lat[:,None,None])*\
                     np.cos(grid.lat[:,None,None]))*Vol[:,None,None]

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
    dz = grid.Altitude[1]-grid.Altitude[0]
    Vol = grid.areasT*dz
    kappa = input.Rd/input.Cp

    potT = temperature*(input.P_Ref/output.Pressure)**kappa
    S = input.Cp * np.log(potT)
    output.Entropy = S*Vol[:,None,None]
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
