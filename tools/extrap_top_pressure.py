#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# ---- code by Russell Deitrick ----------------------------------------------

"""Used with postprocessingsteps.py to extrapolate top boundary for
post-processing with RT modules (Alfrodull, e.g.). This will pad additional
vertical layers onto original grid, assuming isothermal profiles above original
top boundary. Adds layers until a user-set pressure is reached everywhere.
"""
import h5py as h5
import numpy as np
import shutil

def extrap_top_pressure(base_folder, continue_from_original, grid_file, planet_file, x2p):
    #open grid file, get grid resolutions, altitudes
    grid = h5.File(base_folder / grid_file,'r')

    nv = np.int(grid['nv'][0])
    point_num = np.int(grid['point_num'][0])
    z = grid['Altitude'][...]
    zh = grid['Altitudeh'][...]
    dz = zh[1:] - zh[:-1]
    grid.close()

    #open planet file just to get the goddamn gravity
    planet = h5.File(base_folder / planet_file,'r')

    g = planet['Gravit'][0]
    planet.close()

    #open last output file
    output = h5.File(base_folder / continue_from_original,'r')

    #to continue, we need (reshape: Mh, Pressure, Rho, Wh, Rd, Cp),
    #(copy directly: simulation_time, nstep, Tsurface)
    Mh = output['Mh'][...].reshape(point_num,nv,3)
    Pressure = output['Pressure'][...].reshape(point_num,nv)
    Rho = output['Rho'][...].reshape(point_num,nv)
    Rd = output['Rd'][...].reshape(point_num,nv)
    Cp = output['Cp'][...].reshape(point_num,nv)
    Wh = output['Wh'][...].reshape(point_num,nv+1)

    #index of column with max pressure at the top
    index_max_p = np.argmax(Pressure[:,-1])
    if Pressure[index_max_p,-1] <= x2p*1e5:
        print(f"Grid already extends past pressure = {x2p} bar, no extrapolation will be done.")
        continue_from_new = continue_from_original
        nv_new = nv
        top_altitude_new = zh[-1]

    else:
        #temperature
        T = Pressure / Rd / Rho
        if (Rd!=Rd[0]).any() or (Cp!=Cp[0]).any():
            print("WARNING: This extrapolation method cannot account for variations in Rd and Cp.")
            print("         The top layer values of Rd and Cp will be used for extrapolation")

        if (dz!=dz[0]).any():
            print("WARNING: This extrapolation method cannot account for variations in layer thickness")
            print("         The top layer thickness will be used for extrapolation")

        z_max_cols = z[-1] - Rd[:,-1]*T[:,-1]/g * np.log(x2p*1e5/Pressure[:,-1])
        z_max = np.max(z_max_cols)

        #extended section of altitude grid
        zh_ext = np.arange(zh[-1],z_max+dz[-1],dz[-1])
        z_ext = 0.5*(zh_ext[1:]+zh_ext[:-1])
        n_lay_ext = len(z_ext)            #number of additional layers to add
        nv_new = nv + n_lay_ext
        top_altitude_new = zh_ext[-1]

        #extrapolate pressure and density
        Pressure_ext = Pressure[:,-1][:,None] * np.exp(-g / (Rd[:,-1]*T[:,-1])[:,None] *
                                                (z_ext[None,:] - z[-1]))
        Rho_ext = Pressure_ext / Rd[:,-1][:,None] / T[:,-1][:,None]

        #pad wind fields with zero
        Mh_ext = np.zeros((point_num,n_lay_ext,3))
        Wh_ext = np.zeros((point_num,n_lay_ext))

        #fill Rd and Cp previous top values
        Rd_ext = np.zeros((point_num,n_lay_ext)) + Rd[:,-1][:,None]
        Cp_ext = np.zeros((point_num,n_lay_ext)) + Cp[:,-1][:,None]

        #glue the arrays together
        Pressure = np.concatenate([Pressure,Pressure_ext],axis=1)
        Rho = np.concatenate([Rho,Rho_ext],axis=1)
        Mh = np.concatenate([Mh,Mh_ext],axis=1)
        Wh = np.concatenate([Wh,Wh_ext],axis=1)
        Rd = np.concatenate([Rd,Rd_ext],axis=1)
        Cp = np.concatenate([Cp,Cp_ext],axis=1)

        #create h5 file
        new_folder = base_folder / "extrap_top_init"
        if not new_folder.exists():
            new_folder.mkdir()
        continue_from_new = new_folder / continue_from_original
        new_output = h5.File(continue_from_new,'w')
        new_output.create_dataset("Pressure", data=Pressure.ravel())
        new_output.create_dataset("Rho", data=Rho.ravel())
        new_output.create_dataset("Mh", data=Mh.ravel())
        new_output.create_dataset("Wh", data=Wh.ravel())
        new_output.create_dataset("Rd", data=Rd.ravel())
        new_output.create_dataset("Cp", data=Cp.ravel())
        new_output.create_dataset("simulation_time", data=output['simulation_time'])
        new_output.create_dataset("nstep", data=output['nstep'])
        if 'T_surface' in output.keys():
            new_output.create_dataset("T_surface", data=output['T_surface'][...])

        new_output.close()

        #copy planet file
        shutil.copyfile(base_folder / planet_file, new_folder / planet_file)
        planet = h5.File(new_folder / planet_file,'r+')
        planet['vlevel'][0] = nv_new
        planet['Top_altitude'][0] = top_altitude_new
        planet.close()

        #copy grid file
        shutil.copyfile(base_folder / grid_file, new_folder / grid_file)
        grid = h5.File(new_folder / grid_file,'r+')
        grid['nv'][0] = nv_new
        del grid['Altitude']
        del grid['Altitudeh']
        grid.create_dataset("Altitude", data=np.concatenate([z,z_ext]))
        grid.create_dataset("Altitudeh", data=np.concatenate([zh,zh_ext[1:]]))
        grid.close()

    output.close()

    return continue_from_new, nv_new, top_altitude_new
