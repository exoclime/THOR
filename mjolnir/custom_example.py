#---- example of creating custom plots with the mjolnir library----------------
#---- code by Russell Deitrick ------------------------------------------------

import hamarr as ham
import matplotlib.pyplot as plt
import numpy as np

outall = ham.GetOutput('earth_hs','Earth',10,10,openrg=1)

input = outall.input
grid = outall.grid
output = outall.output
rg = outall.rg

sigmaref = np.linspace(input.P_Ref,np.float(50)*100,20)/input.P_Ref

fig, ax = plt.subplots(nrows=2, ncols=2)
fig.set_size_inches(16,14)

# plot zonal mean temperature
z = {'value':rg.Temperature, 'label':r'Temperature (K)', 'name':'temperature', 'cmap':'magma'}
ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=(0,360),save=False,axis=ax[0,0])
ax[0,0].set_ylim(1,0.05)  #manually adjust axis range

# plot zonal mean potential temperature (first we gotta calculate it)
kappa_ad = input.Rd/input.Cp  # adiabatic coefficient
pt = rg.Temperature*(rg.Pressure/input.P_Ref)**(-kappa_ad)
z = {'value':pt, 'label':r'Potential Temperature (K)', 'name':'PT', 'cmap':'magma'}
ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=(0,360),save=False,axis=ax[1,0])
ax[1,0].set_ylim(1,0.05)

# zonal mean zonal wind
z = {'value':rg.U, 'label':r'Zonal Wind Speed (m s$^{-1}$)', 'name':'u', 'cmap':'viridis'}
ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=(0,360),axis=ax[0,1],save=False)
ax[0,1].set_ylim(1,0.05)

# zonal mean vertical wind over one hemisphere
ham.vertical_lat(input,grid,output,rg,sigmaref,z,slice=(-90,90),axis=ax[1,1],save=False)
ax[1,1].set_ylim(1,0.05)

plt.savefig('custom_figure.pdf')
plt.close()
