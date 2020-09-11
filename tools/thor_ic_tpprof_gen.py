"""Quick script to create the vertical altitude/pressure/density
profiles needed for setting initial conditions for a THOR run.  Just
set the several variables below (including specifiying the name of the
desired profile file), run the script, and you'll create an
approximately-valid set of profiles.


2020-06-17 09:45 IJMC: Created and uploaded to THOR.
"""


import numpy as np

#### Set initial parameters:

# Parameters for Earth:
rplanet = 1.0 * 6378136   # meters  
mplanet = 1.0 * 5.974e24 # kg
pmin, psurf = 1e-6, 1    # bars
MMW = 30                  # atomic mass units
temp = 270                # typical atmospheric temperature, in Kelvin
nlayers = 1000
outputfile = 'TPprof_Earth_simple.dat'


# Parameters for a hot Neputune:
rplanet = 4.72 * 6378136   # meters  
mplanet = 29.32 * 5.974e24 # kg
pmin, psurf = 1e-4, 300    # bars
MMW = 2.3                  # atomic mass units
temp = 2000                # typical atmospheric temperature, in Kelvin
nlayers = 1000
outputfile = 'TPprof_ltt9779b_simple_P300-1e-4.dat'


#### Set up grids and constants:
pressure = 1e5 * np.logspace(np.log10(psurf), np.log10(pmin), nlayers+1) # Pa
gsurf = 6.673e-11*mplanet/rplanet**2
scaleheight = 1.38e-23 * temp / (gsurf * MMW * 1.67e-27) # in meters
altitude = np.abs(-scaleheight * np.log(pressure/pressure.max()))  # in meters


#### Calculate density via hydrostatic equilibrium:
dpdz = np.diff(pressure)/np.diff(altitude)
density = -dpdz * (rplanet + altitude[0:-1])**2/(6.673e-11 * mplanet) # SI units


#### Write the file to disk:
f = open(outputfile, 'w')
f.write('Height Pressure Rho\n')
for ii in range(nlayers):
    layervals = (altitude[ii], pressure[ii], density[ii])
    f.write('%1.6e %1.6e %1.6e\n' % layervals)

f.close()

