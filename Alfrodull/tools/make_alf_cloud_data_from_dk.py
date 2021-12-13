"""Prepare a cloud data input fiole from a cloud file from DK"""

import argparse
from pathlib import Path

import h5py

import pathlib
import re

import h5py

import math
import numpy as np
import scipy

import pandas as pd

from scipy.interpolate import interp1d
from scipy.integrate import simps, quad

from astropy import constants as const
from astropy.modeling import models

from astropy import units as u

from alfrodull_input_tools import load_opacities_wavelength_bins


parser = argparse.ArgumentParser(description='Alfrodull cloud data file builder from DK file')

parser.add_argument("source_file", action="store", help = "input file to read data from, in DK units")
parser.add_argument("output_file", action="store", help = "output file to write data in SI units")


parser.add_argument("-o", "--opacities", action="store", default="input/opac_sample_r5.h5", help = "opacities file used as souorce of binning")
parser.add_argument('-c', '--CGS', action='store_true', default=False, help="opacities file is in CGS")


args = parser.parse_args()

cloud_input = Path(args.source_file)
output_file = Path(args.output_file)


# cloud input file from daniel

#cloud_input = pathlib.Path("../input/cross_sections_mgsio3_r1mu.dat")
#opacity_sample = pathlib.Path( "../input/opac_sample_r5.h5" )
#cloud_sample_output = pathlib.Path( "../input/cloud_sample_r5.h5" )


# opacity input to use for
opacity_samples_file = Path( args.opacities )
opacity_samples_is_CGS = args.CGS

if not opacity_samples_file.exists():
    print("Error: need to specify an opacity file with -o to read wavelength binning for Alf")
    exit(-1)

if opacity_samples_is_CGS:
    print(f"Reading opacities in CGS from {opacity_samples_file}")
else:
    print(f"Reading opacities in SI from {opacity_samples_file}")

            
opac_bin_centers_wavelengths, opac_wavelength_interfaces = load_opacities_wavelength_bins(opacity_samples_file)

# change wavelength unit to SI
if opacity_samples_is_CGS:
    opac_wavelengths_centers *= 1e-2
    opac_wavelength_interfaces *= 1e-2

print(f"Found {len(opac_bin_centers_wavelengths)} opacity wavelength bins")
print(f"from {opac_wavelength_interfaces[0]} m to {opac_wavelength_interfaces[-1]} m ")

# load cloud data
#wavelengths (mu)
#size parameter  
#extinction cross section (cm^2)  
#scattering cross section (cm^2)  
#absorption cross section (cm^2)  
#single scattering albedo         
#asymmetry parameter  

names = [
"wavelengths",
"size_parameter",
"extinction_cross_section",
"scattering_cross_section",
"absorption_cross_section",
"single_scattering_albedo",
"asymmetry_parameter" ]

cloud_input_data = pd.read_csv(cloud_input, sep='\s+', header=None, skiprows=1, names=names)

# load and conver to SI
print("Loading cloud data and converting")
cloud_wavelength = cloud_input_data['wavelengths'].values*1e-6
cloud_absorption = cloud_input_data['absorption_cross_section'].values*1e-4
cloud_scattering = cloud_input_data["scattering_cross_section"].values*1e-4

cloud_asymmetry = cloud_input_data['asymmetry_parameter'].values

print(f"Found {len(cloud_wavelength)} bins")
print(f"from {cloud_wavelength[0]} m to {cloud_wavelength[-1]} m")

# interpolating data
interpkind = 'cubic'

cloud_absorption_out = interp1d(cloud_wavelength, 
                                 cloud_absorption, 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(opac_bin_centers_wavelengths)

cloud_absorption_log_out = np.exp(interp1d(cloud_wavelength, 
                                 np.log(cloud_absorption), 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(opac_bin_centers_wavelengths))

cloud_absorption_loglog_out = np.exp(interp1d(np.log(cloud_wavelength), 
                                 np.log(cloud_absorption), 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(np.log(opac_bin_centers_wavelengths)))

cloud_absorption_loglin_out = interp1d(np.log(cloud_wavelength), 
                                 cloud_absorption, 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(np.log(opac_bin_centers_wavelengths))


cloud_scattering_out = interp1d(cloud_wavelength, 
                                 cloud_scattering, 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(opac_bin_centers_wavelengths)

cloud_scattering_log_out = np.exp(interp1d(cloud_wavelength, 
                                 np.log(cloud_scattering), 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(opac_bin_centers_wavelengths))
                                     
cloud_scattering_loglog_out = np.exp(interp1d(np.log(cloud_wavelength), 
                                 np.log(cloud_scattering), 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(np.log(opac_bin_centers_wavelengths)))

cloud_scattering_loglin_out = interp1d(np.log(cloud_wavelength), 
                                 cloud_scattering, 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(np.log(opac_bin_centers_wavelengths))

cloud_asymmetry_out = interp1d(cloud_wavelength, 
                                 cloud_asymmetry, 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(opac_bin_centers_wavelengths)

cloud_asymmetry_loglin_out = interp1d(np.log(cloud_wavelength), 
                                 cloud_asymmetry, 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(np.log(opac_bin_centers_wavelengths))

cloud_asymmetry_loglog_out = np.exp(interp1d(np.log(cloud_wavelength), 
                                 np.log(cloud_asymmetry), 
                                 bounds_error=False,
                                 kind=interpkind,
                                 fill_value='extrapolate' )(np.log(opac_bin_centers_wavelengths)))

with h5py.File(output_file, "w") as f:
    lambda_dset = f.create_dataset("wavelength", (len(opac_bin_centers_wavelengths),), dtype=np.float64)
    absorption_dset = f.create_dataset("absorption", cloud_absorption_loglog_out.shape, dtype=np.float64)
    scattering_dset = f.create_dataset("scattering", cloud_scattering_loglog_out.shape, dtype=np.float64)
    asymmetry_dset = f.create_dataset("asymmetry", cloud_asymmetry_loglog_out.shape, dtype=np.float64)
    lambda_dset[...] = opac_bin_centers_wavelengths
    absorption_dset[...] = cloud_absorption_loglog_out
    scattering_dset[...] = cloud_scattering_loglog_out
    asymmetry_dset[...] = cloud_asymmetry_loglog_out
