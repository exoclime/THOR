"""Command line script to download a phoenix spectrum and output it to a file for Alfrodull"""


from pathlib import Path
import numpy as np

from phoenix import create_interpolated_spectrum
from alfrodull_input_tools import load_opacities_wavelength_bins, rebin_spectrum

import h5py
import argparse

parser = argparse.ArgumentParser(description='Alfrodull unit converter')

parser.add_argument("output_file", action="store", default="output_spectrum.h5", help = "input file to read data from, in CGS units")
parser.add_argument("-i", "--intermediate", action="store", default=None, help = "store intermediate interpolated spectrum with high res binning to this file")

parser.add_argument("-o", "--opacities", action="store", default="input/opac_sample_r5.h5", help = "opacities file used as souorce of binning")

parser.add_argument('-c', '--CGS', action='store_true', default=False, help="opacities file is in CGS")

parser.add_argument('-t', '--teff',
                    default=4970, action='store',
                    type=float,
                    help="Teff")

parser.add_argument('-l', '--logg',
                    default=3.32, action='store',
                    type=float,
                    help="log(g)")

parser.add_argument('-Z', '--FeH',
                    default=0.03, action='store',
                    type=float,
                    help="Fe/H")

parser.add_argument('-a', '--Alpha',
                    default=0.0, action='store',
                    type=float,
                    help="alpha/M")


args = parser.parse_args()



# target data
Teff = args.teff
logg = args.logg
FeH = args.FeH
alphaM = args.Alpha

print("Preparing stellar spectrum from Phoenix with:")
print(f"Teff: {Teff} K")
print(f"log(g): {logg}")
print(f"Fe/H: {FeH}")
print(f"alpha/M: {alphaM}")
print()



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

# output
resampled_stellar_spectrum_path = args.output_file


output_phoenix_interpolated_spectrum = args.intermediate is not None

wavelengths, spectrum = create_interpolated_spectrum(Teff, logg, FeH, alphaM)

# to output the phoenix spectrum
if output_phoenix_interpolated_spectrum:
    print(f"Output intermediate HiRes spectrum to {args.intermediate}")
    with h5py.File(args.intermediate, "w") as f:
        lambda_dset = f.create_dataset("wavelength", (len(wavelengths),), dtype=np.float64)
        lambda_dset.attrs['units'] = 'm'
        flux_dset = f.create_dataset("flux", spectrum.shape, dtype=np.float64)
        flux_dset.attrs['units'] = 'J/s/m^2/m'
        lambda_dset[...] = wavelengths
        flux_dset[...] = spectrum


        
opac_wavelengths_centers, opac_wavelength_interfaces = load_opacities_wavelength_bins(opacity_samples_file)

# change wavelength unit to SI
if opacity_samples_is_CGS:
    opac_wavelengths_centers *= 1e-2
    opac_wavelength_interfaces *= 1e-2

print(f"from {opac_wavelength_interfaces[0]} m to {opac_wavelength_interfaces[-1]} m ")

rebinned_spectrum = rebin_spectrum(opac_wavelength_interfaces, wavelengths, spectrum, Teff)

print(f"Output to {resampled_stellar_spectrum_path}")
with h5py.File(resampled_stellar_spectrum_path, "w") as f:
    lambda_dset = f.create_dataset("wavelength", (len(opac_wavelengths_centers),), dtype=np.float64)
    flux_dset = f.create_dataset("flux", rebinned_spectrum.shape, dtype=np.float64)
    lambda_dset[...] = opac_wavelengths_centers
    flux_dset[...] = rebinned_spectrum

    Tstar_dset = f.create_dataset("T_star", (1,), dtype=np.float64)
    Tstar_dset[...] = Teff
    
    logg_dset = f.create_dataset("logg", (1,), dtype=np.float64)
    logg_dset[...] = logg

    FeH_dset = f.create_dataset("FeH", (1,), dtype=np.float64)
    FeH_dset[...] = FeH

    alphaM_dset = f.create_dataset("alphaM", (1,), dtype=np.float64)
    alphaM_dset[...] = alphaM
