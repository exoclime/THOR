"""Functions to prepare alfrodull input"""

import h5py

import pathlib
import re

import h5py
import imageio
import IPython.display as disp

import math
import numpy as np
import scipy

import pandas as pd

from scipy.interpolate import interp1d
from scipy.integrate import simps, quad

from astropy import constants as const
from astropy.modeling import models

from astropy import units as u


def load_opacities_wavelength_bins(opacity_sample):
    """Load wavelength bins from Alfrodull/HELIOS opacity sample file,
    returns center and interface values"""
    # get wavelength bin centers
    with h5py.File(opacity_sample, "r") as opac_h5:

        # wavelength grid
        try:
            opac_wave = [x for x in opac_h5["center wavelengths"][:]]
        except KeyError:
            opac_wave = [x for x in opac_h5["wavelengths"][:]]
        nbin = len(opac_wave)


        # interface positions of the wavelength bins
        try:
            opac_interwave = [i for i in opac_h5["interface wavelengths"][:]]
        except KeyError:
            # quick and dirty way to get the lamda interface values
            opac_interwave = []
            opac_interwave.append(opac_wave[0] - (opac_wave[1] - opac_wave[0])/2)
            for x in range(len(opac_wave) - 1):
                opac_interwave.append((opac_wave[x+1] + opac_wave[x])/2)
            opac_interwave.append(opac_wave[-1] + (opac_wave[-1] - opac_wave[-2])/2)

            # widths of the wavelength bins
            try:
                opac_deltawave = [w for w in opac_h5["wavelength width of bins"][:]]
            except KeyError:
                opac_deltawave = []
                for x in range(len(opac_interwave) - 1):
                    opac_deltawave.append(opac_interwave[x + 1] - opac_interwave[x])

        return np.array(opac_wave), np.array(opac_interwave)

def convert_helios_opacity_to_SI_units(opacity_samples_in, opacity_samples_out):
    """Load HELIOS opacity file and convert from CGS to SI units"""

    temperatures_unit_conv = 1.0
    pressures_unit_conv    = 1.0e-1
    wavelength_unit_conv   = 1.0e-2
    opacity_unit_conv = 1.0e-1
    scat_unit_conv = 1.0e-4

    with h5py.File(opacity_samples_in, "r") as opac_cgs_h5:
        with h5py.File(opacity_samples_out, "w") as opac_si_h5:
            if 'center wavelengths' in opac_cgs_h5:
                opac_cgs_lambda = opac_cgs_h5['center wavelengths'][...]
            else:
                opac_cgs_lambda = opac_cgs_h5['wavelengths'][...]

            lambda_dset = opac_si_h5.create_dataset("wavelengths", (opac_cgs_lambda.shape), dtype=np.float64)
            lambda_dset[...] = opac_cgs_lambda * wavelength_unit_conv
            lambda_dset.attrs['units'] = 'm'

            if 'ypoints' in opac_cgs_h5:
                opac_cgs_ypoints = opac_cgs_h5['ypoints'][...]

                ypoints_dset = opac_si_h5.create_dataset("ypoints", (opac_cgs_ypoints.shape), dtype=np.float64)
                ypoints_dset[...] = opac_cgs_ypoints
                ypoints_dset.attrs['units'] = '-'

            if "interface wavelengths" in opac_cgs_h5:
                opac_cgs_interwave = opac_cgs_h5["interface wavelengths"][...]

                interwave_dset = opac_si_h5.create_dataset("interface wavelengths", (opac_cgs_interwave.shape), dtype=np.float64)
                interwave_dset[...] = opac_cgs_interwave*wavelength_unit_conv
                interwave_dset.attrs['units'] = 'm'

            if "wavelength width of bins" in opac_cgs_h5:
                opac_cgs_deltawave = opac_cgs_h5["wavelength width of bins"][...]

                deltawave_dset = opac_si_h5.create_dataset("wavelength width of bins", (opac_cgs_deltawave.shape), dtype=np.float64)
                deltawave_dset[...] = opac_cgs_deltawave*wavelength_unit_conv
                deltawave_dset.attrs['units'] = 'm'

            opac_cgs_meanmolmass = opac_cgs_h5["meanmolmass"][...]
            meanmolmass_dset = opac_si_h5.create_dataset("meanmolmass", (opac_cgs_meanmolmass.shape), dtype=np.float64)
            meanmolmass_dset[...] = opac_cgs_meanmolmass
            meanmolmass_dset.attrs['units'] = 'amu'

            opac_cgs_pressures = opac_cgs_h5["pressures"][...]
            pressures_dset = opac_si_h5.create_dataset("pressures", (opac_cgs_pressures.shape), dtype=np.float64)
            pressures_dset[...] = opac_cgs_pressures*pressures_unit_conv
            pressures_dset.attrs['units'] = 'Pa'

            opac_cgs_temperatures = opac_cgs_h5["temperatures"][...]
            temperatures_dset = opac_si_h5.create_dataset("temperatures", (opac_cgs_temperatures.shape), dtype=np.float64)
            temperatures_dset[...] = opac_cgs_temperatures*temperatures_unit_conv
            temperatures_dset.attrs['units'] = 'K'

            opac_cgs_kpoints = opac_cgs_h5["kpoints"][...]
            kpoints_dset = opac_si_h5.create_dataset("kpoints", (opac_cgs_kpoints.shape), dtype=np.float64)
            kpoints_dset[...] = opac_cgs_kpoints*opacity_unit_conv
            kpoints_dset.attrs['units'] = 'm^2.kg^-1'

            opac_cgs_scat = opac_cgs_h5["weighted Rayleigh cross-sections"][...]
            scat_dset = opac_si_h5.create_dataset("weighted Rayleigh cross-sections", (opac_cgs_scat.shape), dtype=np.float64)
            scat_dset[...] = opac_cgs_scat*scat_unit_conv
            scat_dset.attrs['units'] = 'm^2'

def calc_planck(lamda, temp):
    """ calculates the Planckian blackbody function at a given wavelength and temperature """
    C = const.c.value                   # speed of light
    K_B = const.k_B.value               # Boltzmann constant
    H = const.h.value                   # Planck constant


    term1 = 2 * H * C**2 / lamda**5

    term2 = np.exp(H * C / (lamda * K_B * temp)) - 1

    result = term1 * 1 / term2

    return result


# interpolate spectrum
def rebin_spectrum(new_interfaces, orig_lambda, orig_spectrum, Teff):
    """Compute spectrum over binning for Alfrodull by integrating original high-res spectrum to lower res"""

    # integrand to fill missing values at extremes
    def integrand(x):
        return calc_planck(x, Teff)

    # interpolator for original spectrum
    spectrum_interpolator = interp1d(orig_lambda,
                                 orig_spectrum,
                                 bounds_error=False,
                                 fill_value=(orig_spectrum[0], orig_spectrum[-1]) )

    # integrate spectrum in bins

    l_start = new_interfaces[0]

    intg = []

    # avoid the last points that look funny in the spectrum
    for l_end in new_interfaces[1:]:
        integrate_idx_start = np.searchsorted(orig_lambda, l_start, side='right')
        integrate_idx_stop = np.searchsorted(orig_lambda, l_end, side='left')
        if integrate_idx_start == integrate_idx_stop or integrate_idx_stop >= len(orig_lambda[:-1]):
            # use planck functions to fill missing points
            intg.append(quad(integrand, l_start, l_end)[0]/(l_end-l_start)*math.pi)
        else:
            # use simpson integration
            # load datapoints from specctrum
            x_s = orig_lambda[integrate_idx_start:integrate_idx_stop]
            y_s = orig_spectrum[integrate_idx_start:integrate_idx_stop]

            # interpolate end points
            x_s[0] = l_start
            x_s[-1] = l_end
            y_s[0] = spectrum_interpolator(l_start)
            y_s[-1] = spectrum_interpolator(l_end)

            # integrate
            integral = simps(y_s, x_s)
            intg.append(integral/(l_end - l_start))

        l_start = l_end

    return np.array(intg)
