"""Load and interpolate spectrum from PHOENIX. Convert unit to SI.
"""
from ftplib import FTP
import tempfile

from pathlib import Path
import numpy as np

from astropy.io import fits

from scipy.interpolate import interpn

import h5py

# Data ranges from website
Teffs = list(range(2300, 7001, 100)) + list(range(7000, 12001, 200))
loggs = [ x/10 for x in range(0, 61, 5)]  # loggs (unnamed param in filename) [0.0 -> 6.0]
FeHs = [ -4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0] # Z [-4.0 -> +1.0]
alphaMs = [-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2] # [-0.2 -> +1.2]


def get_interp_indices(val, lst):
    """Find indices of elements around value in sorted array, returning same index if value is one of the limit"""
    npl = np.array(lst)
    lo = np.searchsorted(npl, val, side='right', sorter = None)-1
    hi = np.searchsorted(npl, val, side='left', sorter = None)
    return lo, hi

def makeftppath(Z, alpha, logg, teff):
    """Prepare path for download from PHOENIX"""
    if alpha == 0.0: 
        if Z == 0.0:
            basedir = f"Z-{Z:1.1f}"
            filename = f"lte{teff:05d}-{logg:1.2f}-{Z:1.1f}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
        else:
            basedir = f"Z{Z:+1.1f}"
            filename = f"lte{teff:05d}-{logg:1.2f}{Z:+1.1f}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
    else:
        basedir = f"Z-{abs(Z):1.1f}.Alpha={alpha:+1.2f}"
        filename = f"lte{teff:05d}-{logg:1.2f}{Z:+1.1f}.Alpha={alpha:+1.2f}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
    
    return basedir, filename

def get_files(files, tempdir):
    """Get files from PHOENIX ftp site:
    * wavelength grid
    * grid files prepared with makeftppath"""
    with FTP("phoenix.astro.physik.uni-goettingen.de") as ftp:
        ftp.login()

        ftp.cwd("HiResFITS")
        #ftp.retrlines('LIST') 
        wavelength_grid = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
        with open(tempdir / wavelength_grid, "wb") as fp:
            print(f"Retrieving {wavelength_grid} from PHOENIX ftp")
            ftp.retrbinary("RETR " + wavelength_grid, fp.write)
            
            
        datadir = "PHOENIX-ACES-AGSS-COND-2011"

        for d,f in files:
            ftp.cwd(f"{datadir}/"+d)
            with open(tempdir / f, "wb") as fp:
                print(f"Retrieving {d + '/' + f} from PHOENIX ftp")
                ftp.retrbinary("RETR " + f, fp.write)
            ftp.cwd("../../")


def interpolate_dims(spectra_in, rng, var_x, coord_x):
    """Interpolate from list of 'spectra_in' between i and i+rng spectrum. coord_x is the value to interpolate to, var_x the values of both edges of interpolation"""
    spectra_out = []
    
    for i in range(rng):
        if var_x[0] == var_x[1]:
            spectra_out.append(spectra_in[i])
        else:
            spectra_out.append((spectra_in[i+rng] - spectra_in[i])/(var_x[1] - var_x[0])*(coord_x - var_x[0])+spectra_in[i] )

    return spectra_out


def create_interpolated_spectrum(Teff, logg, FeH, alphaM):
    """Create spectrum from PHOENIX, interpolated between stellar characteristics"""
    # get indices to download 
    alphaM_lohi = get_interp_indices(alphaM, alphaMs)
    FeH_lohi = get_interp_indices(FeH, FeHs)
    logg_lohi = get_interp_indices(logg, loggs)
    Teff_lohi = get_interp_indices(Teff, Teffs)

    # files to download, for each pair of ppoints for each dimensions
    files = [
        makeftppath(FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[0]], Teffs[Teff_lohi[0]]),
        makeftppath(FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[0]], Teffs[Teff_lohi[1]]),
        makeftppath(FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[1]], Teffs[Teff_lohi[0]]),
        makeftppath(FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[1]], Teffs[Teff_lohi[1]]),
        makeftppath(FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[0]], Teffs[Teff_lohi[0]]),
        makeftppath(FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[0]], Teffs[Teff_lohi[1]]),
        makeftppath(FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[1]], Teffs[Teff_lohi[0]]),
        makeftppath(FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[1]], Teffs[Teff_lohi[1]]),
        
        makeftppath(FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[0]], Teffs[Teff_lohi[0]]),
        makeftppath(FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[0]], Teffs[Teff_lohi[1]]),
        makeftppath(FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[1]], Teffs[Teff_lohi[0]]),
        makeftppath(FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[1]], Teffs[Teff_lohi[1]]),
        makeftppath(FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[0]], Teffs[Teff_lohi[0]]),
        makeftppath(FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[0]], Teffs[Teff_lohi[1]]),
        makeftppath(FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[1]], Teffs[Teff_lohi[0]]),
        makeftppath(FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[1]], Teffs[Teff_lohi[1]])
    ]

    with tempfile.TemporaryDirectory() as tmpdirname:
        # use this instead of context manager to keep temporary dir
        # tmpdirname = tempfile.mkdtemp()
    
        print(f"Downloading to temp dir: {tmpdirname}")
    
        tempdir = Path(tmpdirname)
        get_files(files, tempdir)
        
        # interpolation points. represents coordinates of hypercube we want to interpolate in 
        grid = [
            (FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[0]], Teffs[Teff_lohi[0]]),
            (FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[0]], Teffs[Teff_lohi[1]]),
            (FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[1]], Teffs[Teff_lohi[0]]),
            (FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[1]], Teffs[Teff_lohi[1]]),
            (FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[0]], Teffs[Teff_lohi[0]]),
            (FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[0]], Teffs[Teff_lohi[1]]),
            (FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[1]], Teffs[Teff_lohi[0]]),
            (FeHs[FeH_lohi[0]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[1]], Teffs[Teff_lohi[1]]),    
            (FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[0]], Teffs[Teff_lohi[0]]),
            (FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[0]], Teffs[Teff_lohi[1]]),
            (FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[1]], Teffs[Teff_lohi[0]]),
            (FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[0]], loggs[logg_lohi[1]], Teffs[Teff_lohi[1]]),
            (FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[0]], Teffs[Teff_lohi[0]]),
            (FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[0]], Teffs[Teff_lohi[1]]),
            (FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[1]], Teffs[Teff_lohi[0]]),
            (FeHs[FeH_lohi[1]], alphaMs[alphaM_lohi[1]], loggs[logg_lohi[1]], Teffs[Teff_lohi[1]])
        ]
        
        
        # load all spectra files from FITS file matching coorddinates in hypercube
        # change unnits to SI (J*s^-1*m^-2*m) from CGS (erg*s^-1*cm^-2*cm^-1)
        spectra = []
        for d, f in files:
            with fits.open(tempdir/files[0][1]) as hdul:
                spectra.append(hdul[0].data * 1e-1)
                
        # interpolate successively along each dimension
        # FeH interpolation:
        spectra2 = interpolate_dims(spectra, 8, (FeHs[FeH_lohi[0]], FeHs[FeH_lohi[1]]), FeH)
                
        # alphaM interpolation:
        spectra3 = interpolate_dims(spectra2, 4, (alphaMs[alphaM_lohi[0]], alphaMs[alphaM_lohi[1]]), alphaM)
        
        # logg interpolation:
        spectra4 = interpolate_dims(spectra3, 2, (loggs[logg_lohi[0]], loggs[logg_lohi[1]]), logg)
        
        # Teff interpolation:
        spectra5 = interpolate_dims(spectra4, 1, (Teffs[Teff_lohi[0]], Teffs[Teff_lohi[1]]), Teff)
        
    
        
        # Load wavelength points
        wavelength_grid = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
        
        # wavelength grid
        # change unit to SI (m) from Angstrom
        with fits.open(tempdir/wavelength_grid) as hdul:
            wavelengths = hdul[0].data*1e-10

    return wavelengths, spectra5[0]

