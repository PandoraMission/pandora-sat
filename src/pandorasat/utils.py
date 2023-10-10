# Standard library
import os

# Third-party
import astropy.units as u
import numpy as np
from astropy.constants import c, h
from astropy.io import fits
from astropy.time import Time
from astroquery import log as asqlog

from . import PACKAGEDIR, __version__

phoenixpath = f"{PACKAGEDIR}/data/phoenix"
os.environ["PYSYN_CDBS"] = phoenixpath

asqlog.setLevel("ERROR")


def photon_energy(wavelength):
    """ Converts photon wavelength to energy. """
    return ((h * c) / wavelength) * 1 / u.photon


def get_flatfield(stddev=0.005, seed=777):
    np.random.seed(seed)
    """ This generates and writes a dummy flatfield file. """
    for detector in ["VISDA", "NIRDA"]:
        hdr = fits.Header()
        hdr["AUTHOR"] = "Christina Hedges"
        hdr["VERSION"] = __version__
        hdr["DATE"] = Time.now().strftime("%d-%m-%Y")
        hdr["STDDEV"] = stddev
        hdu0 = fits.PrimaryHDU(header=hdr)
        hdulist = fits.HDUList(
            [
                hdu0,
                fits.CompImageHDU(
                    data=np.random.normal(1, stddev, (2048, 2048)), name="FLAT"
                ),
            ]
        )
        hdulist.writeto(
            f"{PACKAGEDIR}/data/flatfield_{detector}_{Time.now().strftime('%Y-%m-%d')}.fits",
            overwrite=True,
            checksum=True,
        )
    return


def load_vega():
    wavelength, spectrum = np.loadtxt(f"{PACKAGEDIR}/data/vega.dat").T
    wavelength *= u.angstrom
    spectrum *= u.erg / u.cm**2 / u.s / u.angstrom
    return wavelength, spectrum
