"""Holds metadata and methods on Pandora NIRDA"""

# Standard library
from glob import glob
from dataclasses import dataclass

# Third-party
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io import fits

from . import PACKAGEDIR
from .optics import Optics
from .utils import photon_energy, load_vega


@dataclass
class NIRDetector:
    """Holds information on the Pandora IR detector

    Attributes
    ----------

    name: str
        Name of the detector. This will determine which files are loaded. This
        will be `"nirda"` for this detector
    pixel_scale: float
        The pixel scale of the detector in arcseconds/pixel
    pixel_size: float
        The pixel size in microns/mm
    """
    name: str
    pixel_scale: float
    pixel_size: float

    def _setup(self):
        self.shape = (2048, 512)
        """Some detector specific functions to run on initialization"""
        self.flat = fits.open(
            np.sort(
                np.atleast_1d(glob(f"{PACKAGEDIR}/data/flatfield_NIRDA*.fits"))
            )[-1]
        )[1].data

        # ROW COLUMN JUST LIKE PYTHON
        self.subarray_size = (400, 80)
        self.subarray_center = (self.wcs.wcs.crpix[1], self.wcs.wcs.crpix[0])
        crpix = self.wcs.wcs.crpix
        self.subarray_corner = (
            crpix[1] - self.subarray_size[0] / 2,
            crpix[0] - self.subarray_size[1] / 2,
        )
        # COLUMN, ROW
        self.subarray_row, self.subarray_column = np.meshgrid(
            crpix[1]
            + np.arange(self.subarray_size[0])
            - self.subarray_size[0] / 2,
            crpix[0]
            + np.arange(self.subarray_size[1])
            - self.subarray_size[1] / 2,
            indexing="ij",
        )
        self.trace_range = [-200, 100]

    @property
    def naxis1(self):
        """WCS's are COLUMN major, so naxis1 is the number of columns"""
        return self.shape[1] * u.pixel

    @property
    def naxis2(self):
        """WCS's are COLUMN major, so naxis2 is the number of rows"""
        return self.shape[0] * u.pixel

    @property
    def _dispersion_df(self):
        return pd.read_csv(f"{PACKAGEDIR}/data/pixel_vs_wavelength.csv")

    @property
    def pixel_read_time(self):
        return 1e-5 * u.second / u.pixel

    @property
    def frame_time(self):
        return np.product(self.subarray_size) * u.pixel * self.pixel_read_time

    @property
    def dark(self):
        return 1 * u.electron / u.second

    @property
    def read_noise(self):
        raise ValueError("Not Set")

    @property
    def saturation_limit(self):
        raise ValueError("Not Set")

    @property
    def non_linearity(self):
        raise ValueError("Not Set")

    def throughput(self, wavelength):
        return wavelength.value**0 * 0.61

    def apply_gain(self, values: u.Quantity):
        """Applies a single gain value"""
        return values * 0.5 * u.electron / u.DN

    def qe(self, wavelength):
        """
        Calculate the quantum efficiency of the detector.

        Parameters:
            wavelength (npt.NDArray): Wavelength in microns as `astropy.unit`

        Returns:
            qe (npt.NDArray): Array of the quantum efficiency of the detector
        """
        pass

    def sensitivity(self, wavelength):
        sed = 1 * u.erg / u.s / u.cm**2 / u.angstrom
        E = photon_energy(wavelength)
        telescope_area = np.pi * (Optics.mirror_diameter / 2) ** 2
        photon_flux_density = (
            telescope_area * sed * self.throughput(wavelength) / E
        ).to(u.photon / u.second / u.angstrom) * self.qe(wavelength)
        sensitivity = photon_flux_density / sed
        return sensitivity

    @property
    def midpoint(self):
        """Mid point of the sensitivity function"""
        w = np.arange(0.1, 3, 0.005) * u.micron
        return np.average(w, weights=self.sensitivity(w))

    def _estimate_zeropoint(self):
        """Use Vega SED to estimate the zeropoint of the detector"""
        wavelength, spectrum = load_vega()
        sens = self.sensitivity(wavelength)
        zeropoint = np.trapz(spectrum * sens, wavelength) / np.trapz(
            sens, wavelength
        )
        return zeropoint

    def mag_from_flux(self, flux):
        return -2.5 * np.log10(flux / self.zeropoint)

    def flux_from_mag(self, mag):
        return self.zeropoint * 10 ** (-mag / 2.5)
