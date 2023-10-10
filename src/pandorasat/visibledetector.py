"""Holds metadata and methods on Pandora VISDA"""
# Standard library
from glob import glob
from dataclasses import dataclass

# Third-party
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io import fits

from . import PACKAGEDIR
from .hardware import Optics
from .utils import photon_energy, load_vega


@dataclass
class VisibleDetector:
    """Holds information on the Pandora Visible Detector

    Attributes
    ----------
    name: str
        Name of the detector. This will determine which files are loaded. This
        will be `"visda"` for this detector
    pixel_scale: float
        The pixel scale of the detector in arcseconds/pixel
    pixel_size: float
        The pixel size in microns/mm
    """

    def __post_init__(self):
        """Some detector specific functions to run on initialization"""
        self.shape = (2048, 2048)

        self.flat = fits.open(
            np.sort(
                np.atleast_1d(glob(f"{PACKAGEDIR}/data/flatfield_VISDA*.fits"))
            )[-1]
        )[1].data
        if hasattr(self, "fieldstop_radius"):
            C, R = (
                np.mgrid[
                    : self.shape[0],
                    : self.shape[1],
                ]
                - np.hstack(
                    [
                        self.shape[0],
                        self.shape[1],
                    ]
                )[:, None, None]
                / 2
            )
            r = (self.fieldstop_radius / self.pixel_scale).to(u.pix).value
            self.fieldstop = ~((np.abs(C) >= r) | (np.abs(R) >= r))

    @property
    def pixel_scale(self):
        """Pixel scale of the detector"""
        return 0.78 * u.arcsec / u.pixel

    @property
    def pixel_size(self):
        """Size of a pixel"""
        return 6.5 * u.um / u.pixel

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
        return pd.read_csv(f"{PACKAGEDIR}/data/pixel_vs_wavelength_vis.csv")

    @property
    def dark(self):
        return 1 * u.electron / u.second

    @property
    def read_noise(self):
        return 1.5 * u.electron

    @property
    def bias(self):
        return 100 * u.electron

    @property
    def integration_time(self):
        return 0.2 * u.second

#    @property
#    def fieldstop_radius(self):
#        return 0.3 * u.deg

    def throughput(self, wavelength):
        return wavelength.value**0 * 0.714

    def qe(self, wavelength):
        """
        Calculate the quantum efficiency of the detector.

        Parameters
        ----------
        wavelength : npt.NDArray
            Wavelength in microns as `astropy.unit`

        Returns
        -------
        qe : npt.NDArray
            Array of the quantum efficiency of the detector
        """
        raise NotImplementedError

    def sensitivity(self, wavelength):
        """
        Calulate the sensitivity of the detector.

        Parameters
        ----------
        wavelength : npt.NDArray
            Wavelength in microns as `astropy.unit`

        Returns
        -------
        sensitivity : npt.NDArray
            Array of the sensitivity of the detector
        """
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

    def estimate_zeropoint(self):
        """Use Vega SED to estimate the zeropoint of the detector"""
        wavelength, spectrum = load_vega()
        sens = self.sensitivity(wavelength)
        zeropoint = np.trapz(spectrum * sens, wavelength) / np.trapz(
            sens, wavelength
        )
        return zeropoint

    def mag_from_flux(self, flux):
        """Convert flux to magnitude based on the zeropoint of the detector"""
        return -2.5 * np.log10(flux / self.zeropoint)

    def flux_from_mag(self, mag):
        """Convert magnitude to flux based on the zeropoint of the detector"""
        return self.zeropoint * 10 ** (-mag / 2.5)
