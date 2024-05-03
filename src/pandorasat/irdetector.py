"""Holds metadata and methods on Pandora NIRDA"""

# Standard library
from dataclasses import dataclass

# Third-party
import astropy.units as u
import numpy as np
import pandas as pd

from . import PACKAGEDIR
from .hardware import Hardware
from .utils import load_vega, photon_energy


@dataclass
class NIRDetector:
    """
    Holds information on the Pandora IR detector
    """

    def __post_init__(self):
        """Some detector specific functions to run on initialization"""
        # self.flat = fits.open(
        #     np.sort(
        #         np.atleast_1d(glob(f"{PACKAGEDIR}/data/flatfield_NIRDA*.fits"))
        #     )[-1]
        # )[1].data

    def __repr__(self):
        return "NIRDetector"

    @property
    def shape(self):
        """Shape of the detector in pixels"""
        return (2048, 2048)

    @property
    def pixel_scale(self):
        """Pixel scale of the detector"""
        return 1.19 * u.arcsec / u.pixel

    @property
    def pixel_size(self):
        """Size of a pixel"""
        return 18.0 * u.um / u.pixel

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
    def subarray_size(self):
        return (400, 80)

    def frame_time(self, array_size=None):
        if array_size is None:
            array_size = self.subarray_size
        return np.product(array_size) * u.pixel * self.pixel_read_time

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

    def throughput(self, wavelength: u.Quantity):
        """Throughput at the specified wavelength(s)"""
        df = pd.read_csv(f"{PACKAGEDIR}/data/dichroic-transmission.csv")
        throughput = (
            np.interp(wavelength.to(u.nm).value, *np.asarray(df).T) / 100
        )
        return throughput**2 * 0.71

    def apply_gain(self, values: u.Quantity):
        """Applies a single gain value"""
        return values * 0.5 * u.electron / u.DN

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
        if not hasattr(wavelength, "unit"):
            raise ValueError("Pass a wavelength with units")
        wavelength = np.atleast_1d(wavelength)
        sw_coeffs = np.array([0.65830, -0.05668, 0.25580, -0.08350])
        sw_exponential = 100.0
        sw_wavecut_red = 1.69  # changed from 2.38 for Pandora
        sw_wavecut_blue = 0.75  # new for Pandora
        with np.errstate(invalid="ignore", over="ignore"):
            sw_qe = (
                sw_coeffs[0]
                + sw_coeffs[1] * wavelength.to(u.micron).value
                + sw_coeffs[2] * wavelength.to(u.micron).value ** 2
                + sw_coeffs[3] * wavelength.to(u.micron).value ** 3
            )

            sw_qe = np.where(
                wavelength.to(u.micron).value > sw_wavecut_red,
                sw_qe
                * np.exp(
                    (sw_wavecut_red - wavelength.to(u.micron).value)
                    * sw_exponential
                ),
                sw_qe,
            )

            sw_qe = np.where(
                wavelength.to(u.micron).value < sw_wavecut_blue,
                sw_qe
                * np.exp(
                    -(sw_wavecut_blue - wavelength.to(u.micron).value)
                    * (sw_exponential / 1.5)
                ),
                sw_qe,
            )
        sw_qe[sw_qe < 1e-5] = 0
        return sw_qe * u.electron / u.photon

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
        telescope_area = np.pi * (Hardware().mirror_diameter / 2) ** 2
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
        return -2.5 * np.log10(flux / self.estimate_zeropoint())

    def flux_from_mag(self, mag):
        """Convert magnitude to flux based on the zeropoint of the detector"""
        return self.estimate_zeropoint() * 10 ** (-mag / 2.5)
