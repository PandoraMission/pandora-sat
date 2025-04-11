"""Holds metadata and methods on Pandora VISDA"""

# Standard library
from dataclasses import dataclass

# Third-party
import astropy.units as u
import numpy as np
import pandas as pd

from . import PACKAGEDIR
from .hardware import Hardware
from .mixins import DetectorMixins
from .utils import photon_energy
from .wcs import get_wcs


@dataclass
class VisibleDetector(DetectorMixins):
    """
    Holds information on the Pandora Visible Detector
    """

    def __post_init__(self):
        """Some detector specific functions to run on initialization"""
        # self.flat = fits.open(
        #     np.sort(
        #         np.atleast_1d(glob(f"{PACKAGEDIR}/data/flatfield_VISDA*.fits"))
        #     )[-1]
        # )[1].data
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
            r = (self.fieldstop_radius / self.pixel_size).to(u.pix).value
            self.fieldstop = ~(np.hypot(R, C) > r)
            self._add_trace_params("visda")

        self.zeropoint = self.estimate_zeropoint()

    def __repr__(self):
        return "VisibleDetector"

    @property
    def name(self):
        return "VISDA"

    @property
    def shape(self):
        """Shape of the detector in pixels"""
        return (2048, 2048)

    @property
    def pixel_scale(self):
        """Pixel scale of the detector"""
        return 0.78 * u.arcsec / u.pixel

    @property
    def pixel_size(self):
        """Size of a pixel"""
        return 6.5 * u.um / u.pixel

    @property
    def bits_per_pixel(self):
        """Number of bits per pixel"""
        return 32 * u.bit / u.pixel

    @property
    def naxis1(self):
        """WCS's are COLUMN major, so naxis1 is the number of columns"""
        return self.shape[1] * u.pixel

    @property
    def naxis2(self):
        """WCS's are COLUMN major, so naxis2 is the number of rows"""
        return self.shape[0] * u.pixel

    @property
    def dark_rate(self):
        """Dark Noise"""
        return 1 * u.electron / u.second / u.pixel

    @property
    def read_noise(self):
        """Read Noise"""
        return 1.5 * u.electron / u.pixel

    @property
    def background_rate(self):
        """Detector background rate"""
        return 2 * u.electron / u.second / u.pixel

    @property
    def bias(self):
        """Detector bias"""
        return 100 * u.DN

    @property
    def integration_time(self):
        "Integration time"
        return 0.2 * u.second

    @property
    def fieldstop_radius(self):
        "Radius of the fieldstop"
        return 6.5 * u.mm

    def throughput(self, wavelength: u.Quantity):
        """Optical throughput at the specified wavelength(s)"""
        df = pd.read_csv(f"{PACKAGEDIR}/data/visible_optical_throughput.csv")
        throughput = np.interp(
            wavelength.to(u.nm).value, *np.asarray(df.values).T
        )
        throughput[wavelength.to(u.nm).value < 380] *= 0
        return throughput

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
        df = pd.read_csv(f"{PACKAGEDIR}/data/Pandora_Visible_QE.csv")
        wav, transmission = (
            np.asarray(df.Wavelength) * u.angstrom,
            np.asarray(df.Transmission),
        )
        return (
            np.interp(wavelength, wav, transmission, left=0, right=0)
            * u.electron
            / u.photon
        )

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

    def apply_gain(self, values: u.Quantity):
        """Applies a piecewise gain function"""
        if not isinstance(values, u.Quantity):
            raise ValueError("Must pass a quantity.")
        x = np.atleast_1d(values)
        gain = np.asarray([0.52, 0.6, 0.61, 0.67]) * u.electron / u.DN

        if values.unit == u.electron:
            masks = np.asarray(
                [
                    (x >= 0 * u.electron) & (x < 520 * u.electron),
                    (x >= 520 * u.electron) & (x < 3000 * u.electron),
                    (x >= 3000 * u.electron) & (x < 17080 * u.electron),
                    (x >= 17080 * u.electron),
                ]
            )
            if values.ndim <= 1:
                gain = gain[:, None]
            if values.ndim == 2:
                gain = gain[:, None, None]
            result = u.Quantity(
                (masks * x[None, :] / gain).sum(axis=0), dtype=int, unit=u.DN
            )
            if values.ndim == 0:
                return result[0]
            return result

        elif values.unit == u.DN:
            masks = np.asarray(
                [
                    (x >= 0 * u.DN) & (x < 1e3 * u.DN),
                    (x >= 1e3 * u.DN) & (x < 5e3 * u.DN),
                    (x >= 5e3 * u.DN) & (x < 2.8e4 * u.DN),
                    (x >= 2.8e4 * u.DN),
                ]
            )
            if values.ndim <= 1:
                gain = gain[:, None]
            if values.ndim == 2:
                gain = gain[:, None, None]
            result = u.Quantity(
                (masks * x[None, :] * gain).sum(axis=0),
                dtype=int,
                unit=u.electron,
            )
            if values.ndim == 0:
                return result[0]
            return result

    def get_wcs(
        self,
        ra,
        dec,
        theta=u.Quantity(0, unit="degree"),
        crpix1=1024,
        crpix2=1024,
        order=3,
    ):
        """Returns an astropy.wcs.WCS object"""
        return get_wcs(
            self,
            target_ra=ra,
            target_dec=dec,
            theta=theta,
            crpix1=crpix1,
            crpix2=crpix2,
            order=order,
            distortion_file=f"{PACKAGEDIR}/data/fov_distortion.csv",
        )

    @property
    def info(self):
        zp = self.zeropoint
        return pd.DataFrame(
            {
                "Detector Size": f"({self.naxis1.value.astype(int)}, {self.naxis2.value.astype(int)})",
                "Pixel Scale": f"{self.pixel_scale.value} {self.pixel_scale.unit.to_string('latex')}",
                "Pixel Size": f"{self.pixel_size.value} {self.pixel_size.unit.to_string('latex')}",
                "Read Noise": f"{self.read_noise.value} {self.read_noise.unit.to_string('latex')}",
                "Dark Noise": f"{self.dark_rate.value} {self.dark_rate.unit.to_string('latex')}",
                "Bias": f"{self.bias.value} {self.bias.unit.to_string('latex')}",
                "Wavelength Midpoint": f"{self.midpoint.value:.2f} {self.midpoint.unit.to_string('latex')}",
                "Integration Time": f"{self.integration_time.value} {self.integration_time.unit.to_string('latex')}",
                "Zeropoint": f"{zp.value:.3e}"
                + "$\\mathrm{\\frac{erg}{A\\,s\\,cm^{2}}}$",
            },
            index=[0],
        ).T.rename({0: "VISDA"}, axis="columns")
