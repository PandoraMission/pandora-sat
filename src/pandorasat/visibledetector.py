"""Holds metadata and methods on Pandora VISDA"""

# Standard library
from dataclasses import dataclass

# Third-party
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandoraref as pr

from . import PANDORASTYLE
from .detectormixins import DetectorMixins


@dataclass
class VisibleDetector(DetectorMixins):
    """
    Holds information on the Pandora Visible Detector
    """

    def __post_init__(self):
        """Some detector specific functions to run on initialization"""
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
        self.reference = pr.VISDAReference()

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
    def background_rate(self):
        """Detector background rate"""
        return 2 * u.electron / u.second / u.pixel

    @property
    def integration_time(self):
        "Integration time"
        return 0.2 * u.second

    @property
    def fieldstop_radius(self):
        "Radius of the fieldstop"
        return 6.5 * u.mm

    @property
    def midpoint(self):
        """Mid point of the sensitivity function"""
        w = np.arange(0.1, 3, 0.005) * u.micron
        return np.average(w, weights=self.sensitivity(w))

    # def apply_gain(self, values: u.Quantity):
    #     """Applies a piecewise gain function"""
    #     if not isinstance(values, u.Quantity):
    #         raise ValueError("Must pass a quantity.")
    #     x = np.atleast_1d(values)
    #     gain = np.asarray([0.52, 0.6, 0.61, 0.67]) * u.electron / u.DN

    #     if values.unit == u.electron:
    #         masks = np.asarray(
    #             [
    #                 (x >= 0 * u.electron) & (x < 520 * u.electron),
    #                 (x >= 520 * u.electron) & (x < 3000 * u.electron),
    #                 (x >= 3000 * u.electron) & (x < 17080 * u.electron),
    #                 (x >= 17080 * u.electron),
    #             ]
    #         )
    #         if values.ndim <= 1:
    #             gain = gain[:, None]
    #         if values.ndim == 2:
    #             gain = gain[:, None, None]
    #         result = u.Quantity(
    #             (masks * x[None, :] / gain).sum(axis=0), dtype=int, unit=u.DN
    #         )
    #         if values.ndim == 0:
    #             return result[0]
    #         return result

    #     elif values.unit == u.DN:
    #         masks = np.asarray(
    #             [
    #                 (x >= 0 * u.DN) & (x < 1e3 * u.DN),
    #                 (x >= 1e3 * u.DN) & (x < 5e3 * u.DN),
    #                 (x >= 5e3 * u.DN) & (x < 2.8e4 * u.DN),
    #                 (x >= 2.8e4 * u.DN),
    #             ]
    #         )
    #         if values.ndim <= 1:
    #             gain = gain[:, None]
    #         if values.ndim == 2:
    #             gain = gain[:, None, None]
    #         result = u.Quantity(
    #             (masks * x[None, :] * gain).sum(axis=0),
    #             dtype=int,
    #             unit=u.electron,
    #         )
    #         if values.ndim == 0:
    #             return result[0]
    #         return result

    @property
    def info(self):
        zp = self.zeropoint
        return pd.DataFrame(
            {
                "Detector Size": f"({self.naxis1.value.astype(int)}, {self.naxis2.value.astype(int)})",
                "Pixel Scale": f"{self.pixel_scale.value} {self.pixel_scale.unit.to_string('latex')}",
                "Pixel Size": f"{self.pixel_size.value} {self.pixel_size.unit.to_string('latex')}",
                "Dark Noise": f"{self.dark.value} {self.dark.unit.to_string('latex')}",
                "Bias": f"{self.bias.value.mean()} {self.bias.unit.to_string('latex')}",
                "Read Noise": f"{self.readnoise.value} {self.readnoise.unit.to_string('latex')}",
                "Wavelength Midpoint": f"{self.midpoint.value:.2f} {self.midpoint.unit.to_string('latex')}",
                "Integration Time": f"{self.integration_time.value} {self.integration_time.unit.to_string('latex')}",
                "Zeropoint": f"{zp.value:.3e}"
                + "$\\mathrm{\\frac{erg}{A\\,s\\,cm^{2}}}$",
            },
            index=[0],
        ).T.rename({0: "VISDA"}, axis="columns")

    def plot_sensitivity(self, ax=None):
        """Plot the sensitivity of the detector as a function of wavelength"""
        wavelength = np.linspace(0.1, 1, 1000) * u.micron
        sens = self.reference.get_sensitivity(wavelength=wavelength)
        if ax is None:
            _, ax = plt.subplots()
        with plt.style.context(PANDORASTYLE):
            ax.plot(
                wavelength,
                sens,
                c="k",
            )
            ax.set(
                xticks=np.linspace(*ax.get_xlim(), 9),
                xlabel=f"Wavelength [{wavelength.unit.to_string('latex')}]",
                ylabel=f"Sensitivity [{sens.unit.to_string('latex')}]",
                title=self.name.upper(),
            )
        return ax
