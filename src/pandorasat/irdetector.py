"""Holds metadata and methods on Pandora NIRDA"""

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
class NIRDetector(DetectorMixins):
    """
    Holds information on the Pandora IR detector
    """

    def __post_init__(self):
        """Some detector specific functions to run on initialization"""
        self.reference = pr.NIRDAReference()

    def __repr__(self):
        return "NIRDetector"

    @property
    def name(self):
        return "NIRDA"

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
    def bits_per_pixel(self):
        """Number of bits per pixel"""
        return 16 * u.bit / u.pixel

    @property
    def naxis1(self):
        """WCS's are COLUMN major, so naxis1 is the number of columns"""
        return self.shape[1] * u.pixel

    @property
    def naxis2(self):
        """WCS's are COLUMN major, so naxis2 is the number of rows"""
        return self.shape[0] * u.pixel

    @property
    def pixel_read_time(self):
        """Pixel read time"""
        return 1e-5 * u.second / u.pixel

    @property
    def subarray_size(self):
        """Size of standard subarray"""
        return (400, 80)

    def frame_time(self, array_size=None):
        """Time to read out one frame of the subarray"""
        if array_size is None:
            array_size = self.subarray_size
        return np.prod(array_size) * u.pixel * self.pixel_read_time

    @property
    def zodiacal_background_rate(self):
        "Zodiacal light background rate"
        return 4 * u.electron / u.second / u.pixel

    @property
    def stray_light_rate(self):
        "Stray light rate"
        return 2 * u.electron / u.second / u.pixel

    @property
    def thermal_background_rate(self):
        "NIRDA thermal background rate"
        return 10 * u.electron / u.second / u.pixel

    @property
    def correlated_double_sampling_readnoise(self):
        """This is the read noise obtained when differencing two images."""
        return self.readnoise * np.sqrt(2)

    @property
    def bias_uncertainty(self):
        "Uncertainty in NIRDA detector bias. Every integration has a different bias."
        return (185 * 2) * u.electron

    @property
    def saturation_limit(self):
        "NIRDA saturation limit. Bias contributes to saturation."
        return 80000 * u.electron

    @property
    def midpoint(self):
        """Mid point of the sensitivity function"""
        w = np.arange(0.1, 3, 0.005) * u.micron
        return np.average(w, weights=self.sensitivity(w))

    @property
    def info(self):
        return pd.DataFrame(
            {
                "Detector Size": f"{self.shape}",
                "Subarray Size": f"{self.subarray_size}",
                "Pixel Scale": f"{self.pixel_scale.value} {self.pixel_scale.unit.to_string('latex')}",
                "Pixel Size": f"{self.pixel_size.value} {self.pixel_size.unit.to_string('latex')}",
                "Dark Noise": f"{self.dark.value} {self.dark.unit.to_string('latex')}",
                "Read Noise": f"{self.readnoise.value} {self.readnoise.unit.to_string('latex')}",
                "Wavelength Midpoint": f"{self.midpoint.value:.2f} {self.midpoint.unit.to_string('latex')}",
                "Pixel Read Time": f"{self.pixel_read_time.value:.1e} {self.pixel_read_time.unit.to_string('latex')}",
                "Zeropoint": f"{self.zeropoint.value:.3e}"
                + "$\\mathrm{\\frac{erg}{A\\,s\\,cm^{2}}}$",
                "R @ 1.3$\\mu m$": 65,
            },
            index=[0],
        ).T.rename({0: "NIRDA"}, axis="columns")

    def plot_sensitivity(self, ax=None):
        """Plot the sensitivity of the detector as a function of wavelength"""
        wavelength = np.linspace(0.6, 2, 1000) * u.micron
        pixel = self.reference.get_pixel_position(wavelength=wavelength)
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
            ax.spines[["right", "top"]].set_visible(True)
            ax_p = ax.twiny()
            ax_p.set(xticks=ax.get_xticks(), xlim=ax.get_xlim())
            ax_p.set_xlabel(xlabel="$\delta$ Pixel Position", color="grey")
            ax_p.set_xticklabels(
                labels=list(
                    np.interp(
                        ax.get_xticks(),
                        wavelength.value,
                        pixel.value,
                    ).astype(int)
                ),
                rotation=45,
                color="grey",
            )
        return ax
