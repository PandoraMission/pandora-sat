"""Holds metadata and methods on Pandora NIRDA"""

# Standard library
from glob import glob

# Third-party
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io import fits

from . import PACKAGEDIR
from .detector import Detector


class NIRDetector(Detector):
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