"""Holds metadata and methods on Pandora"""

# Third-party
import astropy.units as u
import numpy as np
import pandas as pd

from .irdetector import NIRDetector
from .optics import Optics
from .orbit import Orbit
from .utils import get_jitter
from .visibledetector import VisibleDetector


# @dataclass
class PandoraSat(object):
    """Holds information and methods for the full Pandora system.

    Args:
        NIRDA (IRDetector): Class of the NIRDA properties
        VISDA (IRDetector): Class of the VISDA properties
        Optics (IRDetector): Class of the Optics properties
        Orbit (IRDetector): Class of the Orbit properties
    """

    def __init__(
        self,
        rowjitter_1sigma: u.Quantity = 0.2 * u.pixel,
        coljitter_1sigma: u.Quantity = 0.2 * u.pixel,
        thetajitter_1sigma: u.Quantity = 0.0005 * u.deg,
        jitter_timescale: u.Quantity = 60 * u.second,
    ):
        self.Orbit = Orbit()
        self.Optics = Optics()
        self.NIRDA = NIRDetector(
            "NIR",
            1.19 * u.arcsec / u.pixel,
            18.0 * u.um / u.pixel,
        )
        self.VISDA = VisibleDetector(
            "Visible",
            0.78 * u.arcsec / u.pixel,
            6.5 * u.um / u.pixel,
        )
        self.rowjitter_1sigma = rowjitter_1sigma
        self.coljitter_1sigma = coljitter_1sigma
        self.thetajitter_1sigma = thetajitter_1sigma
        self.jitter_timescale = jitter_timescale
        self._get_jitter()

        nints = (duration.to(u.s) / self.VISDA.integration_time).value
        self.VISDA.time = (
            self.obstime.jd
            + np.arange(0, nints) / nints * duration.to(u.day).value
        )
        self.VISDA.rowj, self.VISDA.colj, self.VISDA.thetaj = (  # noqa
            np.interp(self.VISDA.time, self.jitter.time, self.jitter.rowj),
            np.interp(self.VISDA.time, self.jitter.time, self.jitter.colj),
            np.interp(self.VISDA.time, self.jitter.time, self.jitter.thetaj),
        )

    def __repr__(self):
        return f"Pandora Observatory (RA: {np.round(self.ra, 3)}, Dec: {np.round(self.dec, 3)}, theta: {np.round(self.theta, 3)})"

    def _repr_html_(self):
        return f"Pandora Observatory (RA: {self.ra._repr_latex_()},  Dec:{self.dec._repr_latex_()}, theta: {self.theta._repr_latex_()})"

    def _get_jitter(self):
        self.jitter = pd.DataFrame(
            np.asarray(
                get_jitter(
                    self.rowjitter_1sigma.value,
                    self.coljitter_1sigma.value,
                    self.thetajitter_1sigma.value,
                    correlation_time=self.jitter_timescale,
                    nframes=int(
                        (
                            5
                            * self.duration.to(u.second)
                            / self.jitter_timescale.to(u.second)
                        ).value
                    ),
                    frame_time=self.jitter_timescale / 5,
                )
            ).T,
            columns=["time", "rowj", "colj", "thetaj"],
        )
        self.jitter.time = (
            self.jitter.time * (u.second).to(u.day)
        ) + self.obstime.jd
