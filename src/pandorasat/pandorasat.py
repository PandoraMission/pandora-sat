"""Holds metadata and methods on Pandora"""

# Third-party
import numpy as np

from .irdetector import NIRDetector
from .hardware import Hardware
from .orbit import Orbit
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

    def __init__(self):
        self.Orbit = Orbit()
        self.Hardware = Hardware()
        self.NIRDA = NIRDetector()
        self.VISDA = VisibleDetector()

    def __repr__(self):
        return f"Pandora Observatory (RA: {np.round(self.ra, 3)}, Dec: {np.round(self.dec, 3)}, theta: {np.round(self.theta, 3)})"

    def _repr_html_(self):
        return f"Pandora Observatory (RA: {self.ra._repr_latex_()},  Dec:{self.dec._repr_latex_()}, theta: {self.theta._repr_latex_()})"
