"""Holds metadata and methods on Pandora"""

from .hardware import Hardware
from .irdetector import NIRDetector
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
        return "Pandora Observatory"

    def _repr_html_(self):
        return "Pandora Observatory"
