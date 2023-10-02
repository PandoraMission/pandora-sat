"""Holds basic metadata on Pandora orbit"""

# Standard library
from dataclasses import dataclass

# Third-party
import astropy.units as u


@dataclass
class Orbit:
    """Holds basic metadata on the orbit of Pandora"""

    period: float = 90 * u.minute

#   orbital velocity vector
#   altitude
#   RA/Dec (position and pointing)

#   equation of ellipse to interpolate?
#   query to JPL Horizons upon launch

    def __repr__(self):
        return "Pandora Orbit"
