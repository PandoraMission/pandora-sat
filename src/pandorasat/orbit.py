"""Holds basic metadata on Pandora orbit"""

# Standard library
from dataclasses import dataclass

# Third-party
import astropy.units as u


@dataclass
class Orbit:
    """Holds basic metadata on the orbit of Pandora"""

    period: float = 96.54 * u.minute  # circular orbit @ 600 km

    def __repr__(self):
        return "Pandora Orbit"
