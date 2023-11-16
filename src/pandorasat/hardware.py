"""Holds basic metadata on the optics of Pandora"""

# Standard library
from dataclasses import dataclass

# Third-party
import astropy.units as u


@dataclass
class Hardware:
    """Holds basic metadata on the optics of Pandora

    Args:
        mirror_diameter (float): Diameter of the Pandora mirror
    """

    def __repr__(self):
        return "Pandora Optics"

    @property
    def mirror_diameter(self):
        """Diameter of Pandora's mirror"""
        return 0.43 * u.m
