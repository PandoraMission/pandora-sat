"""Holds basic metadata on the optics of Pandora"""

# Standard library
from dataclasses import dataclass

# Third-party
import astropy.units as u
import numpy as np


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
        """Diameter of Pandora's mirror. This should be the effective diameter of the primary, removing the secondary diameter."""
        # return 0.43 * u.m
        A_primary = np.pi * (self.primary_mirror_effective_diameter / 2) ** 2
        A_secondary = (
            np.pi * (self.secondary_mirror_physical_diameter / 2) ** 2
        )
        r_primary = ((A_primary - A_secondary) / np.pi) ** 0.5
        return r_primary * 2

    @property
    def primary_mirror_effective_diameter(self):
        """Diameter of Pandora's primary mirror that will reflect light."""
        return (43.5 * u.cm).to(u.m)

    @property
    def secondary_mirror_physical_diameter(self):
        """Diameter of Pandora's secondary mirror that will block light."""
        return (86 * u.mm).to(u.m)
