__version__ = "0.1.2"
# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

from .pandorasat import PandoraSat  # noqa
from .psf import PSF  # noqa
from .targets import Target  # noqa
