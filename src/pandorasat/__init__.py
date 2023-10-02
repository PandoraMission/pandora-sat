__version__ = "0.4.1"
# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
PANDORASTYLE = "{}/data/pandora.mplstyle".format(PACKAGEDIR)

# Standard library
import logging  # noqa: E402


logging.basicConfig()
logger = logging.getLogger("pandorasat")

from .pandorasat import PandoraSat  # noqa