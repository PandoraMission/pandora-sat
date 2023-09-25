__version__ = "0.4.1"
# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
PANDORASTYLE = "{}/data/pandora.mplstyle".format(PACKAGEDIR)

# Standard library
import logging  # noqa: E402
import shutil  # noqa: E402
from glob import glob  # noqa: E402

# Third-party
from astropy.utils.data import download_file  # noqa: E402

from .utils import get_flatfield  # noqa: E402

logging.basicConfig()
logger = logging.getLogger("pandorasat")

from .pandorasat import PandoraSat  # noqa