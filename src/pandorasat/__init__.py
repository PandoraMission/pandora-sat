__version__ = "0.5.7"
# Standard library
import os  # noqa
from glob import glob

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# Standard library
import logging  # noqa: E402

pandorastyle = glob(f"{PACKAGEDIR}/data/pandora.mplstyle")

logging.basicConfig()
logger = logging.getLogger("pandorasat")

from .irdetector import NIRDetector  # noqa: E402, F401
from .pandorasat import PandoraSat  # noqa
from .visibledetector import VisibleDetector  # noqa: E402, F401

# from
# flatnames = glob(f"{PACKAGEDIR}/data/flatfield_*.fits")
# if len(flatnames) == 0:
#     # Make a bogus flatfield
#     logger.warning("No flatfield file found. Generating a random one for you.")
#     get_flatfield()
#     logger.warning(
#         f"Generated flatfield in {PACKAGEDIR}/data/pandora_nir_20220506.fits."
#     )
