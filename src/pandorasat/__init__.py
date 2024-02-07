__version__ = "0.5.0"
# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# Standard library
import logging  # noqa: E402
from glob import glob  # noqa: E402

pandorastyle = glob(f"{PACKAGEDIR}/data/pandora.mplstyle")


from .utils import get_flatfield  # noqa: E402

logging.basicConfig()
logger = logging.getLogger("pandorasat")

from .pandorasat import PandoraSat  # noqa

flatnames = glob(f"{PACKAGEDIR}/data/flatfield_*.fits")
if len(flatnames) == 0:
    # Make a bogus flatfield
    logger.warning("No flatfield file found. Generating a random one for you.")
    get_flatfield()
    logger.warning(
        f"Generated flatfield in {PACKAGEDIR}/data/pandora_nir_20220506.fits."
    )

from .irdetector import NIRDetector  # noqa: E402, F401
from .visibledetector import VisibleDetector  # noqa: E402, F401
