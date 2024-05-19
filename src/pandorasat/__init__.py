__version__ = "0.5.20"
# Standard library
import logging  # noqa: E402
import os  # noqa
import time  # noqa: E402
from glob import glob
from threading import Event, Thread  # noqa: E402

# Third-party
from rich.console import Console  # noqa: E402
from rich.logging import RichHandler  # noqa: E402

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
TESTDIR = "/".join(PACKAGEDIR.split("/")[:-2]) + "/tests/"
PANDORASTYLE = glob(f"{PACKAGEDIR}/data/pandora.mplstyle")


def get_logger(name="pandoralog"):
    """Configure and return a logger with RichHandler."""
    return PandoraLogger(name)


# Custom Logger with Rich
class PandoraLogger(logging.Logger):
    def __init__(self, name, level=logging.INFO):
        super().__init__(name, level)
        console = Console()
        self.handler = RichHandler(
            show_time=False, show_level=False, show_path=False, console=console
        )
        self.handler.setFormatter(
            logging.Formatter(
                "%(asctime)s %(levelname)s: %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        self.addHandler(self.handler)
        self.spinner_thread = None
        self.spinner_event = None

    def start_spinner(self, message="Processing..."):
        if self.spinner_thread is None:
            self.spinner_event = Event()
            self.spinner_thread = Thread(target=self._spinner, args=(message,))
            self.spinner_thread.start()

    def stop_spinner(self):
        if self.spinner_thread is not None:
            self.spinner_event.set()
            self.spinner_thread.join()
            self.spinner_thread = None
            self.spinner_event = None

    def _spinner(self, message):
        with self.handler.console.status(
            "[bold green]" + message
        ) as status:  # noqa
            while not self.spinner_event.is_set():
                time.sleep(0.1)


logger = get_logger("pandorasat")


from .irdetector import NIRDetector  # noqa: E402, F401
from .mixins import DetectorMixins  # noqa: E402, F401
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
