# Standard library
import configparser  # noqa: E402
import logging  # noqa: E402
import os  # noqa
import time  # noqa: E402
from glob import glob
from threading import Event, Thread  # noqa: E402

# Third-party
import numpy as np  # noqa
import pandas as pd  # noqa
from appdirs import user_config_dir, user_data_dir  # noqa: E402
from rich.console import Console  # noqa: E402
from rich.logging import RichHandler  # noqa: E402

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
TESTDIR = "/".join(PACKAGEDIR.split("/")[:-2]) + "/tests/"
DOCSDIR = "/".join(PACKAGEDIR.split("/")[:-2]) + "/docs/"
PANDORASTYLE = glob(f"{PACKAGEDIR}/data/pandora.mplstyle")

# Standard library
from importlib.metadata import PackageNotFoundError, version  # noqa


def get_version():
    try:
        return version("pandorasat")
    except PackageNotFoundError:
        return "unknown"


__version__ = get_version()


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


def get_logger(name="pandorasat"):
    """Configure and return a logger with RichHandler."""
    return PandoraLogger(name)


logger = get_logger("pandorasat")


CONFIGDIR = user_config_dir("pandorasat")
os.makedirs(CONFIGDIR, exist_ok=True)
CONFIGPATH = os.path.join(CONFIGDIR, "config.ini")


def reset_config():
    config = configparser.ConfigParser()
    config["SETTINGS"] = {
        "log_level": "INFO",
        "data_dir": user_data_dir("pandorasat"),
    }

    with open(CONFIGPATH, "w") as configfile:
        config.write(configfile)


def load_config() -> configparser.ConfigParser:
    """
    Loads the configuration file, creating it with defaults if it doesn't exist.

    Returns
    -------
    configparser.ConfigParser
        The loaded configuration.
    """

    config = configparser.ConfigParser()

    if not os.path.exists(CONFIGPATH):
        # Create default configuration
        reset_config()
    config.read(CONFIGPATH)
    return config


def save_config(config: configparser.ConfigParser) -> None:
    """
    Saves the configuration to the file.

    Parameters
    ----------
    config : configparser.ConfigParser
        The configuration to save.
    app_name : str
        Name of the application.
    """
    with open(CONFIGPATH, "w") as configfile:
        config.write(configfile)


config = load_config()

for key in ["data_dir", "log_level"]:
    if key not in config["SETTINGS"]:
        logger.error(
            f"`{key}` missing from the `gaiaoffline` config file. Your configuration is being reset."
        )
        reset_config()
        config = load_config()

logger.setLevel(config["SETTINGS"]["log_level"])
CACHEDIR = config["SETTINGS"]["data_dir"]
os.makedirs(CACHEDIR, exist_ok=True)
PHOENIXPATH = f"{CACHEDIR}/data/phoenix/"
PHOENIXGRIDPATH = f"{CACHEDIR}/data/phoenix/grid/phoenix/phoenixm00/"


def display_config() -> pd.DataFrame:
    dfs = []
    for section in config.sections():
        df = pd.DataFrame(
            np.asarray(
                [(key, value) for key, value in dict(config[section]).items()]
            )
        )
        df["section"] = section
        df.columns = ["key", "value", "section"]
        df = df.set_index(["section", "key"])
        dfs.append(df)
    return pd.concat(dfs)


from .irdetector import NIRDetector  # noqa: E402, F401
from .mixins import DetectorMixins  # noqa: E402, F401
from .pandorasat import PandoraSat  # noqa
from .visibledetector import VisibleDetector  # noqa: E402, F401
