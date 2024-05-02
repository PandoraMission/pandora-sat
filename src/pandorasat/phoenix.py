# Standard library
import contextlib
import functools
import logging
import os
import warnings
from glob import glob
from pathlib import Path

# Third-party
import astropy.units as u
import numpy as np
import requests
from astroquery import log as asqlog

# from . import PACKAGEDIR

logging.basicConfig()
logger = logging.getLogger("pandorasat")

__all__ = [
    "download_phoenix_grid",
    "phoenixcontext",
    "build_phoenix",
    "get_phoenix_model",
]
CACHEDIR = f"{Path.home()}/.pandorasat"
PHOENIXPATH = f"{CACHEDIR}/data/phoenix/"
PHOENIXGRIDPATH = f"{CACHEDIR}/data/phoenix/grid/phoenix/phoenixm00/"


def download_file(file_url, file_path):
    # Download the file from `file_url` and save it locally under `file_path`
    with requests.get(file_url, stream=True) as r:
        r.raise_for_status()
        with open(file_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    # astropy_download_file(file_url, cache=True, show_progress=False, pkgname='pandorasat')


def download_phoenix_grid():
    # Standard library
    import shutil
    import tarfile

    # Third-party
    import numpy as np
    from tqdm import tqdm

    os.makedirs(CACHEDIR, exist_ok=True)
    if os.path.isdir(PHOENIXPATH):
        shutil.rmtree(PHOENIXPATH)
    os.makedirs(PHOENIXGRIDPATH, exist_ok=True)
    url = "https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix/phoenixm00/"
    page = requests.get(url).text
    suffix = '.fits">'
    filenames = np.asarray(
        [
            f"{i.split(suffix)[0]}.fits"
            for i in page.split('<li>&#x1f4c4; <a href="')[2:]
        ]
    )
    temperatures = np.asarray(
        [int(name.split("_")[1].split(".fits")[0]) for name in filenames]
    )
    filenames, temperatures = (
        filenames[np.argsort(temperatures)],
        temperatures[np.argsort(temperatures)],
    )
    filenames = filenames[temperatures < 10000]
    _ = [
        download_file(f"{url}{filename}", f"{PHOENIXGRIDPATH}/{filename}")
        for filename in tqdm(
            filenames,
            desc="Downloading PHOENIX Models",
            leave=True,
            position=0,
        )
    ]
    download_file(
        "http://ssb.stsci.edu/trds/tarfiles/synphot1.tar.gz",
        f"{PHOENIXPATH}synphot1.tar.gz",
    )
    with tarfile.open(f"{PHOENIXPATH}synphot1.tar.gz") as tar:
        tar.extractall(path=f"{PHOENIXPATH}")
    os.remove(f"{PHOENIXPATH}synphot1.tar.gz")
    fnames = glob(f"{PHOENIXPATH}grp/redcat/trds/*")
    _ = [shutil.move(fname, f"{PHOENIXPATH}") for fname in fnames]
    os.removedirs(f"{PHOENIXPATH}grp/redcat/trds/")
    download_file(
        "https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix/catalog.fits",
        f"{PHOENIXPATH}/grid/phoenix/catalog.fits",
    )


def build_phoenix():
    # Check if the directory exists and has any files
    os.makedirs(PHOENIXGRIDPATH, exist_ok=True)
    if (
        len(os.listdir(PHOENIXGRIDPATH)) == 65
    ):  # The directory exists and has files in it
        logger.debug(f"Found PHOENIX data in package in {PHOENIXGRIDPATH}.")
    else:
        logger.warning("No PHOENIX grid found, downloading grid.")
        download_phoenix_grid()
        logger.warning("PHEONIX grid downloaded.")


def phoenixcontext():
    """
    Decorator that temporarily sets the `PYSYN_CDBS` environment variable.

    Parameters
    ----------
    phoenixpath : str
        The value to temporarily set for the `PYSYN_CDBS` environment variable.

    Returns
    -------
    function
        A wrapper function that sets `PYSYN_CDBS` to `phoenixpath` before
        executing the decorated function and restores the original environment
        afterwards.

    Examples
    --------
    Using `set_pysyn_cdbs` to temporarily set `PYSYN_CDBS` for a function:

    >>> @set_pysyn_cdbs()
    ... def my_function():
    ...     # Within this function, os.environ["PYSYN_CDBS"] is set
    ...
    >>> my_function()
    >>> 'PYSYN_CDBS' in os.environ
    False
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with modified_environ(PYSYN_CDBS=PHOENIXPATH):
                return func(*args, **kwargs)

        return wrapper

    return decorator


@contextlib.contextmanager
def modified_environ(**update):
    """
    Temporarily updates the `os.environ` dictionary in-place and restores it upon exit.
    """
    env = os.environ
    original_state = env.copy()

    # Apply updates to the environment
    env.update(update)

    try:
        yield
    finally:
        # Restore original environment
        env.clear()
        env.update(original_state)


asqlog.setLevel("ERROR")


@phoenixcontext()
def get_phoenix_model(teff, logg=4.5, jmag=None, vmag=None):
    build_phoenix()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Extinction files not found in "
        )
        # Third-party
        import pysynphot
    logg1 = logg.value if isinstance(logg, u.Quantity) else logg
    star = pysynphot.Icat(
        "phoenix",
        teff.value if isinstance(teff, u.Quantity) else teff,
        0,
        logg1 if np.isfinite(logg1) else 5,
    )
    if (jmag is not None) & (vmag is None):
        star_norm = star.renorm(
            jmag, "vegamag", pysynphot.ObsBandpass("johnson,j")
        )
    elif (jmag is None) & (vmag is not None):
        star_norm = star.renorm(
            vmag, "vegamag", pysynphot.ObsBandpass("johnson,V")
        )
    else:
        raise ValueError("Input one of either `jmag` or `vmag`")
    star_norm.convert("Micron")
    star_norm.convert("flam")
    mask = (star_norm.wave >= 0.1) * (star_norm.wave <= 3)
    wavelength = star_norm.wave[mask] * u.micron
    wavelength = wavelength.to(u.angstrom)

    sed = star_norm.flux[mask] * u.erg / u.s / u.cm**2 / u.angstrom
    return wavelength, sed
