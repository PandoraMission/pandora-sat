# Standard library
import contextlib
import functools
import os
import shutil
import tarfile
import warnings
from glob import glob

# Third-party
import astropy.units as u
import numpy as np
import requests
from astroquery import log as asqlog
from tqdm import tqdm

from . import CACHEDIR, PHOENIXGRIDPATH, PHOENIXPATH, logger

__all__ = [
    "download_phoenix_grid",
    "phoenixcontext",
    "build_phoenix",
    "get_phoenix_model",
    "download_vega",
]


def download_file(file_url, file_path):
    # Download the file from `file_url` and save it locally under `file_path`
    with requests.get(file_url, stream=True) as r:
        r.raise_for_status()
        with open(file_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    # astropy_download_file(file_url, cache=True, show_progress=False, pkgname='pandorasat')


def download_vega():
    """
    Downloads the Vega calibration file for STSynPhot and moves it to the proper directory, if one does not already exist.
    Ensures the file is set in the config for synphot.
    """
    # Check if the file already exists in the right location
    if os.path.isfile(PHOENIXPATH + "calspec/alpha_lyr_stis_011.fits"):
        logger.debug(f"Found Vega spectrum in {PHOENIXGRIDPATH}calspec.")
    else:
        logger.warning(
            "No Vega spectrum found, downloading from STScI website."
        )
        os.makedirs(PHOENIXPATH + "calspec", exist_ok=True)
        download_file(
            "http://ssb.stsci.edu/cdbs/calspec/alpha_lyr_stis_011.fits",
            PHOENIXPATH + "calspec/alpha_lyr_stis_011.fits",
        )
        logger.warning("Vega spectrum downloaded.")
    # Third-party
    import synphot

    synphot.conf.vega_file = PHOENIXPATH + "calspec/alpha_lyr_stis_011.fits"


def download_phoenix_grid():
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
        download_file(f"{url}{filename}", f"{PHOENIXGRIDPATH}{filename}")
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
        f"{PHOENIXPATH}grid/phoenix/catalog.fits",
    )
    download_vega()


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
    """
    Function that interpolates the PHOENIX grid to a given temperature and surface gravity.
    Returns a SED for a star normalized by its Johnson J or V magnitude via STSynPhot.

    Parameters
    ----------
    teff : float
        The effective temperature of the star (Kelvin).
    logg : float
        The log surface gravity of the star (log cgs).
    jmag : float
        The Johnson J-band magnitude of the star.
    vmag : float
        The Johnson V-band magnitude of the star.

    Returns
    -------
    wavelength : array
        An array of wavelengths from 1,000 to 30,000 Angstroms.
    sed : array
        The SED of the star, in units of ergs s^-1 cm^-2 Angstrom^-1.
    """
    # Third-party
    import stsynphot as stsyn
    from synphot import units as su

    build_phoenix()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Extinction files not found in "
        )
        # Third-party

    logg1 = logg.value if isinstance(logg, u.Quantity) else logg
    star = stsyn.grid_to_spec(
        "phoenix",
        teff.value if isinstance(teff, u.Quantity) else teff,
        0,
        logg1 if np.isfinite(logg1) else 5,
    )
    if (jmag is not None) & (vmag is None):
        star_norm = star.normalize(
            jmag * su.VEGAMAG,
            band=stsyn.band("johnson,j"),
            vegaspec=stsyn.Vega,
        )
    elif (jmag is None) & (vmag is not None):
        star_norm = star.normalize(
            vmag * su.VEGAMAG,
            band=stsyn.band("johnson,v"),
            vegaspec=stsyn.Vega,
        )
    else:
        raise ValueError("Input one of either `jmag` or `vmag`")

    wave = star_norm.waveset.to(u.micron)
    mask = (wave >= 0.1 * u.micron) * (wave <= 3 * u.micron)

    sed = (
        star_norm(wave, flux_unit="flam")[mask]
        / su.FLAM
        * u.erg
        / u.s
        / u.cm**2
        / u.angstrom
    )

    wavelength = wave[mask]
    wavelength = wavelength.to(u.angstrom)

    return wavelength, sed


@phoenixcontext()
def load_vega():
    """Loads a spectrum of Vega using synphot

    Returns
    -------
    wavelength : array
        Wavelength array
    sed : array
        The SED of the Vega, in units of ergs s^-1 cm^-2 Angstrom^-1.
    """
    # Third-party
    from synphot import SourceSpectrum

    download_vega()
    # Third-party

    vega = SourceSpectrum.from_vega()
    wavelength, spectrum = vega.waveset, vega(vega.waveset, flux_unit="flam")

    spectrum = spectrum.value * u.erg / u.cm**2 / u.s / u.angstrom
    return wavelength, spectrum
