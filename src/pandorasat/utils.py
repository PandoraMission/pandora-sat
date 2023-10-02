# Standard library
import os

# Third-party
import astropy.units as u
import numpy as np
from astropy.constants import c, h
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import fits
from astropy.time import Time
from astroquery import log as asqlog

from . import PACKAGEDIR, __version__

phoenixpath = f"{PACKAGEDIR}/data/phoenix"
os.environ["PYSYN_CDBS"] = phoenixpath

asqlog.setLevel("ERROR")
# Third-party


def photon_energy(wavelength):
    return ((h * c) / wavelength) * 1 / u.photon


def get_jitter(
    rowstd: float = 1,
    colstd: float = 0.3,
    thetastd: float = 0.0005,
    correlation_time=1 * u.second,
    nframes=200,
    frame_time=0.2 * u.second,
    seed=None,
):
    """Returns some random, time correlated jitter time-series.

    Parameters:
    ----------
    rowstd: float
        Standard deviation of jitter in pixels in row/y axis
    colstd: float
        Standard deviation of jitter in pixels in col/x axis
    thetastd: float
        Standard deviation of jitter in degrees in y axis
    correlation_time: float
        The timescale over which data is correlated in seconds.
        Increase this value for smoother time-series
    nframes: int
        Number of frames to generate
    frame_time: float
        The time spacing for each frame
    seed: Optional, int
        Optional seed for random walk

    Returns:
    --------
    time : np.ndarray
        Time array in seconds
    row: np.ndarray
        Jitter in the row/y axis in pixels
    col: np.ndarray
        Jitter in the column/x axis in pixels
    theta: np.ndarray
        Jitter in angle in degrees
    """
    time = np.arange(nframes) * frame_time  # noqa:F811
    tstd = (correlation_time / frame_time).value

    def jitter_func(std):
        f = np.random.normal(0, std, size=nframes)
        return convolve(f, Gaussian1DKernel(tstd)) * tstd**0.5

    jitter = []
    for idx, std, unit in zip(
        [0, 1, 2], [rowstd, colstd, thetastd], [u.pixel, u.pixel, u.deg]
    ):
        if seed is not None:
            np.random.seed(seed + idx)
        jitter.append(jitter_func(std) * unit)

    return time, *jitter


def get_flatfield(stddev=0.005, seed=777):
    np.random.seed(seed)
    """ This generates and writes a dummy flatfield file. """
    for detector in ["VISDA", "NIRDA"]:
        hdr = fits.Header()
        hdr["AUTHOR"] = "Christina Hedges"
        hdr["VERSION"] = __version__
        hdr["DATE"] = Time.now().strftime("%d-%m-%Y")
        hdr["STDDEV"] = stddev
        hdu0 = fits.PrimaryHDU(header=hdr)
        hdulist = fits.HDUList(
            [
                hdu0,
                fits.CompImageHDU(
                    data=np.random.normal(1, stddev, (2048, 2048)), name="FLAT"
                ),
            ]
        )
        hdulist.writeto(
            f"{PACKAGEDIR}/data/flatfield_{detector}_{Time.now().strftime('%Y-%m-%d')}.fits",
            overwrite=True,
            checksum=True,
        )
    return
