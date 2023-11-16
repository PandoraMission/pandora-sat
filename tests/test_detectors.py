# Standard library
import pickle

# Third-party
import astropy.units as u
import numpy as np

# First-party/Local
from pandorasat import PACKAGEDIR, PandoraSat, __version__
from pandorasat.irdetector import NIRDetector
from pandorasat.visibledetector import VisibleDetector


def test_version():
    assert __version__ == "0.5.0"


# Check NIR and Visible detectors have the correct sensitivity
def test_pandorasat():
    p = PandoraSat()
    nirda = p.NIRDA
    visda = p.VISDA
    wavelength = np.linspace(0.1, 2, 1000) * u.micron
    nirda.sensitivity(wavelength)
    visda.sensitivity(wavelength)
    return


# Check that NIR and Visible detector SNR mission requirements are met
def test_detector_snr():
    # Fetch test star spectrum to test with
    with open(f"{PACKAGEDIR}/data/test-star.p", "rb") as f:
        model_spec = pickle.load(f)

    # spec_res = 1.3*u.micron / 30
    # wavelength = np.arange(0.1, 3, spec_res.value)*u.micron
    wavelength = (
        np.loadtxt(
            f"{PACKAGEDIR}/data/dichroic-transmission.csv",
            unpack=True,
            skiprows=1,
            delimiter=",",
        )[0]
        * u.nm
    ).to(u.micron)
    spectrum = np.interp(
        wavelength, model_spec["wavelength"], model_spec["spectrum"]
    )

    # Testing NIRDA SNR
    mask = (wavelength.to(u.nm).value > (1300 - 22)) & (
        wavelength.to(u.nm).value < (1300 + 22)
    )
    signal = (
        np.trapz(
            NIRDetector().sensitivity(wavelength) * spectrum * mask, wavelength
        )
        * (2 * u.hour).to(u.second)
    ).value
    assert (signal / np.sqrt(signal)) > 6000

    # Testing VISDA SNR
    signal = (
        np.trapz(
            VisibleDetector().sensitivity(wavelength) * spectrum, wavelength
        )
        * (10 * u.minute).to(u.second)
    ).value
    assert (signal / np.sqrt(signal)) > 1000
    return
