# Third-party
import astropy.units as u
import numpy as np

# First-party/Local
from pandorasat import PandoraSat, __version__
from pandorasat.irdetector import NIRDetector
from pandorasat.visibledetector import VisibleDetector


def test_version():
    assert __version__ == "0.4.1"


# Check NIR and Visible detectors have the correct sensitivity
def test_pandorasat():
    p = PandoraSat()
    nirda = p.NIRDA
    visda = p.VISDA
    wavelength = np.linspace(0.1, 2, 1000) * u.micron
    nirda.sensitivity(wavelength)
    visda.sensitivity(wavelength)
    assert np.isclose(nirda.midpoint.value, 1.29750, atol=1e-3)
    assert np.isclose(visda.midpoint.value, 0.55399, atol=1e-3)
    return


# Check that NIR detector SNR mission requirements are met
def test_NIRDA_snr():
    wavelength = np.arange(0.1, 3, 0.0001)*u.micron
    spectrum = np.loadtxt("/Users/bhord/research/pandora/scratch/nir_spec.txt")
    mask = (wavelength.to(u.nm).value > (1300 - 22)) & (wavelength.to(u.nm).value < (1300 + 22))
    signal = (np.trapz(NIRDetector().sensitivity(wavelength) * spectrum * mask, wavelength) * (2 * u.hour).to(u.second)).value

    assert (signal/np.sqrt(signal)) > 6000
    return


# Check that Visible detector SNR mission requirements are met
def test_VISDA_snr():
    wavelength = np.arange(0.1, 3, 0.0001)*u.micron
    spectrum = np.loadtxt("/Users/bhord/research/pandora/scratch/vis_spec.txt")
    signal = (np.trapz(VisibleDetector().sensitivity(wavelength) * spectrum, wavelength) * (10 * u.minute).to(u.electron)).value

    assert (signal/np.sqrt(signal)) > 1000
    return
