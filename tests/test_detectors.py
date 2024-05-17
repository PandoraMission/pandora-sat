# Standard library
import pickle

# Third-party
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS, Sip

# First-party/Local
from pandorasat import PACKAGEDIR, PANDORASTYLE, TESTDIR, PandoraSat
from pandorasat.irdetector import NIRDetector
from pandorasat.visibledetector import VisibleDetector


# Check NIR and Visible detectors have the correct sensitivity
def test_pandorasat():
    p = PandoraSat()
    nirda = p.NIRDA
    visda = p.VISDA
    wavelength = np.linspace(0.1, 2, 1000) * u.micron
    nirda.sensitivity(wavelength)
    visda.sensitivity(wavelength)

    for detector in [visda, nirda]:
        for attr in ["trace_sensitivity", "trace_wavelength", "trace_pixel"]:
            assert hasattr(detector, attr)
            wcs = detector.get_wcs(150, 0)
            assert isinstance(wcs, WCS)
            assert isinstance(wcs.sip, Sip)
            detector.info


def test_gain():
    for detector in [VisibleDetector(), NIRDetector()]:
        r = detector.apply_gain(1e3 * u.DN)
        assert r.ndim == 0
        assert r.unit == u.electron

        r = detector.apply_gain([1e3] * u.DN)
        assert len(r) == 1
        assert r.ndim == 1

        r = detector.apply_gain([1e3, 1e3] * u.DN)
        assert len(r) == 2
        assert r.ndim == 1

        r = detector.apply_gain(1e3 * u.electron)
        assert r.ndim == 0
        assert r.unit == u.DN

        r = detector.apply_gain([1e3] * u.electron)
        assert len(r) == 1
        assert r.ndim == 1

        r = detector.apply_gain([1e3, 1e3] * u.electron)
        assert len(r) == 2
        assert r.ndim == 1


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


def test_plots():
    """Make some basic files and figures to show the detector information."""
    for detector in [VisibleDetector(), NIRDetector()]:
        with open(f"{TESTDIR}/output/{detector.name}.md", "w") as file:
            file.write(detector.info.to_markdown())
        with plt.style.context(PANDORASTYLE):
            detector.plot_sensitivity()
            plt.savefig(f"{TESTDIR}/output/{detector.name}.png", dpi=150)
            plt.close("all")
