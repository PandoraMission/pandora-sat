# Third-party
import astropy.units as u
import numpy as np
import pytest

# First-party/Local
from pandorasat import PACKAGEDIR, PandoraSat, Target, __version__


def test_version():
    assert __version__ == "0.4.1"


def test_pandorasat():
    p = PandoraSat(ra=180 * u.deg, dec=0 * u.deg, theta=10 * u.deg)
    nirda = p.NIRDA
    visda = p.VISDA
    wavelength = np.linspace(0.1, 2, 1000) * u.micron
    nirda.sensitivity(wavelength)
    visda.sensitivity(wavelength)
    assert np.isclose(nirda.midpoint.value, 1.29750, atol=1e-3)
    assert np.isclose(visda.midpoint.value, 0.55399, atol=1e-3)
    return


@pytest.mark.remote_data
def test_trace():
    p = PandoraSat(ra=180 * u.deg, dec=0 * u.deg, theta=10 * u.deg)
    t = Target.from_gaia("GJ 436")
    nirda = p.NIRDA
    wavelength = np.linspace(0.1, 2, 6000) * u.micron
    spectrum = t.spectrum(wavelength)
    nirda.get_trace(
        wavelength,
        spectrum.value**0 * spectrum.unit,
        target_center=[40, 250],
    )
    return