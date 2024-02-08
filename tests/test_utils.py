# Standard library
import os
from glob import glob

# Third-party
import astropy.units as u
import numpy as np
import pytest

# First-party/Local
from pandorasat import PACKAGEDIR, utils


# test photon_energy
def test_photon_energy():
    assert np.isclose(
        utils.photon_energy(0.21 * u.m).value, 9.41e-25, atol=1e-25
    )
    return


# test get_flatfield
@pytest.mark.skip(
    reason="Creates a large file, we do not need this right now."
)
def test_simulate_flatfield():
    utils.simulate_flatfield()

    # Check that files were made for NIRDA and VISDA
    assert len(glob(f"{PACKAGEDIR}/data/flatfield_NIRDA*.fits")) > 0
    assert len(glob(f"{PACKAGEDIR}/data/flatfield_VISDA*.fits"))
    return


# test load_vega
def test_load_vega():
    # Check that the file exists
    assert os.path.isfile(f"{PACKAGEDIR}/data/vega.dat")

    wav, spec = utils.load_vega()

    # Check that the loaded file is correct
    assert isinstance(wav, u.Quantity)
    assert isinstance(spec, u.Quantity)
    assert len(wav) > 0
    assert len(spec) > 0
    return


def test_get_phoenix():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )

    wavelength, sed = utils.SED(teff=5000, jmag=9)
    assert len(wavelength) == len(sed)
    try:
        u.Quantity(wavelength, u.AA)
    except u.UnitConversionError:
        pytest.fail("Incorrect units")
    try:
        u.Quantity(sed, u.erg / (u.AA * u.s * u.cm**2))
    except u.UnitConversionError:
        pytest.fail("Incorrect units")
