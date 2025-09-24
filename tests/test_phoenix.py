# Standard library
import os

# First-party/Local
from pandorasat import phoenix

PHOENIXPATH = phoenix.PHOENIXPATH
os.environ["PYSYN_CDBS"] = PHOENIXPATH

# Third-party
import astropy.units as u  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
import stsynphot as stsyn  # noqa: E402


# testing get_vega
def test_get_vega():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )

    phoenix.download_vega()
    assert os.path.isfile(PHOENIXPATH + "/calspec/alpha_lyr_stis_011.fits")

    wav, spec = phoenix.load_vega()

    # Check that the loaded file is correct
    assert isinstance(wav, u.Quantity)
    assert isinstance(spec, u.Quantity)
    assert len(wav) > 0
    assert len(spec) > 0
    return


# test get_phoenix_model
def test_spectrum_generation():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    try:
        wavelength, sed = phoenix.get_phoenix_model(teff=5000, jmag=9)
    except Exception:
        pytest.fail("Failed to generate a spectrum")


# test get_phoenix_model units
def test_output_units():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    wavelength, sed = phoenix.get_phoenix_model(teff=5000, jmag=9)
    try:
        u.Quantity(wavelength, u.AA)
    except u.UnitConversionError:
        pytest.fail("Incorrect units")
    try:
        u.Quantity(sed, u.erg / (u.AA * u.s * u.cm**2))
    except u.UnitConversionError:
        pytest.fail("Incorrect units")


# test get_phoenix_model outputs in correct wavelength range and all SED values positive
def test_output_inrange():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    wavelength, sed = phoenix.get_phoenix_model(teff=5000, jmag=9)
    assert np.all(sed.value >= 0)
    assert np.all((wavelength.value >= 1000) & (wavelength.value <= 30000))


# test get_phoenix_model output length
def test_output_length():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    wavelength, sed = phoenix.get_phoenix_model(teff=5000, jmag=9)
    assert len(wavelength) == len(sed)


# test get_phoenix_model errors out of grid range
def test_out_of_grid_T():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    with pytest.raises(UnboundLocalError):
        wavelength, sed = phoenix.get_phoenix_model(teff=10000, jmag=9)


def test_out_of_grid_logg():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    with pytest.raises(stsyn.exceptions.ParameterOutOfBounds):
        wavelength, sed = phoenix.get_phoenix_model(teff=5000, logg=15.0, jmag=9)


def test_get_phoenix():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )

    wavelength, sed = phoenix.SED(teff=5000, jmag=9)
    assert len(wavelength) == len(sed)
    try:
        u.Quantity(wavelength, u.AA)
    except u.UnitConversionError:
        pytest.fail("Incorrect units")
    try:
        u.Quantity(sed, u.erg / (u.AA * u.s * u.cm**2))
    except u.UnitConversionError:
        pytest.fail("Incorrect units")


def test_get_benchmark():
    wavelength, sed = phoenix.load_benchmark()
    assert len(wavelength) == len(sed)
    try:
        u.Quantity(wavelength, u.AA)
    except u.UnitConversionError:
        pytest.fail("Incorrect units")
    try:
        u.Quantity(sed, u.erg / (u.AA * u.s * u.cm**2))
    except u.UnitConversionError:
        pytest.fail("Incorrect units")
