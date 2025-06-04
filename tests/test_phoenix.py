# Standard library
import os
from glob import glob

# Third-party
import astropy.units as u
import numpy as np
import pytest

# First-party/Local
from pandorasat import PACKAGEDIR, phoenix

# testing get_vega
def test_get_vega():
    PHOENIXPATH = phoenix.PHOENIXPATH
    phoenix.get_vega()
    assert os.path.isfile(PHOENIXPATH+'/calspec/alpha_lyr_stis_010.fits')

# test get_phoenix_model
def test_spectrum_generation():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    try:
        wavelength, sed = phoenix.get_phoenix_model(teff=5000, jmag=9)
    except: pytest.fail("Failed to generate a spectrum")

# test get_phoenix_model units
def test_output_units():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    wavelength, sed = phoenix.get_phoenix_model(teff=5000, jmag=9)
    assert len(wavelength) == len(sed)
    try:
        u.Quantity(wavelength, u.AA)
    except u.UnitConversionError:
        pytest.fail("Incorrect units")
    try:
        u.Quantity(sed, u.erg / (u.AA * u.s * u.cm**2))
    except u.UnitConversionError:
        pytest.fail("Incorrect units")

# test get_phoenix_model errors out of grid range
def test_out_of_grid_T():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    with pytest.raises(UnboundLocalError): wavelength, sed = phoenix.get_phoenix_model(teff=10000, jmag=9)

def test_out_of_grid_logg():
    if os.getenv("GITHUB_ACTIONS") == "true":
        pytest.skip(
            "Skipping this test on GitHub Actions this downloads a database of stellar models."
        )
    with pytest.raises(stsyn.exceptions.ParameterOutOfBounds): wavelength, sed = phoenix.get_phoenix_model(teff=5000, logg=15., jmag=9)