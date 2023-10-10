# Standard library
import os
from glob import glob

# Third-party
import astropy.units as u
import numpy as np

# First-party/Local
from pandorasat import utils, PACKAGEDIR


# test photon_energy
def test_photon_energy():
    assert np.isclose(utils.photon_energy(0.21 * u.m).value, 9.41e-25, atol=1e-25)
    return


# test get_flatfield
def test_get_flatfield():
    utils.get_flatfield()

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
