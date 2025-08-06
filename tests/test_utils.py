# Standard library
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


@pytest.mark.remote_data
def test_get_sky_catalog():
    # Works with no units
    cat = utils.get_sky_catalog(ra=210.8023, dec=54.349, radius=0.05)
    assert isinstance(cat, dict)
    assert np.all(
        [
            k in cat.keys()
            for k in [
                "teff",
                "logg",
                "jmag",
                "bmag",
                "RUWE",
                "ang_sep",
                "coords",
                "source_id",
            ]
        ]
    )
    assert len(cat["coords"]) > 1

    # Works with units
    cat = utils.get_sky_catalog(
        ra=210.8023 * u.deg, dec=54.349 * u.deg, radius=0.05 * u.deg
    )
    assert isinstance(cat, dict)
    assert np.all(
        [
            k in cat.keys()
            for k in [
                "teff",
                "logg",
                "jmag",
                "bmag",
                "RUWE",
                "ang_sep",
                "coords",
                "source_id",
            ]
        ]
    )
    assert len(cat["coords"]) > 1

    # Can return the top 1 hit
    cat = utils.get_sky_catalog(
        ra=210.8023 * u.deg, dec=54.349 * u.deg, radius=0.02 * u.deg, limit=1
    )
    assert len(cat["coords"]) == 1
