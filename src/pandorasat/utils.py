# Standard library
import warnings
from functools import lru_cache

# Third-party
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.constants import c, h
from astropy.coordinates import Distance, SkyCoord
from astropy.io import fits
from astropy.time import Time
from astroquery.gaia import Gaia

from . import PACKAGEDIR, __version__
from .phoenix import get_phoenix_model


def SED(teff, logg=4.5, jmag=None, vmag=None):
    """Gives a model SED for a given Teff, logg and magnitude."""
    return get_phoenix_model(teff, logg=logg, jmag=jmag, vmag=vmag)


@lru_cache
def get_sky_catalog(
    ra: float,
    dec: float,
    radius: float = 0.155,
    gbpmagnitude_range: tuple = (-3, 20),
    limit=None,
    gaia_keys: list = [],
    time: Time = Time.now(),
) -> dict:
    """
    Gets a catalog of coordinates on the sky based on an input RA, Dec, and radius as well as
    a magnitude range for Gaia. The user can also specify additional keywords to be grabbed
    from Gaia catalog.

    Parameters
    ----------
    ra : float
        Right Ascension of the center of the query radius in degrees.
    dec : float
        Declination of the center of the query radius in degrees.
    radius : float
        Radius centered on ra and dec that will be queried in degrees.
    gbpmagnitude_range : tuple
        Magnitude limits for the query. Targets outside of this range will not be included in
        the final output dictionary.
    limit : int
        Maximum number of targets from query that will be included in output dictionary. If a
        limit is specified, targets will be included based on proximity to specified ra and dec.
    gaia_keys : list
        List of additional Gaia archive columns to include in the final output dictionary.
    time : astropy.Time object
        Time at which to evaluate the positions of the targets in the output dictionary.

    Returns
    -------
    cat : dict
        Dictionary of values from the Gaia archive for each keyword.
    """

    base_keys = [
        "source_id",
        "ra",
        "dec",
        "parallax",
        "pmra",
        "pmdec",
        "radial_velocity",
        "ruwe",
        "phot_bp_mean_mag",
        "teff_gspphot",
        "logg_gspphot",
        "phot_g_mean_flux",
        "phot_g_mean_mag",
    ]

    all_keys = base_keys + gaia_keys

    query_str = f"""
    SELECT {f'TOP {limit} ' if limit is not None else ''}* FROM (
        SELECT gaia.{', gaia.'.join(all_keys)}, dr2.teff_val AS dr2_teff_val,
        dr2.rv_template_logg AS dr2_logg, tmass.j_m, tmass.j_msigcom, tmass.ph_qual, DISTANCE(
        POINT({u.Quantity(ra, u.deg).value}, {u.Quantity(dec, u.deg).value}),
        POINT(gaia.ra, gaia.dec)) AS ang_sep,
        EPOCH_PROP_POS(gaia.ra, gaia.dec, gaia.parallax, gaia.pmra, gaia.pmdec,
        gaia.radial_velocity, gaia.ref_epoch, 2000) AS propagated_position_vector
        FROM gaiadr3.gaia_source AS gaia
        JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xmatch USING (source_id)
        JOIN gaiadr3.dr2_neighbourhood AS xmatch2 ON gaia.source_id = xmatch2.dr3_source_id
        JOIN gaiadr2.gaia_source AS dr2 ON xmatch2.dr2_source_id = dr2.source_id
        JOIN gaiadr3.tmass_psc_xsc_join AS xjoin USING (clean_tmass_psc_xsc_oid)
        JOIN gaiadr1.tmass_original_valid AS tmass ON
        xjoin.original_psc_source_id = tmass.designation
        WHERE 1 = CONTAINS(
        POINT({u.Quantity(ra, u.deg).value}, {u.Quantity(dec, u.deg).value}),
        CIRCLE(gaia.ra, gaia.dec, {(u.Quantity(radius, u.deg) + 50*u.arcsecond).value}))
        AND gaia.parallax IS NOT NULL
        AND gaia.phot_bp_mean_mag > {gbpmagnitude_range[0]}
        AND gaia.phot_bp_mean_mag < {gbpmagnitude_range[1]}) AS subquery
    WHERE 1 = CONTAINS(
    POINT({u.Quantity(ra, u.deg).value}, {u.Quantity(dec, u.deg).value}),
    CIRCLE(COORD1(subquery.propagated_position_vector), COORD2(subquery.propagated_position_vector), {u.Quantity(radius, u.deg).value}))
    ORDER BY ang_sep ASC
    """
    job = Gaia.launch_job_async(query_str, verbose=False)
    tbl = job.get_results()
    if len(tbl) == 0:
        raise ValueError("Could not find matches.")
    plx = tbl["parallax"].value.filled(fill_value=0)
    plx[plx < 0] = 0
    cat = {
        "jmag": tbl["j_m"].data.filled(np.nan),
        "bmag": tbl["phot_bp_mean_mag"].data.filled(np.nan),
        "gmag": tbl["phot_g_mean_mag"].data.filled(np.nan),
        "gflux": tbl["phot_g_mean_flux"].data.filled(np.nan),
        "ang_sep": tbl["ang_sep"].data.filled(np.nan) * u.deg,
    }
    cat["teff"] = (
        tbl["teff_gspphot"].data.filled(
            tbl["dr2_teff_val"].data.filled(np.nan)
        )
        * u.K
    )
    cat["logg"] = tbl["logg_gspphot"].data.filled(
        tbl["dr2_logg"].data.filled(np.nan)
    )
    cat["RUWE"] = tbl["ruwe"].data.filled(99)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cat["coords"] = SkyCoord(
            ra=tbl["ra"].value.data * u.deg,
            dec=tbl["dec"].value.data * u.deg,
            pm_ra_cosdec=tbl["pmra"].value.filled(fill_value=0)
            * u.mas
            / u.year,
            pm_dec=tbl["pmdec"].value.filled(fill_value=0) * u.mas / u.year,
            obstime=Time.strptime("2016", "%Y"),
            distance=Distance(parallax=plx * u.mas, allow_negative=True),
            radial_velocity=tbl["radial_velocity"].value.filled(fill_value=0)
            * u.km
            / u.s,
        ).apply_space_motion(time)
    cat["source_id"] = np.asarray(
        [f"Gaia DR3 {i}" for i in tbl["source_id"].value.data]
    )
    for key in gaia_keys:
        cat[key] = tbl[key].data.filled(np.nan)
    return cat


def photon_energy(wavelength):
    """Converts photon wavelength to energy."""
    return ((h * c) / wavelength) * 1 / u.photon


def simulate_flatfield(stddev=0.005, seed=777):
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


def load_vega():
    wavelength, spectrum = np.loadtxt(
        f"{PACKAGEDIR}/data/vega.csv", delimiter=","
    ).T
    wavelength *= u.angstrom
    spectrum *= u.erg / u.cm**2 / u.s / u.angstrom
    return wavelength, spectrum


def wavelength_to_rgb(wavelength, gamma=0.8):

    """This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    """

    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    R *= 255
    G *= 255
    B *= 255
    return np.asarray((int(R), int(G), int(B))) / 256


def make_pixel_files():
    from .irdetector import NIRDetector
    from .visibledetector import VisibleDetector

    VISDA = VisibleDetector()
    NIRDA = NIRDetector()
    df = pd.read_csv(f"{PACKAGEDIR}/data/pixel_vs_wavelength.csv")
    pixel = np.round(np.arange(-400, 80, 0.5), 5) * u.pixel
    wav = (
        np.polyval(np.polyfit(df.Pixel, df.Wavelength, 3), pixel.value)
        * u.micron
    )

    sens = NIRDA.sensitivity(wav)
    mask = sens / sens.max() > 0.0001
    corr = np.trapz(sens, wav)
    hdu = fits.TableHDU.from_columns(
        [
            fits.Column(
                name="pixel",
                format="D",
                array=pixel.value[mask],
                unit=pixel.unit.to_string(),
            ),
            fits.Column(
                name="wavelength",
                format="D",
                array=wav.value[mask],
                unit=wav.unit.to_string(),
            ),
            fits.Column(
                name="sensitivity",
                format="D",
                array=(sens[mask] / corr),
                unit=(sens / corr).unit.to_string(),
            ),
        ]
    )
    hdu.header.append(
        fits.Card("SENSCORR", corr.value, "correction to apply to sensitivity")
    )
    hdu.header.append(
        fits.Card("CORRUNIT", corr.unit.to_string(), "units of correction")
    )
    hdu.writeto(f"{PACKAGEDIR}/data/nirda-wav-solution.fits", overwrite=True)

    wav = np.arange(0.25, 1.3, 0.01) * u.micron
    pixel = np.zeros(len(wav)) * u.pixel
    sens = VISDA.sensitivity(wav)
    corr = np.trapz(sens, wav)

    hdu = fits.TableHDU.from_columns(
        [
            fits.Column(
                name="pixel",
                format="D",
                array=pixel.value,
                unit=pixel.unit.to_string(),
            ),
            fits.Column(
                name="wavelength",
                format="D",
                array=wav.value,
                unit=wav.unit.to_string(),
            ),
            fits.Column(
                name="sensitivity",
                format="D",
                array=(sens / corr),
                unit=(sens / corr).unit.to_string(),
            ),
        ]
    )
    hdu.header.append(
        fits.Card("SENSCORR", corr.value, "correction to apply to sensitivity")
    )
    hdu.header.append(
        fits.Card("CORRUNIT", corr.unit.to_string(), "units of correction")
    )
    hdu.writeto(f"{PACKAGEDIR}/data/visda-wav-solution.fits", overwrite=True)
    return
