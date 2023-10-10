"""Holds metadata and methods on Pandora VISDA"""
# Standard library
from glob import glob
from dataclasses import dataclass

# Third-party
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io import fits
from tqdm import tqdm

from . import PACKAGEDIR
from .hardware import Optics
from .utils import photon_energy, load_vega


@dataclass
class VisibleDetector:
    """Holds information on the Pandora Visible Detector

    Attributes
    ----------

    name: str
        Name of the detector. This will determine which files are loaded. This
        will be `"visda"` for this detector
    pixel_scale: float
        The pixel scale of the detector in arcseconds/pixel
    pixel_size: float
        The pixel size in microns/mm
    """
    pixel_scale: float
    pixel_size: float

    def __post_init__(self):
        """Some detector specific functions to run on initialization"""
        self.shape = (2048, 2048)

        self.flat = fits.open(
            np.sort(
                np.atleast_1d(glob(f"{PACKAGEDIR}/data/flatfield_VISDA*.fits"))
            )[-1]
        )[1].data
        if hasattr(self, "fieldstop_radius"):
            C, R = (
                np.mgrid[
                    : self.shape[0],
                    : self.shape[1],
                ]
                - np.hstack(
                    [
                        self.shape[0],
                        self.shape[1],
                    ]
                )[:, None, None]
                / 2
            )
            r = (self.fieldstop_radius / self.pixel_scale).to(u.pix).value
            self.fieldstop = ~((np.abs(C) >= r) | (np.abs(R) >= r))

    @property
    def pixel_scale(self):
        """Pixel scale of the detector"""
        return 0.78 * u.arcsec / u.pixel

    @property
    def pixel_size(self):
        """Size of a pixel"""
        return 6.5 * u.um / u.pixel

    @property
    def naxis1(self):
        """WCS's are COLUMN major, so naxis1 is the number of columns"""
        return self.shape[1] * u.pixel

    @property
    def naxis2(self):
        """WCS's are COLUMN major, so naxis2 is the number of rows"""
        return self.shape[0] * u.pixel

    @property
    def _dispersion_df(self):
        return pd.read_csv(f"{PACKAGEDIR}/data/pixel_vs_wavelength_vis.csv")

    @property
    def dark(self):
        return 1 * u.electron / u.second

    @property
    def read_noise(self):
        return 1.5 * u.electron

    @property
    def bias(self):
        return 100 * u.electron

    @property
    def integration_time(self):
        return 0.2 * u.second

#    @property
#    def fieldstop_radius(self):
#        return 0.3 * u.deg

    def throughput(self, wavelength):
        return wavelength.value**0 * 0.714

    def get_subarray(
        self,
        cat,
        corner,
        nreads=50,
        nframes=10,
        #        freeze_dimensions=[2, 3],
        quiet=False,
        time_series_generators=None,
        include_noise=True,
    ):
        f = np.zeros((nframes, *self.subarray_size), dtype=int)
        time = self.time[: nreads * nframes]
        time = np.asarray([time[idx::nreads] for idx in range(nreads)]).mean(
            axis=0
        )
        for idx, m in cat.iterrows():
            if time_series_generators is None:
                tsgenerator = lambda x: 1  # noqa
            else:
                tsgenerator = time_series_generators[idx]
            if tsgenerator is None:
                tsgenerator = lambda x: 1  # noqa

            prf = self.prf(
                row=m.vis_row,
                col=m.vis_column,
                corner=corner,
                #                freeze_dimensions=freeze_dimensions,
            )
            for tdx in tqdm(
                range(nframes),
                desc="Times",
                position=0,
                leave=True,
                disable=quiet,
            ):
                f[tdx] += self.apply_gain(
                    u.Quantity(
                        np.random.poisson(
                            prf
                            * tsgenerator(time[tdx])
                            * m.vis_counts
                            * nreads
                            * self.integration_time.to(u.second).value
                        ),
                        dtype=int,
                        unit=u.DN,
                    )
                ).value
        if include_noise:
            for tdx in range(nframes):
                f[tdx] += (
                    self.get_background_light_estimate(
                        cat.loc[0, "ra"],
                        cat.loc[0, "dec"],
                        nreads * self.integration_time,
                        self.subarray_size,
                    )
                ).value.astype(int)
                f[tdx] += np.random.normal(
                    loc=self.bias.value,
                    scale=self.read_noise.value,
                    size=self.subarray_size,
                ).astype(int)
                f[tdx] += np.random.poisson(
                    lam=(self.dark * self.integration_time * nreads).value,
                    size=self.subarray_size,
                ).astype(int)

        apers = np.asarray(
            [
                self.get_aper(
                    row=m.vis_row,
                    col=m.vis_column,
                    corner=corner,
                    # freeze_dimensions=freeze_dimensions,
                )
                for idx, m in cat.iterrows()
            ]
        )
        return time, f, apers

    def get_aper(
        self,
        row=0,
        col=0,
        corner=(0, 0),
        shape=None,  # , freeze_dimensions=[2, 3]
    ):
        if shape is None:
            shape = self.subarray_size
        aper = np.zeros(shape)
        for idx in np.arange(0, 1, 0.1):
            for jdx in np.arange(0, 1, 0.1):
                aper += self.prf(
                    row=row + idx,
                    col=col + jdx,
                    corner=corner,
                    shape=shape,
                    #                    freeze_dimensions=freeze_dimensions,
                )
        aper /= 100
        return aper > 0.005

    def qe(self, wavelength):
        """
        Calculate the quantum efficiency of the detector.

        Parameters:
            wavelength (npt.NDArray): Wavelength in microns as `astropy.unit`

        Returns:
            qe (npt.NDArray): Array of the quantum efficiency of the detector
        """
        raise NotImplementedError

    def sensitivity(self, wavelength):
        sed = 1 * u.erg / u.s / u.cm**2 / u.angstrom
        E = photon_energy(wavelength)
        telescope_area = np.pi * (Optics.mirror_diameter / 2) ** 2
        photon_flux_density = (
            telescope_area * sed * self.throughput(wavelength) / E
        ).to(u.photon / u.second / u.angstrom) * self.qe(wavelength)
        sensitivity = photon_flux_density / sed
        return sensitivity

    @property
    def midpoint(self):
        """Mid point of the sensitivity function"""
        w = np.arange(0.1, 3, 0.005) * u.micron
        return np.average(w, weights=self.sensitivity(w))

    def estimate_zeropoint(self):
        """Use Vega SED to estimate the zeropoint of the detector"""
        wavelength, spectrum = load_vega()
        sens = self.sensitivity(wavelength)
        zeropoint = np.trapz(spectrum * sens, wavelength) / np.trapz(
            sens, wavelength
        )
        return zeropoint

    def mag_from_flux(self, flux):
        return -2.5 * np.log10(flux / self.zeropoint)

    def flux_from_mag(self, mag):
        return self.zeropoint * 10 ** (-mag / 2.5)
