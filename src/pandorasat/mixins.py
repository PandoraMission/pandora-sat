# Standard library
import os

# Third-party
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# First-party/Local
from pandorasat import PACKAGEDIR, PANDORASTYLE
from pandorasat.utils import load_vega

__all__ = ["DetectorMixins"]


class DetectorMixins:
    def _add_trace_params(self, detector_name):
        fname = f"{PACKAGEDIR}/data/{detector_name}-wav-solution.fits"
        if not os.path.isfile(fname):
            raise ValueError(f"No wavelength solutions for `{self.name}`.")
        hdu = fits.open(fname)
        for idx in np.arange(1, hdu[1].header["TFIELDS"] + 1):
            name, unit = (
                hdu[1].header[f"TTYPE{idx}"],
                hdu[1].header[f"TUNIT{idx}"],
            )
            setattr(
                self, f"trace_{name}", hdu[1].data[name] * u.Quantity(1, unit)
            )
        self.trace_sensitivity *= hdu[1].header["SENSCORR"] * u.Quantity(
            1, hdu[1].header["CORRUNIT"]
        )

    def plot_sensitivity(self, ax=None):
        if ax is None:
            _, ax = plt.subplots()
        with plt.style.context(PANDORASTYLE):
            ax.plot(
                self.trace_wavelength.value,
                self.trace_sensitivity.value,
                c="k",
            )
            ax.set(
                xticks=np.linspace(*ax.get_xlim(), 9),
                xlabel=f"Wavelength [{self.trace_wavelength.unit.to_string('latex')}]",
                ylabel=f"Sensitivity [{self.trace_sensitivity.unit.to_string('latex')}]",
                title=self.name.upper(),
            )
            ax.spines[["right", "top"]].set_visible(True)
            if (self.trace_pixel.value != 0).any():
                ax_p = ax.twiny()
                ax_p.set(xticks=ax.get_xticks(), xlim=ax.get_xlim())
                ax_p.set_xlabel(xlabel="$\delta$ Pixel Position", color="grey")
                ax_p.set_xticklabels(
                    labels=list(
                        np.interp(
                            ax.get_xticks(),
                            self.trace_wavelength.value,
                            self.trace_pixel.value,
                        ).astype(int)
                    ),
                    rotation=45,
                    color="grey",
                )
        return ax

    def estimate_zeropoint(self):
        """Use Vega SED to estimate the zeropoint of the detector"""
        wavelength, spectrum = load_vega()
        sens = self.sensitivity(wavelength)
        zeropoint = np.trapz(spectrum * sens, wavelength) / np.trapz(
            sens, wavelength
        )
        return zeropoint

    def mag_from_flux(self, flux):
        """Convert flux to magnitude based on the zeropoint of the detector"""
        return -2.5 * np.log10(flux / self.estimate_zeropoint())

    def flux_from_mag(self, mag):
        """Convert magnitude to flux based on the zeropoint of the detector"""
        return self.estimate_zeropoint() * 10 ** (-mag / 2.5)
