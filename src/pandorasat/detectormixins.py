# Standard library
import warnings

# Third-party
import astropy.units as u
import numpy as np

__all__ = ["DetectorMixins"]


class DetectorMixins:
    def qe(self, wavelength):
        """
        Calculate the quantum efficiency of the detector.

        Parameters
        ----------
        wavelength : npt.NDArray
            Wavelength in microns as `astropy.unit`

        Returns
        -------
        qe : npt.NDArray
            Array of the quantum efficiency of the detector
        """
        return self.reference.get_qe(wavelength=wavelength)

    def sensitivity(self, wavelength):
        """
        Calulate the sensitivity of the detector.

        Parameters
        ----------
        wavelength : npt.NDArray
            Wavelength in microns as `astropy.unit`

        Returns
        -------
        sensitivity : npt.NDArray
            Array of the sensitivity of the detector
        """
        return self.reference.get_sensitivity(wavelength=wavelength)

    def throughput(self, wavelength):
        """
        Calulate the throughput of the detector.

        Parameters
        ----------
        wavelength : npt.NDArray
            Wavelength in microns as `astropy.unit`

        Returns
        -------
        sensitivity : npt.NDArray
            Array of the throughput of the detector
        """
        return self.reference.get_throughput(wavelength=wavelength)

    def estimate_zeropoint(self):
        """
        Calulate the zeropoint of the detector.
        """
        warnings.warn(
            "'esimate_zeropoint' is deprecated and will be removed in a future version. "
            "Use 'zeropoint' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.zeropoint

    @property
    def zeropoint(self):
        """
        Calulate the zeropoint of the detector.
        """
        return self.reference.get_zeropoint()

    @property
    def dark(self):
        """Dark Noise"""
        return self.reference.get_dark()

    @property
    def readnoise(self):
        """Read Noise"""
        return self.reference.get_readnoise()

    @property
    def gain(self):
        """Gain"""
        return self.reference.get_gain()

    def apply_gain(self, values: u.Quantity):
        """Applies a single gain value"""
        if not isinstance(values, u.Quantity):
            raise ValueError("Must pass a quantity.")
        if values.unit == u.electron:
            return values / self.gain
        if values.unit == u.DN:
            return values * self.gain

    @property
    def bias(self):
        """Gain"""
        return self.reference.get_bias()

    def get_wcs(
        self,
        ra,
        dec,
        theta=u.Quantity(0, unit="degree"),
        distortion=True,
    ):
        """Get the World Coordinate System for a detector as an astropy.wcs.WCS object, given pointing parameters.
        This method only updates the CRVAL and PC parameters, the rest of the WCS is set by reference products
        within this package.

        Parameters:
        -----------
        target_ra: astropy.units.Quantity
            The target RA in degrees
        target_dec: astropy.units.Quantity
            The target Dec in degrees
        theta: astropy.units.Quantity
            The observatory angle in degrees

        Returns:
        --------
        wcs: astropy.wcs.WCS
            World Coordinate System object
        """
        return self.reference.get_wcs(
            target_ra=u.Quantity(ra, "deg"),
            target_dec=u.Quantity(dec, "deg"),
            theta=u.Quantity(theta, "deg"),
            distortion=distortion,
        )

    def flux_to_mag(self, flux):
        """Convert flux to magnitude based on the zeropoint of the detector"""
        if not isinstance(flux, u.Quantity):
            raise ValueError("Must pass flux as a quantity.")
        if flux.unit == u.electron / u.second:
            # User has passed band pass integrated flux, but this is not normalized correctly
            wavelength = (np.linspace(0.1, 3, 10000) * u.micron).to(u.AA)
            norm = np.trapz(self.sensitivity(wavelength), wavelength)
            return -2.5 * np.log10((flux / norm) / self.zeropoint)
        else:
            raise ValueError(
                f"Must pass units of flux: {(u.electron / u.second).to_string()}."
            )

    def average_flux_density_to_mag(self, average_flux_density):
        """Convert average flux density to magnitude based on the zeropoint of the detector"""
        if not isinstance(average_flux_density, u.Quantity):
            raise ValueError("Must pass flux as a quantity.")
        if average_flux_density.unit == u.erg / u.AA / u.s / u.cm**2:
            return -2.5 * np.log10(average_flux_density / self.zeropoint)
        else:
            raise ValueError(
                f"Must pass units of average flux density: {(u.erg / u.AA / u.s / u.cm).to_string()}."
            )

    def mag_to_flux(self, mag):
        """Convert magnitude to flux based on the zeropoint of the detector"""
        if not isinstance(mag, u.Quantity):
            mag = u.Quantity(mag, u.dimensionless_unscaled)
        if mag.unit != u.dimensionless_unscaled:
            raise ValueError("Magnitude must have dimensionless units.")
        wavelength = (np.linspace(0.1, 3, 10000) * u.micron).to(u.AA)
        norm = np.trapz(self.sensitivity(wavelength), wavelength)
        return norm * self.zeropoint * 10 ** (-mag / 2.5)

    def mag_to_average_flux_density(self, mag):
        """Convert magnitude to average flux density based on the zeropoint of the detector"""
        if not isinstance(mag, u.Quantity):
            mag = u.Quantity(mag, u.dimensionless_unscaled)
        if mag.unit != u.dimensionless_unscaled:
            raise ValueError("Magnitude must have dimensionless units.")
        wavelength = (np.linspace(0.1, 3, 10000) * u.micron).to(u.AA)
        norm = np.trapz(self.sensitivity(wavelength), wavelength)
        return norm * self.zeropoint * 10 ** (-mag / 2.5)
