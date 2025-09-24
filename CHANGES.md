# 0.12.5

- Updated `phoenix.py` context manager so that the vega spectrum is always set right in the config

# 0.12.0

- Updated `pandorasat` to use `pandoraref` as the source of reference data products.
- Moved location of some SED tools to the `phoenix` file
- Added benchmark file back in

# 0.11.0

- @yoavrotman Updated the pysynphot dependency to synphot and stsynphot

# 0.10.0

- Updated the WCS so no reflections in either axis are expected. This is our current expectation, but could change.

# 0.9.2

- Updated throughput to full rolled up expectation

# 0.8.6

- updated `VisibleDetector.get_wcs` to accept optional crpix arguments

# 0.8.5

- updated `get_sky_catalog` to take an epoch argument.

# 0.8.0

- added new sources of noise
- updated magnitude to flux functions
- added documentation
- added tutorial
- updated transmission functions
- added configuration file
- changed phoneix grid location
- added benchmark SED to package files

# 0.5.0

- completed separation of pandora-sat from pandora-sim

# 0.4.1

- Separated all simulation tools into a separate package, pandora-sim
- Renamed optics.py to hardware.py to make its contents more general
- Wrote new tests to reflect the new intent of the package as a container for information

# 0.3

- Changed API to enforce row-major indexing everywhere
- Added `conventions.ipynb` documentation to state conventions

# 0.2.0

- Added ability to create a distorted WCS
- Added jitter in the observatory rotation angle as well as position
- Increased field stop radius
- Added variable gain
- Added a very simple way to "simulate" cosmic rays

# 0.1.2

- Added capability and documentation to make basic visual simulations
