<a href="https://github.com/pandoramission/pandora-sat/actions/workflows/tests.yml"><img src="https://github.com/pandoramission/pandora-sat/workflows/pytest/badge.svg" alt="Test status"/></a> [![Generic badge](https://img.shields.io/badge/documentation-live-blue.svg)](https://pandoramission.github.io/pandora-sat/)

# PandoraSat

This Python package contains metadata for Pandora and **basic** functions describing things such as detector sensitivity and zeropoint estimation.

### Installation

Eventually, to install this package you will be able to use

```
pip install pandora-sat --upgrade
```

For now, clone this repository, enter the `pandora-sat` directory, and execute the command `poetry install` to build the package locally. Make sure to upgrade regularly to get the latest estimates of the Pandora meta data.


### Example Usage

Below is an example usage of some of the functionality in this package. In general, this package will allow you to get metadata from specific subsystems of Pandora.

```python
from pandorasat import PandoraSat
p = PandoraSat()
print(p.NIRDA.pixel_scale)
print(p.Hardware.mirror_diameter)
print(p.VISDA.sensitivity(wavelength))
print(p.Orbit.period)
```

See our API documentation for full details on the metadata available in this package.

To update any of the values or functions contained within `pandora-sat` due to new testing, commissioning, etc., please open a pull request. Update the relevant values or functions on your branch and then the updates will be reviewed prior to being merged into the main branch of `pandora-sat`.
