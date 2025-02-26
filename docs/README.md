<a href="https://github.com/pandoramission/pandora-sat/actions/workflows/tests.yml"><img src="https://github.com/pandoramission/pandora-sat/workflows/tests/badge.svg" alt="Test status"/></a><a href="https://github.com/pandoramission/pandora-sat/actions/workflows/black.yml"><img src="https://github.com/pandoramission/pandora-sat/workflows/black/badge.svg" alt="black status"/></a> <a href="https://github.com/pandoramission/pandora-sat/actions/workflows/flake8.yml"><img src="https://github.com/pandoramission/pandora-sat/workflows/flake8/badge.svg" alt="flake8 status"/></a> [![Generic badge](https://img.shields.io/badge/documentation-live-blue.svg)](https://pandoramission.github.io/pandora-sat/)
[![PyPI - Version](https://img.shields.io/pypi/v/pandorasat)](https://pypi.org/project/pandorasat/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pandorasat)](https://pypi.org/project/pandorasat/)

# PandoraSat

This Python package contains metadata for Pandora and **basic** functions describing things such as detector sensitivity and zeropoint estimation.

### Installation

To install you can use

```
pip install pandorasat --upgrade
```

You should update your package often, as we frequently put out new versions with updated Current Best Estimates, and some limited new functionality. Check your version number using

```
import pandorasat as ps
ps.__version__
```

## Pandora Information

This repository helps you understand what the Pandora SmallSat will be capabale of.

## Example Usage

Below is an example usage of some of the functionality in this package. In general, this package will allow you to get metadata from specific subsystems of Pandora.

```python
from pandorasat import VisibleDetector, NIRDetector
visda = VisibleDetector()
visda.pixel_scale
visda.pixel_size

nirda = NIRDetector()
nirda.plot_sensitivity()
```

See our API documentation for full details on the metadata available in this package.

To update any of the values or functions contained within `pandora-sat` due to new testing, commissioning, etc., please open a pull request. Update the relevant values or functions on your branch and then the updates will be reviewed prior to being merged into the main branch of `pandora-sat`.
