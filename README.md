# pywatershed

[![ci-badge](https://github.com/ec-usgs/pywatershed/workflows/CI/badge.svg?branch=develop)](https://github.com/ec-usgs/pywatershed/actions?query=workflow%3ACI)
[![codecov-badge](https://codecov.io/gh/ec-usgs/pywatershed/branch/main/graph/badge.svg)](https://codecov.io/gh/ec-usgs/pywatershed)
[![Documentation Status](https://readthedocs.org/projects/pywatershed/badge/?version=latest)](https://pywatershed.readthedocs.io/en/latest/?badge=latest)
[![asv](http://img.shields.io/badge/benchmarked%20by-asv-blue.svg?style=flat)](https://github.com/ec-usgs/pywatershed)

[![PyPI Version](https://img.shields.io/pypi/v/pywatershed.png)](https://pypi.python.org/pypi/pywatershed)
[![PyPI Status](https://img.shields.io/pypi/status/pywatershed.png)](https://pypi.python.org/pypi/pywatershed)
[![PyPI Versions](https://img.shields.io/pypi/pyversions/pywatershed.png)](https://pypi.python.org/pypi/pywatershed)

[//]: # (<img src="https://raw.githubusercontent.com/ec-usgs/pywatershed/main/resources/images/prms_flow.png" alt="prms_flow" style="width:50;height:20">)

## Purpose

The purpose of this repository is to refactor and redesign the [PRMS modeling
system](https://www.usgs.gov/software/precipitation-runoff-modeling-system-prms)
while maintaining its functionality. Code modernization is a step towards
unification with [MODFLOW 6 (MF6)](https://github.com/MODFLOW-USGS/modflow6).

The following motivations are taken from our
[AGU poster from December 2022](https://agu2022fallmeeting-agu.ipostersessions.com/default.aspx?s=05-E1-C6-40-DF-0D-4D-C7-4E-DE-D2-61-02-05-8F-0A)
which provides additional details on motivations, project status, and current
directions of this project as of approximately January 2023.

Goals of the USGS Enterprise Capacity (EC) project include:
  * A sustainable integrated, hydrologic modeling framework for the U.S. Geological Survey (USGS)
  * Interoperable modeling across the USGS, partner agencies, and academia

Goals for EC Watershed Modeling:
  * Couple the Precipitation-Runoff Modeling System (PRMS, e.g. Regan et al, 2018)  with MODFLOW 6 (MF6, e.g. Langevin et al, 2017) in a sustainable way
  * Redesign PRMS to be more modern and flexible
  * Prioritize process representations in the current National Hydrological Model (NHM) based on PRMS 5.2.1

Prototype an EC watershed model: "pywatershed"
  * Redesign PRMS quickly in python
  * Couple to MF6 via BMI/XMI interface (Hughes et al, 2021; Hutton et al, 2020)
  * Establish a prototyping ground for EC codes that couples to the compiled framework: low cost proof of concepts (at the price of potentially less computational performance)
  * Enable process representation hypothesis testing
  * Use cutting-edge techniques and technologies to improve models
  * Machine learning, automatic differentiation
  * Address challenges of modeling across space and time scales
  * Transition prototype watershed model to compiled EC code

## Installation

To install the software you will need Python 3.9 or 3.10.

We currently recommend dependencies be installed with [Mamba](https://mamba.readthedocs.io/en/latest/)
which will be much faster than Ananconda. An environment containing all core and optional dependencies
can be created from the project root with:

```
conda env create -f environment.yml
```

The `pywatershed` package is [available on PyPI](https://pypi.org/project/pywatershed/) and can be installed with:

```
pip install pywatershed
```

A number of extra dependencies are needed to run the example notebooks. These can be installed with pip with

```
pip install "pywatershed[optional]"
```

These installation steps are suitable for `pywatershed` end users. See the [developer documentation](./DEVELOPER.md) for detailed instructions on configuring a development environment.

## Contributing

See the [developer documentation](./DEVELOPER.md) for instructions on setting up a development environment. See the [contribution guide](./CONTRIBUTING.md) to contribute to this project.

## Example Notebooks

Jupyter notebooks containing examples are found in the
[examples/](https://github.com/EC-USGS/pywatershed/tree/main/examples) directory. Numbered notebooks are tested.
Notebooks 00 and 01 walk the user through the setting the python environment and running the software tests.
Notebook 02 demonstrates modeling with pywatershed. Non-numbered notebooks cover additional topics. These notebooks
are note yet covered by testing and so may be expected to have some issues until they are added to testing.

## Overview of Repository Contents

The contents of directories at this level is described. Therein you may discover another README.md for more information.

```
.github/    Github actions, scripts and Python environments for continuous integration (CI) and releasing, 
autotest/   pywatershed package testing using pytest
bin/        PRMS executables distributed
doc/        Package/code documentation source code
examples/   How to use the package, mostly jupyter notebooks
prms_src/   PRMS source used for generating executables in bin/
pywatershed/      Package source
reference/  Ancillary materials for development
resources/  Static stuff like images
test_data/  Data used for automated testing
```

## Disclaimer

This information is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The information has not received final approval by the U.S. Geological Survey (USGS) and is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the information.

From: https://www2.usgs.gov/fsp/fsp_disclaimers.asp#5

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits "official USGS copyright policy")

Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.
This software is provided "AS IS."
