# pywatershed

[![ci-badge](https://github.com/ec-usgs/pywatershed/workflows/CI/badge.svg?branch=develop)](https://github.com/ec-usgs/pywatershed/actions?query=workflow%3ACI)
[![codecov-badge](https://codecov.io/gh/ec-usgs/pywatershed/branch/main/graph/badge.svg)](https://codecov.io/gh/ec-usgs/pywatershed)
[![Documentation Status](https://readthedocs.org/projects/pywatershed/badge/?version=latest)](https://pywatershed.readthedocs.io/en/latest/?badge=latest)
[![asv](http://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat)](https://github.com/ec-usgs/pywatershed)
[![Formatted with black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)

[![Available on pypi](https://img.shields.io/pypi/v/pywatershed.svg)](https://pypi.python.org/pypi/pywatershed)
[![PyPI Status](https://img.shields.io/pypi/status/pywatershed.svg)](https://pypi.python.org/pypi/pywatershed)
[![PyPI Versions](https://img.shields.io/pypi/pyversions/pywatershed.svg)](https://pypi.python.org/pypi/pywatershed)

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pywatershed/badges/version.svg)](https://anaconda.org/conda-forge/pywatershed)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pywatershed/badges/platforms.svg)](https://anaconda.org/conda-forge/pywatershed)

[![WholeTale](https://raw.githubusercontent.com/whole-tale/wt-design-docs/master/badges/wholetale-explore.svg)](https://dashboard.wholetale.org/run/64ae29e8a887f48b9f173678?tab=metadata)


<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Purpose](#purpose)
- [Installation](#installation)
- [Contributing](#contributing)
- [Example Notebooks](#example-notebooks)
- [Overview of Repository Contents](#overview-of-repository-contents)
- [Disclaimer](#disclaimer)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Purpose

The purpose of this repository is to refactor and redesign the [PRMS modeling
system](https://www.usgs.gov/software/precipitation-runoff-modeling-system-prms)
while maintaining its functionality. Code modernization is a step towards
unification with [MODFLOW 6 (MF6)](https://github.com/MODFLOW-USGS/modflow6).

The following motivations are taken from our [AGU poster from December
2022](https://agu2022fallmeeting-agu.ipostersessions.com/default.aspx?s=05-E1-C6-40-DF-0D-4D-C7-4E-DE-D2-61-02-05-8F-0A)
which provides additional details on motivations, project status, and current
directions of this project as of approximately January 2023.

Goals of the USGS Enterprise Capacity (EC) project include:

  * A sustainable integrated, hydrologic modeling framework for the U.S.
    Geological Survey (USGS)
  * Interoperable modeling across the USGS, partner agencies, and academia

Goals for EC Watershed Modeling:

  * Couple the Precipitation-Runoff Modeling System (PRMS, e.g. Regan et al,
	2018)  with MODFLOW 6 (MF6, e.g. Langevin et al, 2017) in a sustainable
	way
  * Redesign PRMS to be more modern and flexible
  * Prioritize process representations in the current National Hydrological
    Model (NHM) based on PRMS 5.2.1

Prototype an EC watershed model: "pywatershed"

  * Redesign PRMS quickly in python
  * Couple to MF6 via BMI/XMI interface (Hughes et al, 2021; Hutton et al, 2020)
  * Establish a prototyping ground for EC codes that couples to the compiled
	framework: low cost proof of concepts (at the price of potentially less
    computational performance) * Enable process representation hypothesis testing
  * Use cutting-edge techniques and technologies to improve models 
  * Machine learning, automatic differentiation 
  * Address challenges of modeling across space and time scales 
  * Transition prototype watershed model to compiled EC code


## Installation

`pywatershed` uses Python 3.9 or 3.10.

The `pywatershed` package is [available on
PyPI](https://pypi.org/project/pywatershed/) but installation of all
dependencies sets (lint, test, optional, doc, and all) may not be reliable on
all platforms. 

The `pywatershed` package is [available on
conda-forge](https://anaconda.org/conda-forge/pywatershed). The installation
is the quickest way to get up and running by provides only the minimal set of
dependencies (not including jupyter nor all packages needed for running the
example notebooks, also not suitable for development purposes). 

We recommend the following installation procedures to get fully-functional
environments for running `pywatershed` and its example notebooks. We strongly
recommend using [Mamba](https://mamba.readthedocs.io/en/latest/)to first
instal dependencies from the `environment_y_jupyter.yml` file in the
repository before installing `pywatershed` itself. Mamba will be much faster
than Ananconda (but the conda command could also be used). 

If you wish to use the stable release, you will use `main` in place of 
`<branch>` in the following commands. If you want to follow developemnt, you'll
use `develop` instead.

Without using `git` (directly), you may:
```
curl -L -O https://raw.githubusercontent.com/EC-USGS/pywatershed/<branch>/environment_w_jupyter.yml
mamba env create -f environment_w_jupyter.yml
conda activate pws
pip install git+https://github.com/EC-USGS/pywatershed.git@<branch>
```

Or to use `git` and to be able to develop:

```
git clone https://github.com/EC-USGS/pywatershed.git
cd pywatershed
mamba env create -f environment_w_jupyter.yml
activate pws
pip install -e .
```

(If you want to name the environment other than the default `pws`, use the
command 
`mamba env update --name your_env_name --file environment_w_jupyter.yml --prune`
you will also need to activate this environment by name.)


We install the `environment_w_jupyter.yml` to provide all known dependencies 
including those for running the eample notebooks. (The `environment.yml` 
does not contain jupyter or jupyterlab because this interferes with installation
on WholeTale, see Example Notebooks seection below.)

## Contributing

See the [developer documentation](./DEVELOPER.md) for instructions on setting up
a development environment. See the [contribution guide](./CONTRIBUTING.md) to
contribute to this project.

## Example Notebooks

For introductory example notebooks, look in the
[`examples/`](https://github.com/EC-USGS/pywatershed/tree/main/examples>)
directory in the repository. Numbered starting at 00, these are meant to be
completed in order. Non-numbered notebooks coveradditional topics. These
notebooks are note yet covered by testing and so may be expected to have some
issues until they are added to testing. In `examples/developer/` there are
notebooks of interest to developers who may want to learn about running the
software tests.

Though no notebook outputs are saved in Github, these notebooks can easily
navigated to and run in WholeTale containers (free but sign-up or log-in
required). This is a very easy and quick way to get started without needing to
install pywatershed requirements yourself. WholeTale is an NSF funded project
and supports logins from many institutions, e.g. the USGS, and you may not need
to register.

There are containers for both the `main` and `develop` branches.

[![WholeTale](https://raw.githubusercontent.com/whole-tale/wt-design-docs/master/badges/wholetale-explore.svg)](https://dashboard.wholetale.org)

  * [WholeTale container for latest release (main
	branch)](https://dashboard.wholetale.org/run/64ae29e8a887f48b9f173678?tab=metadata)
  * [WholeTale container for develop
	branch](https://dashboard.wholetale.org/run/64ae25c3a887f48b9f1735c8?tab=metadata)

WholeTale will give you a jupyter-lab running in the root of this
repository. You can navigate to `examples/` and then open and run the notebooks
of your choice.  The develop container may require the user to update the
repository (`git pull origin`) to stay current with development.

## Overview of Repository Contents

The contents of directories at this level is described. Therein you may discover
another README.md for more information.

```
.github/: Github actions, scripts and Python environments for continuous integration (CI) and releasing,
asv_benchmarks/: preformance benchmarking by ASV
autotest/: pywatershed package testing using pytest
autotest_exs/: pywatershed example notebook testing using pytest
bin/:PRMS executables distributed
doc/:Package/code documentation source code
evaluation/: tools for evaluation of pywatershed
examples/:How to use the package, mostly jupyter notebooks
prms_src/:PRMS source used for generating executables in bin/
pywatershed/:Package source
reference/:Ancillary materials for development
resources/:Static stuff like images
test_data/:Data used for automated testing
```

## Disclaimer

This information is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The information has not
received final approval by the U.S. Geological Survey (USGS) and is provided on
the condition that neither the USGS nor the U.S. Government shall be held liable
for any damages resulting from the authorized or unauthorized use of the
information.

From: https://www2.usgs.gov/fsp/fsp_disclaimers.asp#5

This software is in the public domain because it contains materials that
originally came from the U.S. Geological Survey, an agency of the United States
Department of Interior. For more information, see the [official USGS copyright
policy](https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits
"official USGS copyright policy")

Although this software program has been used by the USGS, no warranty, expressed
or implied, is made by the USGS or the U.S. Government as to the accuracy and
functioning of the program and related program material nor shall the fact of
distribution constitute any such warranty, and no responsibility is assumed by
the USGS in connection therewith.  This software is provided "AS IS."
