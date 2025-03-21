# pywatershed

[![ci-badge](https://github.com/ec-usgs/pywatershed/workflows/CI/badge.svg?branch=develop)](https://github.com/ec-usgs/pywatershed/actions?query=workflow%3ACI)
[![codecov-badge](https://codecov.io/gh/ec-usgs/pywatershed/branch/main/graph/badge.svg)](https://codecov.io/gh/ec-usgs/pywatershed)
[![Documentation Status](https://readthedocs.org/projects/pywatershed/badge/?version=latest)](https://pywatershed.readthedocs.io/en/latest/?badge=latest)
[![asv](http://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat)](https://github.com/ec-usgs/pywatershed)
[![Formatted with Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![Available on pypi](https://img.shields.io/pypi/v/pywatershed.svg)](https://pypi.python.org/pypi/pywatershed)
[![PyPI Status](https://img.shields.io/pypi/status/pywatershed.svg)](https://pypi.python.org/pypi/pywatershed)
[![PyPI Versions](https://img.shields.io/pypi/pyversions/pywatershed.svg)](https://pypi.python.org/pypi/pywatershed)

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pywatershed/badges/version.svg)](https://anaconda.org/conda-forge/pywatershed)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pywatershed/badges/platforms.svg)](https://anaconda.org/conda-forge/pywatershed)

[![DOI:10.5066/P13EWPEV](https://img.shields.io/badge/DOI-10.5066/P13EWPEV-b4a9fe.svg)](https://doi.org/10.5066/P13EWPEV)

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [About](#about)
- [Installation](#installation)
- [Getting started / Example notebooks](#getting-started--example-notebooks)
- [Community engagement](#community-engagement)
- [How to Cite](#how-to-cite)
- [Disclaimer](#disclaimer)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## About

Welcome to the pywatershed repository!

Pywatershed is Python package for simulating hydrologic processes motivated by
the need to modernize important, legacy hydrologic models at the USGS,
particularly the
[Precipitation-Runoff Modeling System](https://www.usgs.gov/software/precipitation-runoff-modeling-system-prms)
(PRMS, Markstrom et al., 2015) and its role in
[GSFLOW](https://www.usgs.gov/software/gsflow-coupled-groundwater-and-surface-water-flow-model>)
(Markstrom et al., 2008).
The goal of modernization is to make these legacy models more flexible as process
representations, to support testing of alternative hydrologic process
conceptualizations, and to facilitate the incorporation of cutting edge
modeling techniques and data sources. Pywatershed is a place for experimentation
with software design, process representation, and data fusion in the context
of well-established hydrologic process modeling.

For more information on the goals and status of pywatershed, please see the [pywatershed docs](https://pywatershed.readthedocs.io/).

## Installation

`pywatershed` uses Python 3.10 or 3.11.

The `pywatershed` package is [available on
PyPI](https://pypi.org/project/pywatershed/) but installation of all
dependencies sets (lint, test, optional, doc, and all) may not be reliable on
all platforms.

The `pywatershed` package is [available on
conda-forge](https://anaconda.org/conda-forge/pywatershed). The installation
is the quickest way to get up and running by provides only the minimal set of
dependencies (not including Jupyter nor all packages needed for running the
example notebooks, also not suitable for development purposes).

We recommend the following installation procedures to get fully-functional
environments for running `pywatershed` and its example notebooks. We strongly
recommend using [Mamba](https://mamba.readthedocs.io/en/latest/)to first
instal dependencies from the `environment_y_jupyter.yml` file in the
repository before installing `pywatershed` itself. Mamba will be much faster
than Ananconda (but the conda command could also be used).

If you wish to use the stable release, you will use `main` in place of
`<branch>` in the following commands. If you want to follow development, you'll
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
including those for running the example notebooks.

## Getting started / Example notebooks

Please note that you can browse the API reference, developer info, and index
in the [pywatershed docs](<(https://pywatershed.readthedocs.io/)>). But
_the best way to get started with pywatershed is to dive into the example
notebooks_.

For introductory example notebooks, look in the
[`examples/`](https://github.com/EC-USGS/pywatershed/tree/main/examples>)
directory in the repository. Numbered starting at 00, these are meant to be
completed in order. Numbered starting at 00, these are meant to be completed
in order. Notebook outputs are not saved in Github.

Non-numbered notebooks in `examples/` cover additional topics. These
notebooks are not yet covered by testing and you may encounter some
issues. In `examples/developer/` there are notebooks of interest to
developers who may want to learn about running the software tests.

## Community engagement

We value your feedback! Please use [discussions](https://github.com/EC-USGS/pywatershed/discussions)
or [issues](https://github.com/EC-USGS/pywatershed/issues) on Github.
For more in-depth contributions, please start by reading over
the pywatershed
[DEVELOPER.md](https://github.com/EC-USGS/pywatershed/blob/develop/DEVELOPER.md) and
[CONTRIBUTING.md](https://github.com/EC-USGS/pywatershed/blob/develop/CONTRIBUTING.md)
guidelines.

Thank you for your interest.

## How to Cite

McCreight, J. L., Langevin, C. D., Hughes, J. D., & Bonelli, W. P. (2024). pywatershed (Version 2.0.0) [Computer software]. [https://doi.org/10.5066/P13EWPEV](https://doi.org/10.5066/P13EWPEV)

## Disclaimer

This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of
release constitute any such warranty. The software is provided on the condition
that neither the USGS nor the U.S. Government shall be held liable for any
damages resulting from the authorized or unauthorized use of the software.
