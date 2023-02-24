# pynhm
[![ci-badge](https://github.com/ec-usgs/pynhm/workflows/CI/badge.svg?branch=main)](https://github.com/ec-usgs/pynhm/actions?query=workflow%3ACI)
[![codecov-badge](https://codecov.io/gh/ec-usgs/pynhm/branch/main/graph/badge.svg)](https://codecov.io/gh/ec-usgs/pynhm)
[![Documentation Status](https://readthedocs.org/projects/pynhm/badge/?version=latest)](https://pynhm.readthedocs.io/en/latest/?badge=latest)

[//]: # (<img src="https://raw.githubusercontent.com/ec-usgs/pynhm/main/resources/images/prms_flow.png" alt="prms_flow" style="width:50;height:20">)

Purpose
=========
The purpose of this repository is to refactor and redesign the PRMS modeling system to modernize it and to move towards a unification with
MODFLOW 6.

The following motivations are taken from our [AGU poster from December 2022](https://agu2022fallmeeting-agu.ipostersessions.com/default.aspx?s=05-E1-C6-40-DF-0D-4D-C7-4E-DE-D2-61-02-05-8F-0A)
which provides additional details on motivations, project status, and current directions of this project as of approximately January 2023.


Project goals of the USGS Enterprise Capacity (EC) include:
  * A sustainable integrated, hydrologic modeling framework for the U.S. Geological Survey (USGS)
  * Interoperable modeling across the USGS, partner agencies, and academia

Goals for EC Watershed Modeling:
  * Couple the Precipitation-Runoff Modeling System (PRMS, e.g. Regan et al, 2018)  with MODFLOW 6 (MF6, e.g. Langevin et al, 2017) in a sustainable way
  * Redesign PRMS to be more modern and flexible
  * Prioritize process representations in the current National Hydrological Model (NHM) based on PRMS 5.2.1

Prototype an EC watershed model: "pynhm"
  * Redesign PRMS quickly in python
  * Couple to MF6 via BMI/XMI interface (Hughes et al, 2021; Hutton et al, 2020)
  * Establish a prototyping ground for EC codes that couples to the compiled framework: low cost proof of concepts (at the price of potentially less computational performance)
  * Enable process representation hypothesis testing
  * Use cutting-edge techniques and technologies to improve models
  * Machine learning, automatic differentiation
  * Adress challenges of modeling across space and time scales
  * Transition prototype watershed model to compiled EC code


Overview
==========
The contents of directories at this level is described. Therein you may discover another README.md for more information.

```
.github/    Github actions for deploying continuous integration (CI)
autotest/   pynhm package testing using pytest
bin/        PRMS executables distributed
ci/         Python environments for CI
doc/        Package/code documentation source code
examples/   How to use the package, mostly jupyter notebooks
prms_src/   PRMS source used for generating executables in bin/
pynhm/      Package source
reference/  Ancillary materials for development
resources/  Static stuff like images
test_data/  Data used for automated testing
```

Installation and Python Environments
=====================================
To install the software (without running tests) you will need *Python >= 3.8*.

We recommend installing the python package dependencies using anaconda or miniconda. Generally
most users will likely want to create the `pynhm_nb` conda environment by running.

```conda env create -f examples/examples_env.yml```

One could also do

```pip install -r examples/exampes_env.txt```

but this is not guaranteed.

Once the environment is established,

`cd pynhn; pip install .`

will install the package. If you would like to compile the fortran compuational kernels for
certain physical process representations, you'll need a fortran compiler and you will run

`export PYNHM_FORTRAN=true; cd pynhm;  pip install .`

See Developer Requirements below for more details.


Developer Requirements and Python Environments
====================================
*Git* is helpful (but not strictly required).
Both C and Fortran compilers are required. We are currently using gnu (gcc, gfortran) 11 and 12
as well as intel (icc, ifort) 2021 on Windows, Linux, and MacOS (including Apple Silicon). Both of
these are freely obtainable but the installation process varies widely. We are looking for
a conda-based approach to this problem, but currently do not have a solution. Generally *gcc* and
*gfortran* are recommended. On *Apple Silicon*, the PRMS source code  currently
only compiles with intel while the fortran kernels in pynhm only compile with gnu.

*Python >= 3.8* is required. Three different python environments are specified within the repository.
These are:

* Minimal (for developing/testing), 'pynhm': ci/requirements/environment.yml
* Notebooks (~= minimal + jupyter), 'pynhm_nb': examples/examples_env.yml
* Documentation (only if you want to build the documentation), 'pynhm-docs': ci/requirements/doc.yml

We recommend (because we test it in CI) using anacoda or miniconda to establish these environments
with the following commands

```conda env create -f env_of_choice.yml```

which will create the environment with "name" specified on the first line of the file, given before the path
to the file above.

More detailed installation instructions using conda can be found in `examples/00_python_virtual_env.ipynb`.

There are also .txt equivalents that can be used for installing from pip, like so:

```pip install -r env_of_choice.txt```

though these are not comprehensive installs as with conda.

Once the python envionment and dependencies are established, pynhm is installed for development with the
following command

`cd pynhn; pip install -e .`

The numpy extension F2PY is used to provide fortran compiled kernels of core calculations to boost
performance. F2PY is documented (within numpy)[https://numpy.org/doc/stable/f2py/index.html]. This
repository is configured to compile on install. This may not always be successful and may depend
on your local environment. Currently, we have not established this compilation procedure for Windows.
On linux and MacOS, compilation of fortran kernels on package installation is achieved by the following
code:

```
export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"
export CC=path/to/gcc  # for example
export FC=path/to/gfortran  # for example
export PYNHM_FORTRAN=true
cd path/to/pynhm
pip install -e .
```

To run the tests of the repository, we first need to generate the test data. This consists of running PRMS
and then converting the output to netcdf:

```
cd path/to/pynhm/test_data/scripts
pytest -v -n=4 test_run_domains.py
pytest -v -n=8 test_nc_domains.py
```

Finally, run the tests themselves,

```
cd path/to/pynhm/autotest
pytest -v -n=8
```

All tests should pass, XPASS, or XFAIL. XFAIL is an expected failure.

Contributing
============
We welcome community development! Please file Issues and/or Pull Requests in the appropriate places on github. The continuous
integration (CI) procedure is the first gate keeper for new code contribution. The CI procdure is defined by
`.github/workflows/ci.yaml`. This includes running the formatting and linting packages `isort`, `black`, and
`flake8` in addition to generating the test data and running the tests in `autotest/`. New codes need new tests so they can
be verified moving ahead in time.

Disclaimer
==========

This information is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The information has not received final approval by the U.S. Geological Survey (USGS) and is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the information.

From: https://www2.usgs.gov/fsp/fsp_disclaimers.asp#5

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits "official USGS copyright policy")

Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.
This software is provided "AS IS."
