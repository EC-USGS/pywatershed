# pynhm
[![ci-badge](https://github.com/ec-usgs/pynhm/workflows/CI/badge.svg?branch=main)](https://github.com/ec-usgs/pynhm/actions?query=workflow%3ACI)
[![codecov-badge](https://codecov.io/gh/ec-usgs/pynhm/branch/main/graph/badge.svg)](https://codecov.io/gh/ec-usgs/pynhm)
[![Documentation Status](https://readthedocs.org/projects/pynhm/badge/?version=latest)](https://pynhm.readthedocs.io/en/latest/?badge=latest)
<img src="https://raw.githubusercontent.com/ec-usgs/pynhm/main/resources/images/prms_flow.png" alt="prms_flow" style="width:50;height:20">

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


Requirements and Python Environments
====================================
Git and python 3.8 - 3.10 are required. 

We suggest installing the python dependencies using conda and the following yaml files. 

* Minimal: ci/requirements/environment.yml
* Documentation: ci/requirements/doc.yml
* Notebooks: examples/examples_env.yml

```conda env create -f env_of_choice.yml```

More detailed installation instructions using conda can be found in `examples/00_python_virtual_env.ipynb`.

There are also .txt equivalents that can be used for installing from pip, like so:

```pip install -r env_of_choice.txt```

though these are not comprehensive installs as with conda. 


Compiled Code: Fortran
=============
The numpy extension F2PY is used to provide compiled versions of core calculation codes to boost performance. F2PY is documented (within numpy)[https://numpy.org/doc/stable/f2py/index.html]. This repository is configured to compile on install. This may not always be successful and may depend on your local environment. Currently, compilation often requires

`export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"`

on MacOS. Please report compilation issues.

To turn off compilation of fortran source on install, prior to `pip install .` set the environment variable

`export PYNHM_FORTRAN=False`

where the value on the right side is not case sensitive. This may break tests of fortran source, but will allow the rest of the package to be installed.

Disclaimer
==========

This information is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The information has not received final approval by the U.S. Geological Survey (USGS) and is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the information.

From: https://www2.usgs.gov/fsp/fsp_disclaimers.asp#5

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits "official USGS copyright policy")

Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.
This software is provided "AS IS."

