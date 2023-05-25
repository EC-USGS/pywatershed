# Developing `pywatershed`

This document describes how to set up a development environment, install dependencies, compile C/Fortran code, generate example data, and run the tests and example notebooks.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Requirements](#requirements)
  - [Git](#git)
  - [Compilers](#compilers)
  - [Python](#python)
    - [Creating a virtual environment](#creating-a-virtual-environment)
      - [Conda](#conda)
      - [Pip](#pip)
    - [Installing `pywatershed` in development mode](#installing-pywatershed-in-development-mode)
    - [F2PY](#f2py)
- [Testing](#testing)
- [Branching model](#branching-model)
- [Miscellaneous](#miscellaneous)
  - [Locating the root](#locating-the-root)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Requirements

### Git

[Git](https://git-scm.com/) is required to contribute code back to the repository. GitHub's [Guide to Installing Git](https://help.github.com/articles/set-up-git) is a good source of information.

### Compilers

C and Fortran compilers are required. We are currently using gnu (gcc, gfortran)
11 and 12 as well as intel (icc, ifort) 2021 on Windows, Linux, and MacOS
(including Apple Silicon). Both of these are freely obtainable but the installation
process varies widely. We are looking for a conda-based approach to obtaining
compilers, but currently do not have a solution. Compilers are needed for
two applications:

  1. Compiling and running C/Fortran PRMS code to generate testing/verification data
  2. Compiling (installing) and running fortran backends/kernels for some hydrological
     process representations in pywatershed

On Apple Silicon, the PRMS source code is only currently known to compile with intel while
the fortran kernels in pywatershed only compile with gnu.

### Python

Python >= 3.9 is required.

#### Creating a virtual environment

A virtual Python environment is recommended, and can be created with Python's builtin `venv` or a conda distribution like Anaconda or Miniconda.

##### Conda

To create a conda environment with all core, testing and optional dependencies, run from the project root:

```
conda env create -f environment.yml
```

The conda `environment.yml` contains all dependencies needed to develop, test, and run example notebooks. More detailed Conda environment installation instructions can be found in `examples/00_python_virtual_env.ipynb`.

##### Pip

To install all dependencies with `pip`:

```
pip install ".[all]"
```

Several  dependency groups are defined in `pyproject.toml` and can be selected instead of `all` for a more lightweight environment:

- `lint`
- `test`
- `optional`

#### Installing `pywatershed` in development mode

Once the python environment and dependencies are established and activated (e.g. `conda activate pywatershed` or, assuming a virtual environment in `venv` in the root, `source venv/bin/activate`), `pywatershed` can be installed in [editable development mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html) with:

```
pip install -e .
```

#### F2PY

The numpy extension F2PY is used to provide fortran compiled kernels of core calculations to boost
performance. F2PY is documented [within numpy](https://numpy.org/doc/stable/f2py/index.html). This
repository is configured NOT to compile on install by default. Currently, we have not established
this compilation procedure for Windows. On Linux and MacOS, compilation of fortran kernels on package
installation is achieved by setting several environent variables before installing the `pywatershed` module.
For instance, from the project root:

```
export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"
export CC=path/to/gcc  # for example
export FC=path/to/gfortran  # for example
export PYWS_FORTRAN=true
pip install -e .
```

## Testing

Once the dependencies are available, we want to verify the software by running its test suite. The
following testing procedures are also covered in the notebook `examples/01_automated_testing.ipynb`.
To run the tests, we first need to generate the test data. This consists of running PRMS
and then converting the output to netcdf. From `test_data/scripts`:

```
pytest -v -n=4 test_run_domains.py
pytest -v -n=8 test_nc_domains.py
```

Finally, the tests can be run from the `autotest` directory:

```
pytest -v -n=8
```

All tests should pass, XPASS, or XFAIL. XFAIL is an expected failure. Substitute `-n auto` to automatically use all available cores on your machine.

## Branching model

This project uses the [git flow](https://nvie.com/posts/a-successful-git-branching-model/): development occurs on the `develop` branch, while `main` is reserved for the state of the latest release. Development PRs are typically squashed to `develop`, to avoid merge commits. At release time, release branches are merged to `main`, and then `main` is merged back into `develop`.

## Miscellaneous

### Locating the root

Python scripts often need to reference files elsewhere in the project. To allow scripts to be run from anywhere in the project hierarchy, scripts should locate the project root relative to themselves, then use paths relative to the root for file access, rather than using relative paths (e.g., `../some/path`).

For a script in a subdirectory of the root, for instance, the conventional approach would be:

```Python
project_root_path = Path(__file__).parent.parent
```
