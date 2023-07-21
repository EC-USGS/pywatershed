.. currentmodule:: pywatershed

What's New
==========

.. ipython:: python
    :suppress:

    import numpy as np
    import pywatershed as pws

    np.random.seed(123456)


.. _whats-new.0.2.1:

v0.2.1 (19 July 2023)
---------------------

Bug fixes
~~~~~~~~~
- Package data was not properly installed.
  (:pull:`219`) By `James McCreight <https://github.com/jmccreight>`_.
- Small addition to notebook 02
  (:pull:`219`) By `James McCreight <https://github.com/jmccreight>`_.


.. _whats-new.0.2.0:

v0.2.0 (18 July 2023)
---------------------

New features
~~~~~~~~~~~~
- New example notebooks. Moved old notebooks to `examples/developer`.
  (:pull:`204`)
  By `James McCreight <https://github.com/jmccreight>`_.
- New way to specify `Model` instantiation either in-memory or from yaml files
  using a model dictionary. The approach is loosely based on MODFLOW 6 input
  organization. See `Model` documentation. Introduced the concept of
  discretizations for PRMS, defining "dis_hru" and "dis_seg". These are
  components of how model dictionaries are specified.
  (:pull:`188`) By `James McCreight <https://github.com/jmccreight>`_.
- New `Control.from_yaml()` method. (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.
- What's new workflow (behold!) per :issue:`180` and :pull:`181`
  By `James McCreight <https://github.com/jmccreight>`_.
- Add automatic release workflow to PyPi as per :issue:`178`. Associated
  implementation of gitflow and semver conventions. Overhauled
  `CONTRIBUTING.md`, `DEVELOPER.md`, `README.md`, and `.github/RELEASE.md`
  to document adopted practices. Adoption of `git-cliff` to generate change
  logs by filtering comitt messages, see `cliff.toml`. Clean up of environment
  files and streamlining against `pyproject.toml`. Symlink gfortran dylibs to
  `/usr/local/lib` on macOS CI so PRMS binaries included in this repo can find
  them. (:pull:`179`)
  By `Wes Bonelli <https://github.com/w-bonelli>`_.


Breaking changes
~~~~~~~~~~~~~~~~
- Move Control attribute "config" to "options" for handling global options.
  (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.
- Remove arguments from `Model` initialization. Options pass via control, new
  `set_options()` method on Process and ConservativeProcess
  (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.
- `Control` no longer takes a `Parameter` object as an initialization argument.
  `Process` subclasses now require arguments `discretization` and `parameters`.
  The firstargument of `Model` not a indefinite number of processes, it is now
  either a list of `Process` subclasses or a model dictionary (see `Model`
  documentation. (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.


Deprecations
~~~~~~~~~~~~


Performance
~~~~~~~~~~~
- Introduce ASV performance benchmarks for import and various NHM configurations
  in pywatershed. (:issue:`170` and :pull:`184`)
  By `James McCreight <https://github.com/jmccreight>`_.


Bug fixes
~~~~~~~~~
- Remove non-pep-compliant post-release reset PR steps. (:pull:`203`)
  By `Wes Bonelli <https://github.com/w-bonelli>`_.
- Add doc building requirements to environment.yml (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.
- Revive fortran compiling for editable installs (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.
- Made the Parameter class data completely private by converting dicts to
  MappingProxyTypes and setting numpy.ndarrays to read-only. (:issue:`177`
  and :pull:`183`)
  By `James McCreight <https://github.com/jmccreight>`_.
- ModelGraph improvements and fixes (however result are platform dependent)
  (:pull:`162`) By `James McCreight <https://github.com/jmccreight>`_.


Documentation
~~~~~~~~~~~~~
- Model class, DatasetDict and general documentation overhaul (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.


Internal changes
~~~~~~~~~~~~~~~~
- Introduce precommit hooks: `nbstripout`, `blackdoc`, and `doctoc`.
  (:pull:`197`)
  By `James McCreight <https://github.com/jmccreight>`_.
- Rename StorageUnit to ConservativeProcess that subclasses from a new Process
  class that contains most of the StorageUnit functionality. (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.
- New set_options() method on Process and ConservativeProcess to set their
  initialization options as '_` atrributes. (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.
- Clean up of how the `calc_method` option assigns function names to reduce
  the total amount of code and do it upon initialization. (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.
- Rename many modules to use lower-snake-case names including those in base/,
  atmoshpere/, and hydrology/ (:pull:`188`)
  By `James McCreight <https://github.com/jmccreight>`_.
- NHM "self-driving" tests
  (:pull:`160`)
  By `James McCreight <https://github.com/jmccreight>`_.
- Refactor dependencies for standard pypi installation. (:pull:`164`,
  :issue:`178`)
  By `Joseph Hughes <https://github.com/jdhughes-usgs>`_.


.. _whats-new.0.1.1:

v0.1.1 (27 April 2023)
----------------------

Initial release.
