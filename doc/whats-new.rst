.. currentmodule:: pywatershed

What's New
==========

.. ipython:: python
    :suppress:

    import numpy as np
    import pywatershed as pws

    np.random.seed(123456)


.. _whats-new.0.1.2:

v0.1.2 (unreleased)
-----------------------

New Features
~~~~~~~~~~~~

- Add automatic release workflow to PyPi as per :issue:`178`. Associated
  implementation of gitflow and semver conventions. Overhauled
  `CONTRIBUTING.md`, `DEVELOPER.md`, `README.md`, and `.github/RELEASE.md`
  to document adopted practices. Adoption of `git-cliff` to generate change
  logs by filtering comitt messages, see `cliff.toml`. Clean up of environment
  files and streamlining against `pyproject.toml`. Symlink gfortran dylibs to
  `/usr/local/lib` on macOS CI so PRMS binaries included in this repo can find
  them. (:pull:`179`)
  By `Wes Bonelli <https://github.com/w-bonelli>`_.
- What's new workflow (behold!) per :issue:`180` and :pull:`181`
  By `James McCreight <https://github.com/jmccreight>`_.


Breaking changes
~~~~~~~~~~~~~~~~


Deprecations
~~~~~~~~~~~~

Performance
~~~~~~~~~~~

- Introduce ASV performance benchmarks for import and various NHM configurations
  in pywatershed. (:issue:`170` and :pull:`184`)
  By `James McCreight <https://github.com/jmccreight>`_.


Bug fixes
~~~~~~~~~


Documentation
~~~~~~~~~~~~~


Internal Changes
~~~~~~~~~~~~~~~~

.. _whats-new.0.1.1:

v0.1.1 (27 April 2023)
--------------------

Initial release.
