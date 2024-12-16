<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [ASV Benchmarking](#asv-benchmarking)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# ASV Benchmarking

[https://github.com/airspeed-velocity/asv/](https://github.com/airspeed-velocity/asv/)
[https://asv.readthedocs.io/en/stable/](https://asv.readthedocs.io/en/stable/)

__Notes__

This suite of performance benchmarks tests import and various aspects of
running the nhm configuration(s) in pywatershed.

The benchmark suite runs on local hardware but could be implemented in CI
(though consistent performance can be tricky in CI). To run on your machine,
follow the readthedocs link above and register your hardware.

The .asv driectory is not particularly sacred. A typical command for doing
performance regression is

ASV apparently copies your repo (including .git) and checks stuff out.
Commits requested have to be on the branch specified in this repo (not
just on any branch).

```
asv continuous --verbose --show-stderr --factor 1.5 53a2e51929  dd99eb7922
```

Currently, [ASV is not getting the version of python correct](https://github.com/airspeed-velocity/asv/issues/1294)
in the environment it is building (with conda). I managed to get it to use
python 3.10 for the `./environment.yml` file. When I tried reducing
dependencies needed for benchmarking, it kept installing python 3.11 which wont
install pywatershed per the `pyproject.toml` requires-python field. For this
reason, in `__import__` i'm printing `sys.version` and since each test is run
in a separate subprocess, the actual python version number will be littered
around your screen when using the `--verbose` flag.

