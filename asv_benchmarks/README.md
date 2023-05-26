# ASV Benchmarking

[https://github.com/airspeed-velocity/asv/](https://github.com/airspeed-velocity/asv/)
[https://asv.readthedocs.io/en/stable/](https://asv.readthedocs.io/en/stable/)

This suite of performance benchmarks tests import and various aspects of
running the nhm configuration(s) in pywatershed.

Currently, ASV is not getting the version of python correct in the environment
it is building (with conda). I managed to get it to use python 3.10 for the
`./environment.yml` file. When I tried reducing dependencies needed for
benchmarking, it kept installing python 3.11 which wont install pywatershed
per the `pyproject.toml` requires-python field.

The .asv driectory is not particularly sacred. A typical command for doing
performance regression is

```
asv continuous --verbose --show-stderr --factor 1.5 53a2e51929  dd99eb7922
```

ASV apparently copies your repo (including .git) and checks stuff out.
Commits requested have to be on the branch specified in this repo (not
just on any branch).

