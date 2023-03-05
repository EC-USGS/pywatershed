# Test data

The `test_data` directory contains
  * domain directories for domain tests
  * `scripts/`, which contains scripts to run domain simulations and generate data.

## Domain directories

Some domain directories only contain `control.test` and `domain.yaml` files. These domains 
are not to be run in a self-contained way within the repository. They rely on other files on 
other systems. The main example is `conus_2yr/`. It is to be run on Denali/Tallgrass but its
input files are too large to be included in this repo and it is desirable to write model 
output to a disk faster than the one where the repository lives. However, we do want to 
version control certain files in the directory. See the `conus_2yr/README.md` for details 
on how it is setup.

Other domain directories can be run from pytest given the repository. The pytest 
flag `--all_domains` will detect these repos and run them.

Original source data can be found on Denali for most (all?) domains in `/home/jmccreight/pywatershed_data`.

# Generating data

The `test_data/scripts` subdirectory contains code for reproducing test data in the domains. Importantly, `test_run_domains.py` should be run occasionally to update the test data.  After running the domains, NetCDF files can be created from simulation outputs by running the tests in `test_nc_domains.py`. E.g.,

```shell
pytest -v -n auto test_run_domains.py
pytest -v -n auto test_nc_domains.py
```

NetCDF dependencies are encoded implicitly into the `pytest` fixture system: `test_nc_domains.py` uses a custom test parametrization with `pytest_generate_tests` to map each CSV file created by the domain simulation to one or more NetCDF files, which are then aggregated into further files on session teardown by [yield fixtures](https://docs.pytest.org/en/7.2.x/how-to/fixtures.html#teardown-cleanup-aka-fixture-finalization). A [filelock](https://pytest-xdist.readthedocs.io/en/latest/how-to.html#making-session-scoped-fixtures-execute-only-once) is used to ensure aggregate files are only created once, even with multiple `pytest-xdist` workers.
