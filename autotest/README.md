# Autotest

## Usage

```
cd autotest
pytest
```
Pytest options can be explored via `pytest --help`.


## Developer

This is how the pynhm package tests itself.

The test suite consists both
	* stand alone tests, and
	* domain tests

Stand alone tests do not require any input files or output files besides what
is supplied in the testing framework. These tests generally run quickly
(e.g. testing basic logic of a class or type requirements or results).

Domain tests require input or output files for an NHM domain. These tests
scale with the size or number of HRUs in a domain. The test suite takes
arguments related to domain tests:

```
 --domain_yaml=DOMAIN_YAML
                        YAML file(s) for indiv domain tests. You can pass multiples of this argument. Default
                        value (not shown here) is --domain_yaml=../test_data/drb_2yr/drb_2yr.yaml.
  --print_ans           Print results and assert False for all domain tests
```

Note that the `--domain_yaml` argument may be repeated within a single call to test multiple
domains in a single pytest.

The `domain_yaml` file is an *evolving* set of data that includes static domain
inputs (e.g. CBH forcing files, parameter files), static or reference model
output (from PRMS/NHM), and the answers to domain tests.

Examples of `domain_yaml` files can be found in, for example, in
`pynhm/test_data/drb_2yr/drb_2yr.yaml`
and
`pynhm/test_data/conus_2yr/conus_2yr.yaml`.


### Domain inputs
These are specific files for a specific model domain (e.g. the Delaware River Basin,
or the CONUS NHM, etc). The files are paths and are to be specified as either
  * relative paths: relative to the location of the yaml file
    in which the path is appearing
  * absolute paths: use absolute paths (for some machine and potentially user)
The examples listed above demonstrate use of both relative and absolute paths.

Domain inputs are generally created by scripts run in `test_data/scripts`. The
`drb_2yr` case shows the relationship between the prms binary, the input files,
and the output files. These scripts help maintain sanity when generating new
files for domain tests.


### Answers for domain tests
Tests have objectively correct answers. The answer key is stored in the yaml in
top-level attribute/key called `test_ans`. This can be seen in the examples
listed above. This attibute has a sub-attributes for each file with domain
tests in `autotest/test_*py` (the text `test_` and `.py` is dropped in the
yaml reference). Below this, each file has named tests attributes which in turn
potentially have named cases/keys when iterated with other fixtures (such as
variables or types).

The tests answers for a given test are output of some kind of summary statistic
or other kind of reduction performed on data in memory. The answers vary with
domain and are most easily collected by running the tests themselves while
verifying accuracy manually when putting the answers into the domain yaml file.

For convenience, a user can select to print the answers at run/test time using
the `--print_ans` argument when running ptest. This prints the answer for
domain tests and always asserts False after doing so for each test in a file.
This option aids in conctruction of the tests_ans section of the domain
yaml file for new domains or when extending tests on existing domains. Running
`pytest -s ...` prints output in a format that can be copied to the yaml and
edited quickly to get new test results enshrined in the answer key.
