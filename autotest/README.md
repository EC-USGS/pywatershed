# Autotest

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
  --print_ans=PRINT_ANS
                        Print results and assert False for all domain tests
```

Note that the `--domain_yaml` argument may be repeated within a single call to test multiple 
domains in a single pytest.

The domain_yaml file is an *evolving* set of data that includes static domain
inputs (e.g. CBH forcing files, parameter files), static or reference model 
output (from PRMS/NHM), and the answers to domain tests. The answers for domain
tests can be printed using the `--print_ans` argument which prints the test
answers (up to user to make sure these are correct) and makes the test fail. 
This option aids in conctruction of the tests_ans section of the domain
yaml file. 




