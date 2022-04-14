# test_data/

This directory contains
  * domain directories for domain tests
  * scripts/ which document how domain tests are established.


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



# scripts/

This contains code for reproducing test data in the domains. Importantly, `run_domains.py` should be run
occasionally to update the test data.

Original source data can be found on Denali for most (all?) domains in `/home/jmccreight/pynhm_data`.
