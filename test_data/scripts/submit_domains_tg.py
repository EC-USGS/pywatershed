#!/bin/bash
#SBATCH --job-name=conus_2yr
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=1
#SBATCH --account=impd
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jmccreight@usgs.gov
#SBATCH -o conus_tg_job--%j.out

# This is for TALLGRASS

# conda envs are documented in /home/jmccreight/conda_installs.sh

# make sure we have conda in the scheduler
source /home/jmccreight/.bashrc

#  of gcc/gfort that PRMS 5.2.1 requires
conda activate gnu_tg
which gcc
which gfortran

conda activate --stack pynhm_tg
conda list

python run_domains.py

exit $?
