#!/bin/bash
# Purpose:
# Generate a 2 year test case of DRB from the full 40yr

home_dir=~

source_dir=$home_dir/pywatershed_data/20220209_gm_byHWobs_CONUS
target_dir=$home_dir/pywatershed_data/conus_2yr
pywatershed_repo_dir=$home_dir/pywatershed

mkdir -p $target_dir

# we'll take the control file from the repo
cp $pywatershed_repo_dir/test_data/conus_2yr/control.test $target_dir/.

invariant_files=(myparam.param)
cbh_files=(prcp.cbh  rhavg.cbh  tmax.cbh  tmin.cbh)

# (2 years daily) + (1980 is leap year) + (three header lines)
n_lines=$(expr 365 \* 2 + 1 + 3)

for ii in "${invariant_files[@]}"
do
    echo $ii
    cp $source_dir/$ii $target_dir/$ii
done

for cc in "${cbh_files[@]}"
do
    echo $cc
    head -n $n_lines $source_dir/$cc > $target_dir/$cc
done

# write protect the 2yr subset
chmod -R 555 $target_dir

exit 0
