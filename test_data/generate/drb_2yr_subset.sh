#!/bin/bash
# Purpose:
# Generate a 2 year test case of DRB from the full 40yr

data_dir=../../../../data
source_dir=$data_dir/prms_nhm_applications/drb
target_dir=$data_dir/prms_nhm_applications/drb_2yr

mkdir -p $target_dir

invariant_files=(control.default  myparam.param)
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

chmod -R 555 $target_dir

exit 0
