#!/bin/bash

source_dir=/Users/jamesmcc/usgs/data/prms_nhm_applications/drb
target_dir=/Users/jamesmcc/usgs/data/prms_nhm_applications/drb_2yr

run_dir=/Users/jamesmcc/usgs/data/prms_nhm_applications/drb_2yr_check
source_run_dir=$run_dir/run_drb
target_run_dir=$run_dir/run_drb_2yr
rm -rf $source_run_dir
rm -rf $source_target_dir
mkdir -p $source_run_dir
mkdir -p $target_run_dir

prms=/Users/jamesmcc/usgs/prmsNHMpy/prms5.2.1/bin/prms


# The OG
cd $source_run_dir
cp $prms .
for ii in control.default myparam.param sf_data *cbh
do
    cp $source_dir/$ii .
done
# change the start and end dates by line number
# change the output option
# osx specific sed
sed -i '' "6s/.*/1979/" control.default
sed -i '' "7s/.*/1/" control.default
sed -i '' "16s/.*/1981/" control.default
sed -i '' "17s/.*/1/" control.default
sed -i '' "136s/.*/2/" control.default
./prms control.default -set param_file ./myparam.param


# The Subset
cd $target_run_dir
cp $prms .
for ii in control.default myparam.param sf_data *cbh
do
    cp $source_dir/$ii .
done
# change the start and end dates by line number
# change the output option
# osx specific sed
sed -i '' "6s/.*/1979/" control.default
sed -i '' "7s/.*/1/" control.default
sed -i '' "16s/.*/1981/" control.default
sed -i '' "17s/.*/1/" control.default
sed -i '' "136s/.*/2/" control.default
./prms control.default -set param_file ./myparam.param

result=`\diff $source_run_dir/stats.csv $target_run_dir/stats.csv`

if [[ $result -eq 0 ]]; then
    echo 'Success'
else
    echo 'Failure'
fi

exit $result
