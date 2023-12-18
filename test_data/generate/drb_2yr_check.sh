#!/bin/bash

# config
data_dir=../../../../data
source_dir=$data_dir/prms_nhm_applications/drb
target_dir=$data_dir/prms_nhm_applications/drb_2yr
run_dir=$data_dir/prms_nhm_applications/drb_2yr_check
repo_dir=../../../
prms=$repo_dirprms_src/prms5.2.1/bin/prms

# ---------------------------------

source_run_dir=$run_dir/run_drb
target_run_dir=$run_dir/run_drb_2yr
# The control file for this is under version control
control_file=$repo_dir/test_data/drb_2yr/control.test

rm -rf $source_run_dir
rm -rf $source_target_dir
mkdir -p $source_run_dir/output
mkdir -p $target_run_dir/output


# The OG
cd $source_run_dir
cp $prms .
for ii in myparam.param sf_data *cbh
do
    cp $source_dir/$ii .
done
cp $control_file .
./prms control.test -set param_file ./myparam.param || exit 2


# The Subset
cd $target_run_dir
cp $prms .
for ii in myparam.param sf_data *cbh
do
    cp $target_dir/$ii .
done
cp $control_file .
./prms control.test -set param_file ./myparam.param || exit 3

# Checks
diff_files () {
    file1=$1
    file2=$2
    if [[ $file1 == $file2 ]];
    then
	echo
	echo "Diffing the same file: $file1 == $file2"
	echo
	return 1
    fi
    ls $file1 > /dev/null 2>&1 || return 1
    ls $file2 > /dev/null 2>&1 || return 1
    \diff $file1 $file2
    result=$?
    if [[ $result -eq 0 ]]
    then
	echo success: `basename $file1`
    else
	echo FAILURE: `basename $file1`
    fi
    return $result
}

for ff in $target_run_dir/output/*
do
    ff2=$source_run_dir/output/$(basename $ff)
    diff_files $ff $ff2
    result=$?
    if [[ $result -ne 0 ]]
    then
	exit 11
    fi
done

controlled_files=(
    control.test \
    myparam.param \
    prcp.cbh \
    rhavg.cbh \
    tmax.cbh \
    tmin.cbh \
    output/nhru_tmaxf.csv \
    output/nhru_tminf.csv \
    output/nhru_hru_ppt.csv \
    output/nhru_hru_rain.csv \
    output/nhru_hru_snow.csv \
    output/soltab_debug )

echo "You can optionally copy the version controlled files "
echo "to the repo using the following commands: "
echo
for ff in "${controlled_files[@]}"
do
    echo "cp $target_run_dir/$ff $repo_dir/test_data/drb_2yr/$ff"
done
echo


exit 0
