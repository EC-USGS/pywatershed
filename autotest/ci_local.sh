#!/bin/bash

# This is a local version of CI testing.
# Unfortunately it has to be kept up to date with ci.yaml.

# Notes:
# * This is developed only on my Mac M1, so there are implicit assumptions
#   around that which may not work on other platforms or on M1 machines
#   which are not set up identically.



# local configuration
pytest_n=8
# should probably clone mf6 locally and checkout latest develop
modflow_repo_location=../../modflow6_for_pws_ci

# options
# all "no data" options. if passed, these turn OFF sections of the tests.
while getopts 'hilmtosrdug' opt; do
  case "$opt" in
    h)
      h=h
      echo "Printing HELP:"
      ;;
    i)
      i=i
      echo "Not testing or re-installing pywatershed"
      ;;
    l)
      l=l
      echo "Not linting pywatershed"
      ;;
    m)
      m=m
      echo "Not updating or building Modflow6"
      ;;
    t)
      t=t
      echo "Not running the tests"
      ;;
    o)
      o=o
      echo "Not running the domainless tests"
      ;;
    s)
      s=s
      echo "Not running the sagehen_5yr tests"
      ;;
    r)
      r=r
      echo "Not running the hru_1 tests"
      ;;
    d)
      d=d
      echo "Not running the drb_2yr tests"
      ;;
    u)
      u=u
      echo "Not running the ucb_2yr tests"
      ;;
    g)
      g=g
      echo "Not generating test data for any run tests"
      ;;
  esac
done
shift "$(($OPTIND -1))"


if [ ! -z "${h}" ]; then
    echo "Using '-' infront of a letter turns that section off"
    echo
    echo "Structure of options"
    echo "--------------------"
    echo "i: installation"
    echo "l: linting"
    echo "m: modflow update and build"
    echo "t: tests"
    echo "  o: domainless tests"
    echo "  s: sagehen_5yr"
    echo "    g: generate sagehen data"
    echo "  r: hru_1"
    echo "    g: generate hru_1 data"
    echo "  d: drb_2yr"
    echo "    g: generate drb_2yr data"
    echo "  u: ucb_2yr"
    echo "    g: generate ucb_2yr data"

    exit 0
fi
    
echo ""
echo ""


start_dir=`pwd`

# name: Set environment variables
export PYWS_FORTRAN=false
export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"
export PYNHM_FORTRAN=false
export `head -n1 ../.mf6_ci_ref_remote`
export `tail -n1 ../.mf6_ci_ref_remote`


if [ -z "${i}" ]; then
    echo
    echo
    echo "******************************"
    echo "Installation"
    echo "******************************"
    echo

    # run from repository root
    cd ..

    pip uninstall -y pywatershed || exit 1

    ## name: Upgrade pip and install build and twine
    python -m pip install --upgrade pip || exit 1
    pip install wheel build 'twine<5.0.0' 'importlib_metadata<=7.0.1' || exit 1

    ## name: Base installation
    pip --verbose install . || exit 1

    ## name: Print pyhmn version
    python -c "import pywatershed; print(pywatershed.__version__)" || exit 1

    ## name: Build pywatershed, check dist outputs
    python -m build || exit 1
    twine check --strict dist/* || exit 1

    cd $start_dir || exit 1
fi


if [ -z "${l}" ]; then
    echo
    echo
    echo "******************************"
    echo "Linting: check and format"
    echo "******************************"
    echo

    # run from repository root
    cd ..
    ruff check . || exit 1
    ruff format --check . || exit 1
    cd $start_dir || exit 1
fi


if [ -z "${m}" ]; then
    echo
    echo
    echo "******************************"
    echo "Modflow6 Update and Build"
    echo "******************************"
    echo

    # name: Enforce MF6 ref and remote merge to main
    req_ref=develop  # if not develop, submit an issue
    echo $MF6_REF
    if [[ "$MF6_REF" != "$req_ref" ]]; then exit 1; fi
    req_remote=MODFLOW-USGS/modflow6
    echo $MF6_REMOTE
    if [[ "$MF6_REMOTE" != "$req_remote" ]]; then
	echo "bad mf6 remote in .mf6_ci_ref_remote"
	exit 1
    fi

    # Checkout MODFLOW 6 (from $start_dir)
    if [ ! -d $modflow_repo_location ]; then
	git clone git@github.com:$req_remote $modflow_repo_location || exit 1
    fi
    cd "${modflow_repo_location}" || exit 1
    git checkout $req_ref || exit 1
    git fetch origin || exit 1
    git merge origin/$req_ref || exit 1

    # Update flopy MODFLOW 6 classes in the current environment
    cd autotest || exit 1
    python update_flopy.py


    # Build mf6 locally instead of installing mf6 nightly build
    # install conda env for mf6
    cd "${modflow_repo_location}" || exit 1
    env_name=mf64ci
    # only necessary the first time ?
    # env_file=environment.yml
    # mamba remove -y --name $env_name --all || exit 1
    # mamba create -y --name $env_name || exit 1
    # mamba env update --name $env_name --file $env_file --prune  || exit 1

    # conda activate $env_name
    source /Users/jmccreight/mambaforge/bin/activate $env_name
    # only necessary the first time
    # meson setup --prefix=$(pwd) --libdir=bin builddir
    meson install -C builddir
    conda deactivate

    cd $start_dir
fi

export PATH=$PATH:$modflow_repo_location/bin

# Use the installation above if performed, else use an existing installation
# - name: Install pywatershed
#   run: |
#     pip install .

# - name: Version info
#   run: |
#     pip -V
#     pip list

if [ -z "${t}" ]; then
    echo
    echo
    echo "******************************"
    echo "TESTS"
    echo "******************************"
    echo

   cd ..

   echo
   echo "Get GIS files for tests"
   python pywatershed/utils/gis_files.py || exit 1

   cd autotest

   if [ -z "${o}" ]; then
       echo
       echo
       echo "===================="
       echo "DOMAINLESS"
       echo "===================="
       echo
       echo "domainless - run tests not requiring domain data"
       pytest -m domainless -n=$pytest_n -vv || exit 1
   fi

   if [ -z "${s}" ]; then
       echo
       echo
       echo "===================="
       echo "DOMAIN: sagehen_5yr"
       echo "===================="
       echo
       if [ -z "${g}" ]; then
	   echo
	   echo ".........."
	   echo "sagehen_5yr_no_cascades - generate and manage test data domain, "
	   echo "  run PRMS and convert csv output to NetCDF"
	   python generate_test_data.py \
		  -n=$pytest_n --domain=sagehen_5yr \
		  --control_pattern=sagehen_no_cascades.control \
		  --remove_prms_csvs --remove_prms_output_dirs || exit 1
       fi

       # - name: sagehen_5yr_no_cascades - list netcdf input files
       #   working-directory: test_data
       #   run: |
       #     find sagehen_5yr/output_no_cascades -name '*.nc'

       echo
       echo ".........."
       echo "sagehen_5yr_no_cascades - pywatershed tests"
       echo ".........."
       echo
       pytest \
	   -vv \
	   -rs \
	   -n=$pytest_n \
	   -m "not domainless" \
	   --domain=sagehen_5yr \
	   --control_pattern=sagehen_no_cascades.control \
	   --durations=0  || exit 1
   fi

   if [ -z "${r}" ]; then
       echo
       echo
       echo "===================="
       echo "DOMAIN: hru_1"
       echo "===================="
       echo
       if [ -z "${g}" ]; then

	   echo
	   echo ".........."
	   echo "hru_1_nhm - generate and manage test data domain, run PRMS "
	   echo "  and convert csv output to NetCDF"
	   echo ".........."
	   echo
	   python generate_test_data.py \
	       -n=$pytest_n --domain=hru_1 --control_pattern=nhm.control \
	       --remove_prms_csvs --remove_prms_output_dirs || exit 1
       fi

       # - name: hru_1_nhm - list netcdf input files
       #   working-directory: test_data
       #   run: |
       #     find hru_1/output -name '*.nc'

       echo
       echo ".........."
       echo "hru_1_nhm - pywatershed tests"
       echo ".........."
       echo
       pytest \
	   -vv \
	   -rs \
	   -n=$pytest_n \
	   -m "not domainless" \
	   --domain=hru_1 \
	   --control_pattern=nhm.control \
	   --durations=0 || exit 1

   fi

   if [ -z "${d}" ]; then
       echo
       echo
       echo "===================="
       echo "DOMAIN: drb_2yr"
       echo "===================="
       echo
       if [ -z "${g}" ]; then
	   echo
	   echo ".........."
	   echo "drb_2yr with and without dprst and obsin - generate and "
	   echo "  manage test data"
	   echo ".........."
	   echo
	   python generate_test_data.py \
	      -n=$pytest_n --domain=drb_2yr \
	      --remove_prms_csvs --remove_prms_output_dirs || exit 1
       fi

       # - name: drb_2yr_nhm - list netcdf input files
       #   working-directory: test_data
       #   run: |
       #     find drb_2yr/output -name '*.nc'

       # - name: drb_2yr_no_dprst - list netcdf input files
       #   working-directory: test_data
       #   run: |
       #     find drb_2yr/output_no_dprst -name '*.nc'

       echo ".........."
       echo "drb_2yr_nhm - pywatershed tests"
       echo ".........."
       echo
       echo "Running:"
       echo "    pytest -vv -rs -n=$pytest_n -m 'not domainless' --domain=drb_2yr "
       echo "       --control_pattern=nhm.control   --durations=0"
       pytest \
           -vv \
           -rs \
           -n=$pytest_n \
           -m "not domainless" \
           --domain=drb_2yr \
           --control_pattern=nhm.control \
           --durations=0 || exit 1

       # Specific tests not redundant with dprst
       echo ".........."
       echo "drb_2yr_no_dprst - pywatershed tests"
       echo ".........."
       echo
       pytest \
           test_prms_runoff.py \
           test_prms_soilzone.py \
           test_prms_groundwater.py \
           test_prms_above_snow.py \
           test_prms_below_snow.py \
           -vv \
           -rs \
           -n=$pytest_n \
           -m "not domainless" \
           --domain=drb_2yr \
           --control_pattern=no_dprst \
           --durations=0 || exit 1

       # # Specific tests not redundant with dprst
       echo ".........."
       echo "drb_2yr_obsin - pywatershed tests"
       echo ".........."
       echo
       pytest \
           test_obsin_flow_node.py \
           -vv \
           -n=0 \
           -m "not domainless" \
           --domain=drb_2yr \
           --control_pattern=nhm_obsin.control \
           --durations=0 || exit 1

   fi

   if [ -z "${u}" ]; then
       echo
       echo
       echo "===================="
       echo "DOMAIN: ucb_2yr"
       echo "===================="
       echo
       if [ -z "${g}" ]; then
	   echo
	   echo ".........."
	   echo "ucb_2yr_nhm - generate and manage test data"
	   echo ".........."
	   echo
	   python generate_test_data.py \
	      -n=$pytest_n --domain=ucb_2yr --control_pattern=nhm.control \
	      --remove_prms_csvs --remove_prms_output_dirs || exit 1
       fi

       # - name: ucb_2yr_nhm - list netcdf input files
       #   working-directory: test_data
       #   run: |
       #     find ucb_2yr/output -name '*.nc'

       echo ".........."
       echo "ucb_2yr_nhm - pywatershed tests"
       echo ".........."
       echo
       pytest \
           -vv \
	   -n=$pytest_n \
           -m "not domainless" \
           --domain=ucb_2yr \
           --control_pattern=nhm.control \
           --durations=0 || exit 1
   fi
       
fi

if [ -z "${i}" ]; then
    # If install was done, put the install back to its original, editable state.
    # Do it here so the tests use the test install if it is done.
    pip uninstall -y pywatershed || exit 1
    cd .. || exit 1
    pip install -e . || exit 1

    cd $start_dir || exit 1
fi

exit 0
