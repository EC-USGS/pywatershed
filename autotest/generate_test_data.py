#!/usr/bin/env python
"""Generate pywatershed test data.

usage examples:
    python:  `python generate_test_data.py --domains=hru_1,drb_2yr -n5`
    ipython: `ipython generate_test_data.py -- -n=5 --domains=hru_1,drb_2yr`
    command line: `./generate_test_data.py --domains=hru_1,drb_2yr -n5`

For help, call with -h or --help. Note that arguments not listed by
generate_test_data.py are passed directly to pytest. See `pytest --help`
for its extensive list of options.

For details look at DEVELOPER.md, autotest/README.md, and test_data/README.md.

This wrapper calls pytest to create and destroy test data in the test_data
directory. It also creates and destroys dot files which are looked for and at
by autotest to get a rough sense if the test data are up to date. These dot
files are `test_data/.test_data_version_{domain_name}.txt` for each domain.

When in doubt about the freshness of test data, remake it!

Note that this wrapper was create because the test version file
tracking was practically impossible with pytest-xdist.
"""

import argparse
import os
import pathlib as pl
import shutil
import sys
from copy import deepcopy
from pprint import pprint

import pytest


def parse_args():
    """Parse the arguments.

    Args: None.
    Returns: A list of arguments parsed in a way to pass to run_prms_domains.py
        and convert_prms_output_to_nc.py.

    Special handling of --domains here to kind of cleanup inconsistencies
    betweenrun_prms_domains.py and convert_prms_output_to_nc.py and also to
    control recording of which domains have available test_data.
    """
    desc = (
        "Generate pywatershed test data. "
        "*** Arguments not listed by generate_test_data.py are passed "
        "directly to pytest. See `pytest --help` for its extensive list of "
        "options. ***"
    )
    parser = argparse.ArgumentParser(
        description=desc,
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="Help on calling generate_test_data.py.",
    )
    parser.add_argument(
        "--domains",
        default="hru_1,drb_2yr,ucb_2yr",
        help=(
            "Domain name values for which to generate test data, comma "
            "separate multiple values."
        ),
        type=str,
    )
    parser.add_argument(
        "--remove_prms_csvs",
        help=(
            "Option to remove PRMS csv output in domain/output/ directories "
            "after running"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--remove_prms_output_dirs",
        help=(
            "Option to remove existing PRMS output directories before running"
        ),
        action="store_true",
    )

    known, unknown = parser.parse_known_args()

    domains = known.domains
    if "," in domains:
        domains_split = domains.split(",")
        domain_list = []
        for dd in domains_split:
            domain_list += [f"--domain={dd}"]
    else:
        domain_list = [f"--domain={domains}"]

    remove_prms_csvs = [f"--remove_prms_csvs={known.remove_prms_csvs}"]
    remove_prms_output_dirs = [
        f"--remove_prms_output_dirs={known.remove_prms_output_dirs}"
    ]

    arg_list = (
        domain_list + unknown + remove_prms_csvs + remove_prms_output_dirs
    )

    return arg_list


def main():
    arg_list = parse_args()

    print("\nGenerating test data. Running in ../test_data/generate")
    os.chdir("../test_data/generate")

    # collect info on domains to be run_domains
    domains_to_run = []
    control_patterns = []
    remove_prms_csvs = False
    remove_prms_output_dirs = False
    rm_vals = []
    for ii in arg_list:
        if "--domain" in ii:
            domains_to_run += [ii.split("--domain")[-1].split("=")[-1]]
        if "--control_pattern" in ii:
            control_patterns += [
                ii.split("--control_pattern")[-1].split("=")[-1]
            ]
        elif "--remove_prms_csvs" in ii:
            remove_prms_csvs = ii.split("=")[-1] == "True"
            rm_vals += [ii]
        elif "--remove_prms_output_dirs" in ii:
            remove_prms_output_dirs = ii.split("=")[-1] == "True"
            rm_vals += [ii]

    for ii in rm_vals:
        arg_list.remove(ii)

    simulations = []
    for domain in domains_to_run:
        dom_dir = pl.Path(f"../{domain}")
        candidates = sorted(dom_dir.glob("*.control"))

        if not len(control_patterns):
            simulations += candidates
            continue
        else:
            for candidate in candidates:
                for pattern in control_patterns:
                    if pattern in candidate.name:
                        simulations += [str(candidate)]

    n_simulations = len(simulations)

    print("\nGenerating test data for the following domain/control_files:")
    pprint(simulations)
    print()

    for ss in simulations:
        ss_pl = pl.Path(ss)
        domain_tag = ss_pl.parent.name
        control_tag = ss_pl.with_suffix("").name
        test_data_version_file = pl.Path(
            f"../.test_data_version_{domain_tag}_{control_tag}.txt"
        )
        test_data_version_file.unlink(missing_ok=True)

    if remove_prms_output_dirs:
        # todo: make this take the domain argument
        print("\nRemoving existing PRMS output directories ...")
        rm_output_dirs_list = arg_list + ["remove_output_dirs.py"]
        retcode_output_dirs = pytest.main(rm_output_dirs_list)
        if retcode_output_dirs != pytest.ExitCode.OK:
            print("\nRemoving existing PRMS output dirs failed.")
            sys.exit(retcode_output_dirs.value)

    print("\nRunning PRMS domains ...")
    # if -n is in the list and numeric and > ndomains+1
    #     then set it to ndomains + 1 just for running the domains
    n_orig = None
    for ii in arg_list:
        if "-n" in ii:
            arg_list_no_n = deepcopy(arg_list)
            _ = arg_list_no_n.remove(ii)
            n_orig = ii.split("-n")[-1].split("=")[-1]
            n_orig = int(n_orig) if n_orig.isdigit() else n_orig

    if (
        (n_orig is None)
        or (not isinstance(n_orig, int))
        or (n_orig <= (n_simulations + 1))
    ):
        run_arg_list = arg_list + ["run_prms_domains.py"]
        conv_arg_list = arg_list + ["convert_prms_output_to_nc.py"]
    else:
        run_arg_list = arg_list_no_n + [
            f"-n={n_simulations + 1}",
            "run_prms_domains.py",
        ]
        conv_arg_list = arg_list_no_n + [
            f"-n={n_orig}",
            "convert_prms_output_to_nc.py",
        ]

    retcode_run = pytest.main(run_arg_list)
    if retcode_run != pytest.ExitCode.OK:
        print("Running PRMS domains failed.")
        sys.exit(retcode_run.value)

    print("\nConverting PRMS output to NetCDF ...")
    retcode_conv = pytest.main(conv_arg_list)
    if retcode_conv != pytest.ExitCode.OK:
        print("\nConverting PRMS output to NetCDF failed.")
        sys.exit(retcode_run.value)

    for ss in simulations:
        ss_pl = pl.Path(ss)
        domain_tag = ss_pl.parent.name
        control_tag = ss_pl.with_suffix("").name
        test_data_version_file = pl.Path(
            f"../.test_data_version_{domain_tag}_{control_tag}.txt"
        )
        shutil.copy("../../version.txt", test_data_version_file)

    if remove_prms_csvs:
        # todo: make this take the domain argument
        print("\nRemoving PRMS CSV output file ...")
        rm_csv_arg_list = arg_list + ["remove_prms_csvs.py"]
        retcode_rm_csvs = pytest.main(rm_csv_arg_list)
        if retcode_rm_csvs != pytest.ExitCode.OK:
            print("\nRemoving PRMS CSV output files failed.")
            sys.exit(retcode_rm_csvs.value)

    print(
        "\nSuccess generating pywatershed test data for domains "
        f"{domains_to_run}\n"
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
