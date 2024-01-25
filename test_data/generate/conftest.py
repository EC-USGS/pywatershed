import os
import pathlib as pl
import sys
from fnmatch import fnmatch
from platform import processor
from shutil import copy2
import time
from typing import List
from warnings import warn

import pytest

rootdir = ".."

# Subset this to speed up tests by eliminating some domains
test_dirs = sorted(
    [path for path in pl.Path(rootdir).iterdir() if path.is_dir()]
)

# This would change to handle other/additional schedulers
domain_globs_schedule = ["*conus*"]

final_var_names = ["through_rain", "seg_lateral_inflow"]


def pytest_addoption(parser):
    parser.addoption(
        "--force",
        required=False,
        action="store_true",
        default=False,
        help=("Force run without scheduler present"),
    )

    parser.addoption(
        "--domain",
        required=False,
        action="append",
        default=[],
        help=(
            "Domain(s) to run (name of domain dir and NOT path to it). "
            "You can pass multiples of this argument. If not used, "
            "all self-contained domains in test_data/ will be run."
        ),
    )


@pytest.fixture()
def exe():
    platform = sys.platform.lower()
    if platform == "win32":
        exe_name = "prms_win_gfort_dbl_prec.exe"
    elif platform == "darwin":
        if processor() == "arm":
            exe_name = "prms_mac_m1_ifort_dbl_prec"
        else:
            exe_name = "prms_mac_intel_gfort_dbl_prec"
    elif platform == "linux":
        exe_name = "prms_linux_gfort_dbl_prec"
    exe_pth = pl.Path(f"../../bin/{exe_name}").resolve()
    return exe_pth


def scheduler_active():
    slurm = os.getenv("SLURM_JOB_ID") is not None
    pbs = os.getenv("PBS_JOBID") is not None
    # add more scheduler checks
    return slurm or pbs


def enforce_scheduler(test_dir):
    """Enforce the use of scheduler

    Args:
        test_dir: the domain directory to run or schedule

    Return:
        True if must use a scheduler, False if not
    """

    if scheduler_active():
        return False
    glob_match = list(
        fnmatch(str(test_dir), gg) for gg in domain_globs_schedule
    )
    if any(glob_match):
        msg = (
            f"Skipping domain '{test_dir}' which must be scheduled or use "
            "--force to override skip"
        )
        warn(msg, UserWarning)
        return True

    return False


def collect_simulations(
    domain_list: list, force: bool = True, verbose: bool = False
):
    simulations = {}
    for test_dir in test_dirs:
        # ensure this is a self-contained run (all files in repo)
        if not (test_dir / "prcp.cbh").exists():
            continue

        # filter selected domains
        if len(domain_list) and (test_dir.name not in domain_list):
            continue

        # optionally enforce scheduler
        if not force:
            skip = enforce_scheduler(test_dir)
            if skip:
                continue

        # if control file is found, add simulation
        control_files = test_dir.glob("*.control")
        for control in control_files:
            id = f"{test_dir.name}_{control.with_suffix('').name}"
            simulations[id] = {
                "ws": test_dir,
                "control_file": control,
            }

    # make sure all requested domains were found
    if len(domain_list) and (len(simulations) < len(domain_list)):
        requested = set(domain_list)
        found = [pl.Path(dd).name for dd in simulations.keys()]
        requested_not_found = requested.difference(found)
        msg = (
            f"The following requested domains were not found: "
            f"{requested_not_found}"
        )
        pytest.exit(msg)

    if verbose:
        print("\nrun_domains.py found the following domains to run:\n")
        print(f"{list(simulations.keys())}")

    return simulations


def collect_csv_files(simulations: list) -> List[pl.Path]:
    csv_files = []
    for key, value in simulations.items():
        output_pth = pl.Path(key) / "output"
        csv_files_dom = sorted(output_pth.glob("*.csv"))
        csv_files += [ff for ff in csv_files_dom if ff.name != "stats.csv"]
    return csv_files


def pytest_generate_tests(metafunc):
    domain_list = metafunc.config.getoption("domain")
    force = metafunc.config.getoption("force")
    simulations = collect_simulations(domain_list, force)
    csv_files = collect_csv_files(simulations)

    if "csv_file" in metafunc.fixturenames:
        ids = [ff.parent.parent.name + ":" + ff.name for ff in csv_files]
        metafunc.parametrize("csv_file", csv_files, ids=ids)

    if "soltab_file" in metafunc.fixturenames:
        soltab_files = [
            pl.Path(kk) / "soltab_debug" for kk in simulations.keys()
        ]
        ids = [ff.parent.name + ":" + ff.name for ff in soltab_files]
        metafunc.parametrize(
            "soltab_file", soltab_files, ids=ids, scope="session"
        )

    if "simulation" in metafunc.fixturenames:
        metafunc.parametrize(
            "simulation",
            list(simulations.values()),
            ids=list(simulations.keys()),
            scope="session",
        )

    if "final_file" in metafunc.fixturenames:
        final_files = [
            pl.Path(kk) / var
            for kk in simulations.keys()
            for var in final_var_names
        ]
        ids = [ff.parent.name + ":" + ff.name for ff in final_files]
        # these are not really file names they are domain/key_var
        metafunc.parametrize(
            "final_file", final_files, ids=ids, scope="session"
        )
