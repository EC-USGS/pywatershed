import os
import pathlib as pl
import sys
from fnmatch import fnmatch
from platform import processor
from typing import List
import numpy as np

import pytest

from pywatershed import CsvFile, Soltab


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
    exe_name = "prms"
    platform = sys.platform.lower()
    if platform == "win32":
        exe_name += "_win.exe"
    elif platform == "darwin":
        if processor() == "arm":
            exe_name += "_mac_m1_intel"
        else:
            exe_name += "_mac"
    elif platform == "linux":
        exe_name += "_linux"
    exe_pth = pl.Path(f"../../bin/{exe_name}").resolve()
    return exe_pth


rootdir = ".."
test_dirs = sorted(
    [path for path in pl.Path(rootdir).iterdir() if path.is_dir()]
)


# This would change to handle other/additional schedulers
domain_globs_schedule = ["*conus*"]


def scheduler_active():
    slurm = os.getenv("SLURM_JOB_ID") is not None
    pbs = os.getenv("PBS_JOBID") is not None
    # add more scheduler checks
    return slurm or pbs


def enforce_scheduler(test_dir):
    if scheduler_active():
        return None
    glob_match = list(
        fnmatch(str(test_dir), gg) for gg in domain_globs_schedule
    )
    if any(glob_match):
        raise RuntimeError(
            f"Domain '{test_dir}' must be scheduled (use --force to override)"
        )
    return None


def collect_simulations(domain_list: list, force: bool = True, verbose: bool = False):
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
            enforce_scheduler(test_dir)

        # if control file is found, add simulation
        ctrl_file = next(iter([p for p in test_dir.iterdir() if p.is_file() and p.name == "control.test"]), None)
        if ctrl_file:
            simulations[str(test_dir)] = ctrl_file.name

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


def collect_nc_files(simulations: list, var_list: list):
    sim_dirs = list(simulations.keys())
    nc_files = []
    for var in var_list:
        for sim in sim_dirs:
            nc_files += [(pl.Path(sim) / f"output/{var}.nc")]
    return nc_files


def pytest_generate_tests(metafunc):
    domain_list = metafunc.config.getoption("domain")
    force = metafunc.config.getoption("force")
    simulations = collect_simulations(domain_list, force)
    csv_files = collect_csv_files(simulations)

    key = "simulation"
    if key in metafunc.fixturenames:
        sims = [
            {"ws": key, "control_file": val} for key, val in simulations.items()
        ]
        ids = [pl.Path(k).name for k in simulations.keys()]
        metafunc.parametrize(key, sims, ids=ids, scope="session")

    key = "soltab_file"
    if key in metafunc.fixturenames:
        soltab_files = [pl.Path(k) / "soltab_debug" for k in simulations.keys()]
        ids = [f.parent.name + ":" + f.name for f in soltab_files]
        metafunc.parametrize(key, soltab_files, ids=ids, scope="session")

    key = "csv_file"
    if key in metafunc.fixturenames:
        ids = [f.parent.name + ":" + f.name for f in csv_files]
        metafunc.parametrize(key, csv_files, ids=ids)

    
