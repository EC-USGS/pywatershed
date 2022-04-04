from fnmatch import fnmatch

import os
import pathlib as pl
import shutil
import sys

import pytest


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


def collect_simulations(domain_list: list, force: bool):
    simulations = {}
    for test_dir in test_dirs:

        # If the test_dir is required to be scheduled, check if we are
        # running in the scheduler

        for pth in test_dir.iterdir():
            # checking for prcp.cbh ensure this is a self-contained run (all files in repo)
            if (
                (test_dir / "prcp.cbh").exists()
                and pth.is_file()
                and pth.name == "control.test"
            ):

                if len(domain_list) and (test_dir.name not in domain_list):
                    continue

                if not force:
                    enforce_scheduler(test_dir)

                # add simulation
                simulations[str(test_dir)] = pth.name

                # delete the existing output dir and re-create it
                output_dir = pl.Path(test_dir) / "output"
                if output_dir.exists():
                    output_rm = pl.Path(test_dir) / "rm_output"
                    output_dir.rename(output_rm)
                    shutil.rmtree(output_rm)
                output_dir.mkdir(parents=True)

    if len(domain_list) and (len(simulations) < len(domain_list)):
        requested = set(domain_list)
        found = [pl.Path(dd).name for dd in simulations.keys()]
        requested_not_found = requested.difference(found)
        msg = f"The following requested domains were not found: {requested_not_found}"
        pytest.exit(msg)

    print("\nrun_domains.py found the following domains to run:\n")
    print(f"{list(simulations.keys())}")
    return simulations


def pytest_generate_tests(metafunc):
    if "simulation" in metafunc.fixturenames:
        domain_list = metafunc.config.getoption("domain")
        force = metafunc.config.getoption("force")
        simulations = collect_simulations(domain_list, force)
        sim_list = [
            {"ws": key, "control_file": val}
            for key, val in simulations.items()
        ]
        ids = [pl.Path(ss).name for ss in simulations.keys()]
        metafunc.parametrize("simulation", sim_list, ids=ids)
