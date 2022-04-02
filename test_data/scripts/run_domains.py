# This requires flopy to be installed
# which is not in the requirements.

import os
import pathlib as pl
import sys
from fnmatch import fnmatch

import pytest
from flopy import run_model

exe_name = "prms"
platform = sys.platform.lower()
if platform == "win32":
    exe_name += "_win.exe"
elif platform == "darwin":
    exe_name += "_mac"
elif platform == "linux":
    exe_name += "_linux"
exe_pth = pl.Path(f"../../bin/{exe_name}").resolve()


rootdir = ".."
test_dirs = sorted(
    [path for path in pl.Path(rootdir).iterdir() if path.is_dir()]
)

simulations = {}


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
            f"Domain '{test_dir}' must be scheduled (use -f to force)"
        )
    return None


def collect_simulations(force):
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

                if not force:
                    enforce_scheduler(test_dir)

                # add simulation
                simulations[str(test_dir)] = pth.name

                # delete the existing output dir and re-create it
                output_dir = pl.Path(test_dir) / "output"
                if output_dir.exists():
                    shutil.rmtree(output_dir)
                output_dir.mkdir(parents=True)

    print("\nrun_domains.py found the following domains to run:\n")
    print(f"{list(simulations.keys())}")


def test_exe_available():
    assert exe_pth.is_file(), f"'{exe_pth}'...does not exist"
    assert os.access(exe_pth, os.X_OK)


@pytest.mark.parametrize(
    "ws,control_file",
    simulations.items(),
    ids=simulations.keys(),
)
def test_prms_run(ws, control_file):
    print(f"\n\n\n{'*' * 70}\n{'*' * 70}")
    print(f"run_domains.py: Running '{control_file}' in {ws}\n\n", flush=True)
    success, buff = run_model(
        exe_pth,
        control_file,
        model_ws=ws,
        cargs=[
            "-MAXDATALNLEN",
            "60000",
        ],
        normal_msg="Normal completion of PRMS",
    )
    assert success, f"could not run prms model in '{ws}'"

    print(f"run_domains.py: End of domain {ws}\n", flush=True)


def parse_arguments():

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--force",
        required=False,
        action="store_true",
        default=False,
        help=("Force run without scheduler present"),
    )

    args = parser.parse_args()
    return args.force


if __name__ == "__main__":

    force = parse_arguments()
    collect_simulations(force)

    test_exe_available()

    for key, value in simulations.items():
        test_prms_run(key, value)
