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

# This would change to handle other/additional schedulers
domain_globs_schedule = {"*conus*": "slurm"}


def is_scheduled(scheduler):
    if scheduler == "slurm":
        return os.getenv("SLURM_JOB_ID") is not None
    # add more scheduler checks here


simulations = {}
for test_dir in test_dirs:

    # If the test_dir is required to be scheduled, check if we are
    # running in the scheduler
    doms_sched = list(domain_globs_schedule.keys())
    glob_match = list(fnmatch(str(test_dir), kk) for kk in doms_sched)
    if any(glob_match):
        the_glob = doms_sched[glob_match.index(True)]
        scheduler = domain_globs_schedule[the_glob]
        if not is_scheduled(scheduler):
            raise RuntimeError(
                f"Domain '{test_dir}' must be scheduled using {scheduler}"
            )

    for pth in test_dir.iterdir():
        # checking for prcp.cbh ensure this is a self-contained run (all files in repo)
        if (
            (test_dir / "prcp.cbh").exists()
            and pth.is_file()
            and pth.name == "control.test"
        ):
            # add simulation
            simulations[str(test_dir)] = pth.name

            # create output directory if it doesn't exist
            # JLM: should we force remove existing outputs?
            output_dir = pl.Path(test_dir) / "output"
            output_dir.mkdir(parents=True, exist_ok=True)

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


if __name__ == "__main__":
    test_exe_available()

    for key, value in simulations.items():
        test_prms_run(key, value)
