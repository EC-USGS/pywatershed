import os
import pathlib as pl
import sys

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
exe_pth = pl.Path(f"../../bin/{exe_name}")


rootdir = ".."
test_dirs = [path for path in pl.Path(rootdir).iterdir() if path.is_dir()]

simulations = {}
for test_dir in test_dirs:
    for pth in test_dir.iterdir():
        if pth.is_file() and pth.name == "control.test":
            # add simulation
            simulations[str(test_dir)] = pth.name

            # create output directory if it doesn't exist
            output_dir = pl.Path(test_dir) / "output"
            output_dir.mkdir(parents=True, exist_ok=True)


def test_exe_available():
    assert exe_pth.is_file(), f"'{exe_pth}'...does not exist"

    assert os.access(exe_pth, os.X_OK)


@pytest.mark.parametrize(
    "ws,control_file",
    simulations.items(),
    ids=simulations.keys(),
)
def test_prms_run(ws, control_file):
    print(f"running '{control_file}' in {ws}")
    success, buff = run_model(
        exe_pth,
        control_file,
        model_ws=ws,
        normal_msg="Normal completion of PRMS",
    )
    assert success, f"could not run prms model in '{ws}'"


if __name__ == "__main__":
    test_exe_available()

    for key, value in simulations.items():
        test_prms_run(key, value)
