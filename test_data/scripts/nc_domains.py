import os
import pathlib as pl
import sys

import pytest

from pynhm import CsvFile

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

csv_files = []
for key, value in simulations.items():
    output_pth = pl.Path(key) / "output"
    csv_files += [
        str(file_pth)
        for file_pth in output_pth.iterdir()
        if file_pth.is_file()
        and file_pth.suffix.lower() == ".csv"
        and file_pth.name != "stats.csv"
    ]


@pytest.mark.parametrize(
    "file_path",
    csv_files,
    ids=csv_files,
)
def test_csv_to_netcdf(file_path):
    nc_pth = str(pl.Path(file_path)).replace(".csv", ".nc")
    CsvFile(file_path).to_netcdf(nc_pth)
