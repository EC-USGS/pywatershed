import os
import pathlib as pl
import sys

import pytest

from pynhm import CsvFile


def test_csv_to_netcdf(csv_files):
    nc_path = csv_files.with_suffix(".nc")
    CsvFile(csv_files).to_netcdf(nc_path)
    assert nc_path.exists()
