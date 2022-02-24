import os
import pathlib as pl
import sys

import pandas as pd
import pytest

cwd = os.getcwd()
if cwd.endswith("autotest"):
    sys.path.append("..")
    rel_path = ".."
elif cwd.endswith("pynhm"):
    sys.path.append(".")
    rel_path = "."

from pynhm.preprocess import CsvFile

test_dirs = (
    pl.Path(f"{rel_path}/test_data/box_01/output"),
    pl.Path(f"{rel_path}/test_data/drb_2yr/output"),
)
csv_files = {}
for test_dir in test_dirs:
    if test_dir.is_dir():
        csv_files[test_dir] = list(test_dir.glob("**/*.csv"))
csv_keys = tuple(csv_files.keys())


@pytest.mark.parametrize("csv_key", csv_keys)
def test_single_csv(csv_key):
    csv = CsvFile(name=csv_files[csv_key][0])

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


@pytest.mark.parametrize("csv_key", csv_keys)
def test_multiple_csv(csv_key):
    csv = CsvFile()
    for name in csv_files[csv_key]:
        csv.add_path(name)

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


if __name__ == "__main__":
    for csv_key in csv_keys:
        test_single_csv(csv_key)

    for csv_key in csv_keys:
        test_multiple_csv(csv_key)
