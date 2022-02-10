import os
import sys

import pytest

cwd = os.getcwd()
if cwd.endswith("autotest"):
    sys.path.append("..")
    rel_path = ".."
elif cwd.endswith("prmsNHMpy"):
    sys.path.append(".")
    rel_path = "."

from pynhm import CbhInput

base_path = os.path.join(rel_path, "test_models")
cbh_files = (
    os.path.join(rel_path, "prms_models", "box01", "input", "tmax.cbh"),
    os.path.join(rel_path, "prms_models", "box01", "input", "tmin.cbh"),
    os.path.join(rel_path, "prms_models", "box01", "input", "prcp.cbh"),
    os.path.join(base_path, "merced", "input", "mercd.data"),
    os.path.join(base_path, "acf", "input", "ACF_current_Ptet.data"),
)

number_columns = (1, 1, 1, 79, 258)
data_types = (
    ("tmax",),
    ("tmin",),
    ("prcp",),
    ("tmax", "tmin", "precip", "runoff"),
    ("ptet",),
)


@pytest.mark.parametrize("idx, cbh_file", enumerate(cbh_files))
def test_cbh_read(idx, cbh_file):
    print(f"parsing...'{cbh_file}'")

    cbh = CbhInput(cbh_file)
    columns = cbh.get_dataframe_columns

    assert len(columns) == number_columns[idx], (
        f"test {idx + 1}: "
        + f"number of columns ({len(columns)}) "
        + f"does not equal {number_columns[idx]}"
    )

    for data_type in data_types[idx]:
        check = len([column for column in columns if data_type in column]) > 0

        assert check, f"test {idx + 1}: {data_type} not in cbh file"

    print(f"success parsing...'{cbh_file}'")

    return


if __name__ == "__main__":

    for idx, cbh_file in enumerate(cbh_files):
        test_cbh_read(idx, cbh_file)
