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
cbh_files = (os.path.join(base_path, "merced", "input", "mercd.data"),)

number_columns = (79,)


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

    print(f"success parsing...'{cbh_file}'")

    return


if __name__ == "__main__":

    for idx, cbh_file in enumerate(cbh_files):
        test_cbh_read(idx, cbh_file)
