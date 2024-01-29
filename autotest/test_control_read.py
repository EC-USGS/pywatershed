import pathlib as pl

import numpy as np
import pytest

from pywatershed.utils import ControlVariables, compare_control_files
from utils import assert_or_print

test_answers = {
    "start_time": np.datetime64("1979-01-01T00:00:00"),
    "end_time": np.datetime64("1980-12-31T00:00:00"),
    "initial_deltat": np.timedelta64(24, "h"),
}


@pytest.fixture
def control_keys():
    return tuple(
        (
            "start_time",
            "end_time",
            "initial_deltat",
        )
    )


def test_control_read(simulation, control_keys):
    control_file = simulation["control_file"]
    print(f"parsing...'{control_file}'")
    control = ControlVariables.load(control_file)

    results = {
        key: val
        for key, val in control.control.items()
        if key in test_answers.keys()
    }
    assert_or_print(results, test_answers, print_ans=simulation["print_ans"])
    print(f"success parsing...'{control_file}'")
    return


def test_control_compare():
    common_dir = pl.Path("../test_data/common")
    n_diffs = compare_control_files(
        common_dir / "control.single_hru",
        common_dir / "control.multi_hru",
        silent=True,
    )
    assert n_diffs == 5
    return
