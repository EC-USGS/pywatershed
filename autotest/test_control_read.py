import pathlib as pl

import numpy as np
import pytest

from pywatershed.utils import ControlVariables, compare_control_files
from utils import assert_or_print


@pytest.fixture
def control_keys():
    return tuple(
        (
            "start_time",
            "end_time",
            "initial_deltat",
        )
    )


def test_control_read(domain, control_keys):
    control_file = domain["control_file"]
    print(f"parsing...'{control_file}'")
    control = ControlVariables.load(control_file)
    answers = domain["test_ans"]["control_read"]

    for key, value in answers.items():
        if key in (
            "start_time",
            "end_time",
        ):
            answers[key] = np.datetime64(value)
        elif key in ("initial_deltat",):
            answers[key] = np.timedelta64(int(value), "h")

    results = {
        key: val
        for key, val in control.control.items()
        if key in answers.keys()
    }
    assert_or_print(results, answers, print_ans=domain["print_ans"])
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
