import pathlib as pl

import numpy as np
import pytest
from utils import assert_or_print

from pywatershed.utils import ControlVariables, compare_control_files

test_answers = {
    "hru_1": {
        "start_time": np.datetime64("1979-01-01T00:00:00"),
        "end_time": np.datetime64("2019-12-31T00:00:00"),
        "initial_deltat": np.timedelta64(24, "h"),
    },
    "drb_2yr": {
        "start_time": np.datetime64("1979-01-01T00:00:00"),
        "end_time": np.datetime64("1980-12-31T00:00:00"),
        "initial_deltat": np.timedelta64(24, "h"),
    },
    "ucb_2yr": {
        "start_time": np.datetime64("1979-01-01T00:00:00"),
        "end_time": np.datetime64("1980-12-31T00:00:00"),
        "initial_deltat": np.timedelta64(24, "h"),
    },
}


control_keys = tuple(
    (
        "start_time",
        "end_time",
        "initial_deltat",
    )
)


@pytest.mark.domain
def test_control_read(simulation):
    domain_name = simulation["name"].split(":")[0]
    if domain_name not in test_answers.keys():
        pytest.skip(f"Answers not provided for domain {domain_name}")
    domain_answers = test_answers[domain_name]

    control_file = simulation["control_file"]
    print(f"parsing...'{control_file}'")
    control = ControlVariables.load(control_file)

    results = {
        key: val
        for key, val in control.control.items()
        if key in domain_answers.keys()
    }
    assert_or_print(results, domain_answers, print_ans=simulation["print_ans"])
    print(f"success parsing...'{control_file}'")
    return


@pytest.mark.domainless
def test_control_compare():
    common_dir = pl.Path("../test_data/common")
    n_diffs = compare_control_files(
        common_dir / "control.single_hru",
        common_dir / "control.multi_hru",
        silent=True,
    )
    assert n_diffs == 5
    return
