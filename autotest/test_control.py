from datetime import datetime

import numpy as np
import pytest

from pynhm.base.control import Control

time_dict = {
    "start_time": np.datetime64("1979-01-03T00:00:00.00"),
    "end_time": np.datetime64("1979-01-06T00:00:00.00"),
    "time_step": np.timedelta64(1, "D"),
}


@pytest.fixture(scope="function")
def control_simple():
    return Control(**time_dict)


def test_control_simple(control_simple):
    assert control_simple.config is None
    ts = time_dict["time_step"]
    assert control_simple.time_step == ts
    assert control_simple.start_time == time_dict["start_time"]
    assert control_simple.end_time == time_dict["end_time"]
    assert control_simple.current_time is None
    prev_time = control_simple.current_time
    n_times = control_simple.n_times
    assert n_times == 4

    for ii in range(n_times):
        control_simple.advance()
        assert prev_time == control_simple.previous_time
        assert control_simple.i_time == ii
        assert control_simple.current_time == time_dict["start_time"] + ii * ts

        current_time = control_simple.current_time

        # This constitutes a test of utils/time_utils.py
        fmt_var = {
            "%Y": "current_year",
            "%m": "current_month",
            "%j": "current_doy",
        }
        for fmt, var in fmt_var.items():
            check = int(datetime.strftime(current_time.astype(datetime), fmt))
            assert check == control_simple[var]

        # test dowy
        year = control_simple.current_year
        month = control_simple.current_month
        year = year if month >= 10 else year - 1
        wy_start = np.datetime64(f"{year}-10-01")
        dowy = (current_time - wy_start).astype("timedelta64[D]")
        assert dowy == control_simple.current_dowy

        prev_time = control_simple.current_time

    with pytest.raises(ValueError):
        control_simple.advance()
