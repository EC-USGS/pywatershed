from datetime import datetime

import numpy as np
import pytest

from pywatershed.base.control import Control
from pywatershed.hydrology.prms_canopy import PRMSCanopy
from pywatershed.parameters import PrmsParameters  # # TODO: too specific

time_dict = {
    "start_time": np.datetime64("1979-01-03T00:00:00.00"),
    "end_time": np.datetime64("1979-01-06T00:00:00.00"),
    "time_step": np.timedelta64(1, "D"),
}

nhru = 2


@pytest.fixture(scope="function")
def params_simple():
    prms_params = {
        "dims": {"nhru": nhru},
        "data_vars": {
            "hru_area": np.array(nhru * [1.0]),
            "covden_sum": np.array(nhru * [0.5]),
            "covden_win": np.array(nhru * [0.5]),
            "srain_intcp": np.array(nhru * [1.0]),
            "wrain_intcp": np.array(nhru * [1.0]),
            "snow_intcp": np.array(nhru * [1.0]),
            "epan_coef": np.array(nhru * [1.0]),
            "potet_sublim": np.array(nhru * [1.0]),
            "cov_type": np.array(nhru * [1]),
        },
    }
    prms_params["metadata"] = {}
    for kk in prms_params["data_vars"].keys():
        prms_params["metadata"][kk] = {"dims": ("nhru",)}

    return PrmsParameters(**prms_params)


@pytest.fixture(scope="function")
def control_simple():
    return Control(**time_dict)


def test_control_simple(control_simple):
    assert control_simple.options == {}
    ts = time_dict["time_step"]
    assert control_simple.time_step == ts
    assert control_simple.start_time == time_dict["start_time"]
    assert control_simple.end_time == time_dict["end_time"]
    assert control_simple.current_time is None
    assert control_simple.itime_step == -1
    prev_time = control_simple.current_time
    n_times = control_simple.n_times
    assert n_times == 4

    for ii in range(n_times):
        control_simple.advance()
        assert prev_time == control_simple.previous_time
        assert control_simple.itime_step == ii
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
        assert dowy == (control_simple.current_dowy - 1)

        prev_time = control_simple.current_time

    with pytest.raises(ValueError):
        control_simple.advance()


def test_control_advance(control_simple, params_simple):
    # common inputs for 2 canopies
    input_variables = {}
    for key in PRMSCanopy.get_inputs():
        input_variables[key] = np.ones([nhru])

    # todo: this is testing instantiation, but not physics
    # ntimes = control.n_times
    cnp1 = PRMSCanopy(
        control=control_simple,
        discretization=None,
        parameters=params_simple,
        **input_variables,
        verbose=True,
    )
    cnp1.name = "cnp1"

    cnp2 = PRMSCanopy(
        control=control_simple,
        discretization=None,
        parameters=params_simple,
        **input_variables,
        verbose=True,
    )
    cnp2.name = "cnp2"

    # Advance correctly
    control_simple.advance()
    for cnp in [cnp1, cnp2]:
        cnp.advance()
        assert control_simple.itime_step == cnp._itime_step
        # for ii in cnp.inputs:
        #     assert (
        #         cnp._input_variables_dict[ii]._itime_step
        #         == control_simple.itime_step
        #     )

    # This is unnecessary?
    cnp1.calculate(time_length=1.0)
    cnp2.calculate(time_length=1.0)

    # dont advance controller
    for cnp in [cnp1, cnp2]:
        cnp.advance()
        assert control_simple.itime_step == cnp._itime_step
        # for ii in cnp.inputs:
        #     assert (
        #         cnp._input_variables_dict[ii]._itime_step
        #         == control_simple.itime_step
        #     )

    # Advance correctly
    control_simple.advance()
    for cnp in [cnp1, cnp2]:
        cnp.advance()
        assert control_simple.itime_step == cnp._itime_step
        # for ii in cnp.inputs:
        #     assert (
        #         cnp._input_variables_dict[ii]._itime_step
        #         == control_simple.itime_step
        #     )


def test_init_load(domain):
    control = Control.load(domain["control_file"])
    return None
