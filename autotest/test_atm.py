import pathlib as pl
from copy import deepcopy
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from pynhm.atmosphere.AtmBoundaryLayer import AtmBoundaryLayer
from pynhm.atmosphere.NHMBoundaryLayer import NHMBoundaryLayer

from pynhm.utils.parameters import PrmsParameters
from pynhm.utils.prms5util import load_prms_statscsv

test_time = np.arange(
    datetime(1979, 1, 1), datetime(1979, 7, 1), timedelta(days=1)
).astype(np.datetime64)

test_time_match_start = np.arange(
    datetime(1979, 1, 10), datetime(1979, 7, 1), timedelta(days=1)
).astype(np.datetime64)

atm_init_test_dict = {
    "start_time": np.datetime64("1979-01-10T00:00:00.00"),
    "time_step": np.timedelta64(1, "D"),
    "verbose": 3,
    "height_m": 5,
}

test_time_steps = [atm_init_test_dict["time_step"], np.timedelta64(1, "h")]


@pytest.fixture
def atm_init(scope="function"):
    atm = AtmBoundaryLayer(**atm_init_test_dict)
    return atm


@pytest.fixture
def atm_nhm_init(domain, scope="function"):
    atm = NHMBoundaryLayer(domain["cbh_nc"], **atm_init_test_dict)
    return atm


class TestAtmBoundaryLayer:
    def test_init(self, atm_init):
        for key, val in atm_init_test_dict.items():
            key_check = key
            if key not in list(atm_init.__dict__.keys()):
                key_check = f"_{key}"
            assert getattr(atm_init, key_check) == val
        return

    @pytest.mark.parametrize(
        "the_state", ["datetime", "foo"], ids=["valid", "invalid"]
    )
    def test_state_roundtrip(self, atm_init, the_state):
        # Has to be a round trip because there's no defined time or state
        # get and set state via __setitem__ and __getitem__
        # Could test other errors of set and get
        try:
            atm_init[the_state] = test_time
            if the_state in ["foo"]:
                assert False
        except KeyError:
            if the_state in ["foo"]:
                assert True
            else:
                assert False

        try:
            assert (atm_init[the_state] == test_time).all()
            if the_state in ["foo"]:
                assert False
        except KeyError:
            if the_state in ["foo"]:
                assert True
            else:
                assert False

        return

    def test_time_step(self, atm_init):
        assert atm_init.time_step == atm_init_test_dict["time_step"]
        return

    @pytest.mark.parametrize(
        "time_step", test_time_steps, ids=["valid", "invalid"]
    )
    def test_advance(self, atm_init, time_step):
        # the easy way is to hack private data
        atm_init._time_step = time_step
        atm_init["datetime"] = test_time
        try:
            atm_init.advance()
            assert atm_init.current_time == (
                atm_init_test_dict["start_time"] + time_step
            )
            if time_step == test_time_steps[1]:
                assert False
        except ValueError:
            if time_step == test_time_steps[1]:
                assert True
            else:
                assert False
        return

    def test_current_time_index(self, atm_init):
        # Test start time index by construction that
        # start_time is the 10th day in the test_time
        atm_init["datetime"] = test_time
        assert atm_init.current_time_index == np.array([9])
        # Test current time index after advance
        atm_init["datetime"] = test_time_match_start
        n_adv = 5
        for aa in range(n_adv):
            atm_init.advance()
        assert atm_init.current_time_index == np.array([n_adv])
        return

    def test_current_time(self, atm_init):
        # Test start time index by construction that
        # start_time is the 10th day in the test_time
        atm_init["datetime"] = test_time
        assert atm_init.current_time == atm_init_test_dict["start_time"]
        # Test current time index after advance
        atm_init["datetime"] = test_time_match_start
        n_adv = 5
        for aa in range(n_adv):
            atm_init.advance()
        assert atm_init.current_time == test_time_match_start[n_adv]
        return

    def test_current_state_roundtrip(self, atm_init):
        # Must be a round trip since there's no defined time/state
        atm_init["datetime"] = test_time_match_start
        n_adv = 5
        for aa in range(n_adv):
            atm_init.advance()
        assert (
            atm_init.get_current_state("datetime")
            == test_time_match_start[n_adv]
        )
        return


class TestNHMBoundaryLayer:
    def test_init(self, domain, atm_nhm_init):
        init_test_dict = deepcopy(atm_init_test_dict)
        init_test_dict["nc_file"] = domain["cbh_nc"]
        # Check all the init arguments
        for key, val in init_test_dict.items():
            key_check = key
            if key not in list(atm_nhm_init.__dict__.keys()):
                key_check = f"_{key}"
            assert getattr(atm_nhm_init, key_check) == val
        # Check state - this tests get_state.
        results = {
            var: atm_nhm_init[var].mean() for var in atm_nhm_init.variables
        }
        ds = xr.open_dataset(domain["cbh_nc"])
        # preference for adjusted over raw fields in the code
        file_vars = [
            f"{var}_adj" if f"{var}_adj" in ds.variables else var
            for var in atm_nhm_init.variables
        ]
        answers = {var: ds[var].mean().values.tolist() for var in file_vars}
        for key_res in results.keys():
            key_ans = key_res
            if key_ans not in answers:
                key_ans = f"{key_res}_adj"
            assert np.isclose(results[key_res], answers[key_ans])
        return

    def test_state_roundtrip(self, atm_nhm_init):
        # Roundtrip and a half just to use existing values: get, set, get
        mult = 2
        values = atm_nhm_init["prcp"] * mult
        assert (values / atm_nhm_init.prcp == mult).all()  # while we are here
        atm_nhm_init["prcp"] = values
        # One better, these are actually the same object, could test "is"
        assert (atm_nhm_init["prcp"] == values).all()
        return

    def test_current_state(self, atm_nhm_init):
        n_adv = 5
        start_time_index = atm_nhm_init.current_time_index
        answer = deepcopy(atm_nhm_init["prcp"])[start_time_index + n_adv, :]
        for aa in range(n_adv):
            atm_nhm_init.advance()
        assert (atm_nhm_init.get_current_state("prcp") == answer).all()
        return

    @pytest.mark.parametrize(
        "read_vars",
        [["prcp", "tmax", "tmin"], ["foo"], ["tmin_adj"]],
        ids=["toadj", "invalid", "preadj"],
    )
    def test_nc_read_vars(self, atm_nhm_init, read_vars):
        atm_read_dict = deepcopy(atm_init_test_dict)
        atm_read_dict["nc_read_vars"] = read_vars
        try:
            atm_read = NHMBoundaryLayer(atm_nhm_init._nc_file, **atm_read_dict)
            if read_vars[0] in ["foo"]:
                assert False
            else:
                assert set(atm_read.variables) == set(
                    [
                        var.split("_")[0]
                        for var in atm_read_dict["nc_read_vars"]
                    ]
                )
        except ValueError:
            if read_vars[0] in ["foo"]:
                assert True
            else:
                assert False
        return

    # JLM: def test_nc_read_vars

    def test_state_adj(self, domain, atm_nhm_init):
        atm_adj_dict = deepcopy(atm_init_test_dict)
        atm_adj_dict["nc_read_vars"] = ["prcp", "tmax", "tmin"]
        atm_to_adj = NHMBoundaryLayer(atm_nhm_init._nc_file, **atm_adj_dict)
        params = PrmsParameters(domain["param_file"])
        atm_to_adj.param_adjust(params)

        for vv in atm_to_adj.variables:
            # JLM: we loose rhavg in the adjustments... it's a bit of a mystery
            # if that should also be adjusted. checking on that w parker
            if vv in atm_nhm_init.variables:
                assert np.isclose(
                    atm_to_adj[vv],
                    atm_nhm_init[vv],
                    atol=1e-06,
                ).all()

        return

    def test_solar_radiation(self, domain, atm_nhm_init):
        params = PrmsParameters(domain["param_file"])
        # assert output matches output on file.
        prms_output_file = domain["prms_outputs"]["swrad"]
        prms_output = load_prms_statscsv(prms_output_file)
        swrad_ans_dates = prms_output.index.values
        swrad_ans_array = prms_output.to_numpy()
        wh_dates = np.where(np.isin(atm_nhm_init["datetime"], swrad_ans_dates))

        swrad = atm_nhm_init.calculate_sw_rad_degree_day(params)

        result = np.isclose(
            swrad_ans_array,
            swrad[wh_dates, :],
            rtol=1e-121,
            atol=1e-04,  # Only the atol matters here, if atol < 1e-4 fails
        )
        assert result.all()

        return
