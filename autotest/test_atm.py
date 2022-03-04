import pathlib as pl
from copy import deepcopy
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from pynhm.atmosphere.AtmBoundaryLayer import AtmBoundaryLayer
from pynhm.atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from utils import assert_or_print

test_time = np.arange(
    datetime(2000, 1, 1), datetime(2000, 7, 1), timedelta(days=1)
).astype(np.datetime64)

atm_init_test_dict = {
    "start_time": np.datetime64("2000-01-10T00:00:00.00"),
    "time_step": np.timedelta64(1, "D"),
    "verbose": 3,
    "height_m": 5,
}

test_time_steps = [atm_init_test_dict["time_step"], np.timedelta64(1, "h")]


@pytest.fixture
def atm_init():
    atm = AtmBoundaryLayer(**atm_init_test_dict)
    return atm


@pytest.fixture
def atm_nhm_init(domain):
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

    def test_state_roundtrip(self, atm_init):
        # use datetime as state
        atm_init.set_state("datetime", test_time)
        assert (atm_init.get_state("datetime") == test_time).all()
        return

    def test_current_state_roundtrip(self, atm_init):
        # use datetime as state
        atm_init.set_state("datetime", test_time)
        assert (
            atm_init.get_current_state("datetime")[0]
            == atm_init_test_dict["start_time"]
        )
        return

    def test_current_time_index(self, atm_init):
        atm_init.set_state("datetime", test_time)
        result = atm_init.current_time_index
        assert result == np.array([9])
        return

    def test_current_time(self, atm_init):
        atm_init.set_state("datetime", test_time)
        assert atm_init.current_time == atm_init_test_dict["start_time"]
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
        atm_init.set_state("datetime", test_time)
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
            var: atm_nhm_init.get_state(var).mean()
            for var in atm_nhm_init.state_vars_in_file
        }
        ds = xr.open_dataset(domain["cbh_nc"])
        # preference for adjusted over raw fields in the code
        file_vars = [
            f"{var}_adj" if f"{var}_adj" in ds.variables else var
            for var in atm_nhm_init.state_vars_in_file
        ]
        answers = {var: ds[var].mean().values.tolist() for var in file_vars}
        for key_res in results.keys():
            key_ans = key_res
            if key_ans not in answers:
                key_ans = f"{key_res}_adj"
            assert np.isclose(results[key_res], answers[key_ans])
        return

    def test_state_roundtrip(self, atm_nhm_init):
        # use datetime as state
        atm_nhm_init.set_state("datetime", test_time)
        assert (atm_nhm_init.get_state("datetime") == test_time).all()
        return

    # def test_current_state_roundtrip(self, atm_init):
    #     # use datetime as state
    #     atm_init.set_state("datetime", test_time)
    #     assert (
    #         atm_init.get_current_state("datetime")[0]
    #         == atm_init_test_dict["start_time"]
    #     )
    #     return

    # def test_current_time_index(self, atm_init):
    #     atm_init.set_state("datetime", test_time)
    #     result = atm_init.current_time_index
    #     assert result == np.array([9])
    #     return

    # def test_current_time(self, atm_init):
    #     atm_init.set_state("datetime", test_time)
    #     assert atm_init.current_time == atm_init_test_dict["start_time"]
    #     return

    # def test_time_step(self, atm_init):
    #     assert atm_init.time_step == atm_init_test_dict["time_step"]
    #     return

    # @pytest.mark.parametrize(
    #     "time_step", test_time_steps, ids=["valid", "invalid"]
    # )
    # def test_advance(self, atm_init, time_step):
    #     # the easy way is to hack private data
    #     atm_init._time_step = time_step
    #     atm_init.set_state("datetime", test_time)
    #     try:
    #         atm_init.advance()
    #         assert atm_init.current_time == (
    #             atm_init_test_dict["start_time"] + time_step
    #         )
    #         if time_step == test_time_steps[1]:
    #             assert False
    #     except ValueError:
    #         if time_step == test_time_steps[1]:
    #             assert True
    #         else:
    #             assert False
    #     return
