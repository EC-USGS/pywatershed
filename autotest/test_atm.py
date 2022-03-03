import pathlib as pl
from copy import deepcopy
from datetime import datetime, timedelta

import numpy as np
import pytest

from pynhm.atmosphere.AtmBoundaryLayer import AtmBoundaryLayer

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


# Parameterize options to init?
@pytest.fixture
def atm_init():
    atm = AtmBoundaryLayer(**atm_init_test_dict)
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
    def test_advance(self, time_step):
        init_dict = deepcopy(atm_init_test_dict)
        init_dict["time_step"] = time_step
        atm = AtmBoundaryLayer(**init_dict)
        atm.set_state("datetime", test_time)
        try:
            atm.advance()
            assert atm.current_time == (
                init_dict["start_time"] + atm.time_step
            )
            if time_step == test_time_steps[1]:
                assert False
        except ValueError:
            if time_step == test_time_steps[1]:
                assert True
            else:
                assert False
        return
