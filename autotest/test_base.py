from datetime import datetime, timedelta

import numpy as np
import pytest

from pynhm.base.accessor import Accessor
from pynhm.base.control import Control
from pynhm.base.storageUnit import StorageUnit
from pynhm.base.Time import Time
from pynhm.utils.parameters import PrmsParameters


class TestAccessor:
    def test_accessor(self):
        da = Accessor()
        assert isinstance(da, Accessor)
        key, value = ("a", 0)
        da[key] = value
        assert da[key] is value
        del da[key]
        with pytest.raises(AttributeError):
            da[key]
        return


time_dict_g = {
    "start_time": np.datetime64("1979-01-03T00:00:00.00"),
    "end_time": np.datetime64("1979-01-04T00:00:00.00"),
    "time_step": np.timedelta64(1, "D"),
}


@pytest.fixture(scope="function")
def time_dict():
    return time_dict_g


class TestControl:
    def test_control(self, time_dict):
        control = Control(**time_dict)
        assert control.start_time == time_dict_g["start_time"]
        assert control.end_time == time_dict_g["end_time"]
        assert control.time_step == time_dict_g["time_step"]
        assert control.current_time == time_dict_g["start_time"]
        assert control.previous_time is None
        assert control.n_times == 2
        assert control.i_time == 0
        control.advance()
        assert control.start_time == time_dict_g["start_time"]
        assert control.end_time == time_dict_g["end_time"]
        assert control.time_step == time_dict_g["time_step"]
        assert control.current_time == time_dict_g["end_time"]
        assert control.previous_time == time_dict_g["start_time"]
        assert control.n_times == 2
        assert control.i_time == 1
        with pytest.raises(ValueError):
            control.advance()
        return None

    def test_init_load(self, domain):
        control = Control.load(domain["control_file"])
        return None
