import pathlib as pl

import numpy as np
import pytest

from pynhm.atmosphere.AtmBoundaryLayer import AtmBoundaryLayer


class TestAtmBoundaryLayer:
    def test_init(self):
        test_dict = {
            "start_time": np.datetime64("2000-01-01T12:00:00.00"),
            "time_step": np.timedelta64(1, "h"),
            "verbose": 3,
            "height_m": 5,
        }
        atm = AtmBoundaryLayer(**test_dict)
        for key, val in test_dict.items():
            key_check = key
            if key not in list(atm.__dict__.keys()):
                key_check = f"_{key}"
            assert getattr(atm, key_check) == val
