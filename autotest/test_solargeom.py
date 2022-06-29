from copy import deepcopy
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from pynhm.atmosphere.PRMSSolarGeometry import PRMSSolarGeometry
from pynhm.base.control import Control
from pynhm.utils.parameters import PrmsParameters
from pynhm.utils.prms5util import load_soltab_debug

test_time = np.arange(
    datetime(1979, 1, 1), datetime(1979, 1, 7), timedelta(days=1)
).astype(np.datetime64)

atm_init_test_dict = {
    "start_time": np.datetime64("1979-01-03T00:00:00.00"),
    "time_step": np.timedelta64(1, "D"),
    "verbosity": 3,
    "height_m": 5,
}

test_time_steps = [atm_init_test_dict["time_step"], np.timedelta64(1, "h")]


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


class TestPRMSSolarGeometry:
    def test_init(self, domain, control, params):
        solar_geom = PRMSSolarGeometry(control, params)

        # answers
        (
            soltab_potsw_ans,
            soltab_horad_potsw_ans,
            soltab_sunhrs_ans,
        ) = load_soltab_debug(domain["prms_run_dir"] / "soltab_debug")

        # check the shapes
        assert soltab_potsw_ans.shape == solar_geom._soltab_potsw.shape
        assert soltab_horad_potsw_ans.shape == solar_geom._soltab_potsw.shape
        assert soltab_sunhrs_ans.shape == solar_geom._soltab_sunhrs.shape

        # check the values
        atol = 1e-4
        assert np.isclose(
            solar_geom._soltab_sunhrs, soltab_sunhrs_ans, atol=atol
        ).all()
        assert np.isclose(
            solar_geom._soltab_potsw, soltab_potsw_ans, atol=atol
        ).all()
        assert np.isclose(
            solar_geom._soltab_horad_potsw,
            soltab_horad_potsw_ans,
            atol=atol,
        ).all()

        data_size = np.product(soltab_potsw_ans.shape)

        # advance/calculate the state
        sunhrs_id = id(solar_geom.soltab_sunhrs)

        for ii in range(4):
            control.advance()
            solar_geom.advance()
            solar_geom.calculate(1.0)
            assert np.allclose(
                solar_geom.soltab_potsw, soltab_potsw_ans[ii, :]
            )
            assert np.allclose(
                solar_geom.soltab_horad_potsw, soltab_horad_potsw_ans[ii, :]
            )
            assert np.allclose(
                solar_geom.soltab_sunhrs, soltab_sunhrs_ans[ii, :]
            )
            assert id(solar_geom.soltab_sunhrs) == sunhrs_id

        return
