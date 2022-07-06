from datetime import datetime, timedelta

import numpy as np
import pytest

from pynhm.atmosphere.PRMSSolarGeometry import PRMSSolarGeometry
from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.utils.parameters import PrmsParameters

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
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


@pytest.mark.parametrize(
    "from_prms_file", (True, False), ids=("from_prms_file", "compute")
)
def test_solar_geom(domain, control, from_prms_file, tmp_path):

    prms_soltab_file = domain["prms_run_dir"] / "soltab_debug"
    if from_prms_file:
        from_prms_file = prms_soltab_file
    else:
        from_prms_file = None

    solar_geom = PRMSSolarGeometry(
        control, from_prms_file=from_prms_file, netcdf_output_dir=tmp_path
    )

    ans = {}
    for vv in solar_geom.variables:
        nc_path = domain["prms_output_dir"] / f"{vv}.nc"
        ans[vv] = adapter_factory(nc_path, variable_name=vv, control=control)

    # check the shapes
    for vv in solar_geom.variables:
        assert ans[vv]._dataset.dataset[vv].shape == solar_geom[f"_{vv}"].shape

    # check the 2D values
    atol = 1e-5
    for vv in solar_geom.variables:
        assert np.allclose(
            solar_geom[f"_{vv}"], ans[vv]._dataset.dataset[vv], atol=atol
        )

    # check the advance/calculate the state
    sunhrs_id = id(solar_geom.soltab_sunhrs)

    for ii in range(4):
        control.advance()
        solar_geom.advance()
        solar_geom.calculate(1.0)

        for vv in solar_geom.variables:
            ans[vv].advance()

            assert solar_geom[vv].shape == ans[vv].current.shape
            assert np.allclose(solar_geom[vv], ans[vv].current)

        assert id(solar_geom.soltab_sunhrs) == sunhrs_id

    return
