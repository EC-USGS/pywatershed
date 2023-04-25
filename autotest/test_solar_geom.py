import numpy as np
import pytest

from pynhm.atmosphere.PRMSSolarGeometry import PRMSSolarGeometry
from pynhm.base.control import Control
from pynhm.parameters import PrmsParameters


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


# @pytest.mark.xfail
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

    ans = PRMSSolarGeometry(control, from_prms_file=prms_soltab_file)

    # check the shapes
    for vv in solar_geom.variables:
        assert ans[vv].data.shape == solar_geom[vv].data.shape

    # check the 2D values
    atol = np.finfo(np.float32).resolution
    rtol = atol

    for vv in solar_geom.variables:
        assert np.allclose(
            solar_geom[vv].data,
            ans[vv].data,
            atol=atol,
            rtol=rtol,
        )

    # check the advance/calculate the state
    sunhrs_id = id(solar_geom.soltab_sunhrs)

    for ii in range(4):
        control.advance()
        solar_geom.advance()
        solar_geom.calculate(1.0)

        for vv in solar_geom.variables:
            ans[vv].advance()

            assert solar_geom[vv].current.shape == ans[vv].current.shape
            assert np.allclose(solar_geom[vv].current, ans[vv].current)
            assert np.allclose(
                solar_geom[vv].current, ans[vv].current, atol=atol, rtol=rtol
            )

        assert id(solar_geom.soltab_sunhrs) == sunhrs_id

    return
