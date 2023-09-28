import numpy as np
import pytest

from pywatershed.atmosphere.prms_solar_geometry import PRMSSolarGeometry
from pywatershed.base.control import Control
from pywatershed.base.parameters import Parameters
from pywatershed.parameters import PrmsParameters

params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(domain):
    return Control.load_prms(domain["control_file"], warn_unused_options=False)


@pytest.fixture(scope="function")
def discretization(domain):
    dis_hru_file = domain["dir"] / "parameters_dis_hru.nc"
    return Parameters.from_netcdf(dis_hru_file, encoding=False)


@pytest.fixture(scope="function", params=params)
def parameters(domain, request):
    if request.param == "params_one":
        params = PrmsParameters.load(domain["param_file"])
    else:
        param_file = domain["dir"] / "parameters_PRMSSolarGeometry.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.parametrize(
    "from_prms_file", (True, False), ids=("from_prms_file", "compute")
)
def test_compare_prms(
    domain, control, discretization, parameters, tmp_path, from_prms_file
):
    prms_soltab_file = domain["prms_run_dir"] / "soltab_debug"
    if from_prms_file:
        from_prms_file = prms_soltab_file
    else:
        from_prms_file = None

    solar_geom = PRMSSolarGeometry(
        control,
        discretization=discretization,
        parameters=parameters,
        from_prms_file=from_prms_file,
        netcdf_output_dir=tmp_path,
    )
    solar_geom.finalize()

    ans = PRMSSolarGeometry(
        control,
        discretization=discretization,
        parameters=parameters,
        from_prms_file=prms_soltab_file,
    )

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
