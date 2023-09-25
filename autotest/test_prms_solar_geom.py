import numpy as np
import pytest

from pywatershed.atmosphere.prms_solar_geometry import PRMSSolarGeometry
from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.base.parameters import Parameters
from pywatershed.parameters import PrmsParameters

from utils_compare import compare_in_memory, compare_netcdfs

# in this case we'll compare netcdf files and in memory
atol = rtol = np.finfo(np.float32).resolution

params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


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
    output_dir = domain["prms_output_dir"]
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
    solar_geom.output()
    solar_geom.finalize()

    compare_netcdfs(
        PRMSSolarGeometry.get_variables(),
        tmp_path,
        output_dir,
        atol=atol,
        rtol=rtol,
    )

    answers = {}
    for var in PRMSSolarGeometry.get_variables():
        var_pth = output_dir / f"{var}.nc"
        answers[var] = adapter_factory(
            var_pth, variable_name=var, control=control
        )

    # check the advance/calculate the state
    sunhrs_id = id(solar_geom.soltab_sunhrs)

    for ii in range(control.n_times):
        control.advance()
        solar_geom.advance()
        solar_geom.calculate(1.0)

        compare_in_memory(solar_geom, answers, atol=atol, rtol=rtol)
        assert id(solar_geom.soltab_sunhrs) == sunhrs_id

    return
