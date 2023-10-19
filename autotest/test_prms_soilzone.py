import pathlib as pl

import pytest

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.hydrology.prms_soilzone import PRMSSoilzone
from pywatershed.parameters import Parameters, PrmsParameters
from utils_compare import compare_in_memory, compare_netcdfs

# compare in memory (faster) or full output files? or both!
do_compare_output_files = False
do_compare_in_memory = True
rtol = atol = 1.0e-6

calc_methods = ("numpy", "numba")
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
        param_file = domain["dir"] / "parameters_PRMSSoilzone.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(
    domain, control, discretization, parameters, tmp_path, calc_method
):
    tmp_path = pl.Path(tmp_path)

    comparison_var_names = list(
        set(PRMSSoilzone.get_variables())
        # these are not prms variables per se
        # the hru ones have non-hru equivalents being checked
        # soil_zone_max and soil_lower_max would be nice to check but
        # prms5.2.1 wont write them as hru variables
        - set(
            [
                "perv_actet_hru",
                "soil_lower_change_hru",
                "soil_lower_max",
                "soil_rechr_change_hru",
                "soil_zone_max",  # not a prms variable?
            ]
        )
    )

    output_dir = domain["prms_output_dir"]

    input_variables = {}
    for key in PRMSSoilzone.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = nc_path

    soil = PRMSSoilzone(
        control=control,
        discretization=discretization,
        parameters=parameters,
        **input_variables,
        budget_type="error",
        calc_method=calc_method,
    )

    if do_compare_output_files:
        nc_parent = tmp_path / domain["domain_name"]
        soil.initialize_netcdf(nc_parent)

    if do_compare_in_memory:
        answers = {}
        for var in comparison_var_names:
            var_pth = output_dir / f"{var}.nc"
            answers[var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    for istep in range(control.n_times):
        control.advance()
        soil.advance()
        soil.calculate(1.0)
        if do_compare_in_memory:
            compare_in_memory(
                soil, answers, atol=atol, rtol=rtol, skip_missing_ans=True
            )

    soil.finalize()

    if do_compare_output_files:
        compare_netcdfs(
            comparison_var_names,
            tmp_path / domain["domain_name"],
            output_dir,
            atol=atol,
            rtol=rtol,
        )

    return
