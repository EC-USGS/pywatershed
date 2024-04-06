import pathlib as pl

import pytest
from utils_compare import compare_in_memory, compare_netcdfs

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.hydrology.prms_snow import PRMSSnow
from pywatershed.parameters import Parameters, PrmsParameters

# compare in memory (faster) or full output files? or both!
do_compare_output_files = False
do_compare_in_memory = True
rtol = atol = 1.0e-3

calc_methods = ("numpy", "numba")
params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(simulation):
    return Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )


@pytest.fixture(scope="function")
def discretization(simulation):
    dis_hru_file = simulation["dir"] / "parameters_dis_hru.nc"
    return Parameters.from_netcdf(dis_hru_file, encoding=False)


@pytest.fixture(scope="function", params=params)
def parameters(simulation, control, request):
    if request.param == "params_one":
        param_file = simulation["dir"] / control.options["parameter_file"]
        params = PrmsParameters.load(param_file)
    else:
        param_file = simulation["dir"] / "parameters_PRMSSnow.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.xfail
@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(
    simulation, control, discretization, parameters, tmp_path, calc_method
):
    tmp_path = pl.Path(tmp_path)

    # get the answer data
    comparison_var_names = [
        # 'ai',
        # 'albedo',
        # 'frac_swe',
        # 'freeh2o',
        # 'freeh2o_change',
        # 'freeh2o_prev',
        # 'iasw',
        # 'int_alb',
        "iso",
        # 'lso',
        # 'lst',
        # 'mso',
        # 'newsnow',
        # 'pk_def',
        # 'pk_den',
        # 'pk_depth',
        # 'pk_ice',
        # 'pk_ice_change',
        # 'pk_ice_prev',
        # 'pk_precip',
        # 'pk_temp',
        # 'pksv',
        # 'pkwater_ante',
        "pkwater_equiv",
        # 'pkwater_equiv_change',
        # 'pptmix_nopack',
        # 'pss',
        # 'pst',
        # 'salb',
        # 'scrv',
        # 'slst',
        "snow_evap",
        # 'snowcov_area',
        # 'snowcov_areasv',
        # 'snowmelt',
        # 'snsv',
        "tcal",
        "through_rain",
    ]

    output_dir = simulation["output_dir"]

    input_variables = {}
    for key in PRMSSnow.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = nc_path

    snow = PRMSSnow(
        control,
        discretization,
        parameters,
        **input_variables,
        budget_type="error",
        calc_method=calc_method,
    )

    if do_compare_output_files:
        nc_parent = tmp_path / simulation["name"]
        snow.initialize_netcdf(nc_parent)

    if do_compare_in_memory:
        answers = {}
        for var in comparison_var_names:
            var_pth = output_dir / f"{var}.nc"
            answers[var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    for istep in range(control.n_times):
        control.advance()
        snow.advance()
        snow.calculate(1.0)
        snow.output()
        if do_compare_in_memory:
            compare_in_memory(
                snow, answers, atol=atol, rtol=rtol, skip_missing_ans=True
            )

    snow.finalize()

    if do_compare_output_files:
        compare_netcdfs(
            comparison_var_names,
            tmp_path / simulation["name"],
            output_dir,
            atol=atol,
            rtol=rtol,
        )

    return
