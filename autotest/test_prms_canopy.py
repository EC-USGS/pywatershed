import pathlib as pl

import pytest
from utils_compare import compare_in_memory, compare_netcdfs

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.hydrology.prms_canopy import PRMSCanopy
from pywatershed.parameters import Parameters, PrmsParameters

# compare in memory (faster) or full output files? or both!
do_compare_output_files = True
do_compare_in_memory = False
rtol = atol = 1e-12

calc_methods = ("numpy", "numba")
params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(simulation):
    ctl = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    del ctl.options["netcdf_output_dir"]
    del ctl.options["netcdf_output_var_names"]
    return ctl


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
        param_file = simulation["dir"] / "parameters_PRMSCanopy.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(
    simulation, control, discretization, parameters, tmp_path, calc_method
):
    tmp_path = pl.Path(tmp_path)

    # get the answer data
    comparison_var_names = [
        "net_rain",
        "net_snow",
        "net_ppt",
        "intcp_stor",
        "intcp_evap",
        "hru_intcpstor",
        "hru_intcpevap",
        "intcp_changeover",
        "intcp_form",
        # "intcp_transp_on",  # not a prms variable
        "hru_intcpstor_change",
        "hru_intcpstor_old",
    ]

    output_dir = simulation["output_dir"]

    input_variables = {}
    for key in PRMSCanopy.get_inputs():
        nc_pth = output_dir / f"{key}.nc"
        input_variables[key] = nc_pth

    canopy = PRMSCanopy(
        control=control,
        discretization=discretization,
        parameters=parameters,
        **input_variables,
        budget_type="error",
        calc_method=calc_method,
    )

    if do_compare_output_files:
        nc_parent = tmp_path / simulation["name"].replace(":", "_")
        canopy.initialize_netcdf(nc_parent)

    if do_compare_in_memory:
        answers = {}
        for var in comparison_var_names:
            var_pth = output_dir / f"{var}.nc"
            answers[var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    for istep in range(control.n_times):
        control.advance()
        canopy.advance()
        canopy.calculate(1.0)
        canopy.output()
        if do_compare_in_memory:
            for var in answers.values():
                var.advance()
            compare_in_memory(
                canopy,
                answers,
                atol=atol,
                rtol=rtol,
                skip_missing_ans=True,
                fail_after_all_vars=False,
            )

    canopy.finalize()

    if do_compare_output_files:
        compare_netcdfs(
            comparison_var_names,
            tmp_path / simulation["name"].replace(":", "_"),
            output_dir,
            atol=atol,
            rtol=rtol,
        )

    return
