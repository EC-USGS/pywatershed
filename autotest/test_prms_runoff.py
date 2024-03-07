import pathlib as pl

import pytest

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.hydrology.prms_runoff import PRMSRunoff
from pywatershed.parameters import Parameters, PrmsParameters
from utils_compare import compare_in_memory, compare_netcdfs

# compare in memory (faster) or full output files? or both!
do_compare_output_files = False
do_compare_in_memory = True
rtol = atol = 1.0e-10

calc_methods = ("numpy", "numba")
params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(simulation):
    control = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    control.options["netcdf_output_var_names"] = PRMSRunoff.get_variables()

    return control


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
        param_file = simulation["dir"] / "parameters_PRMSRunoff.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(
    simulation,
    control,
    discretization,
    parameters,
    tmp_path,
    calc_method,
):
    tmp_path = pl.Path(tmp_path)

    comparison_var_names = set(PRMSRunoff.get_variables())
    # TODO: this is hacky, improve the design
    if not control.options["dprst_flag"]:
        comparison_var_names = {
            vv for vv in comparison_var_names if "dprst" not in vv
        }

    # TODO: get rid of this exception
    comparison_var_names -= set(
        [
            "dprst_area_open",
        ]
    )

    output_dir = simulation["output_dir"]

    input_variables = {}
    for key in PRMSRunoff.get_inputs():
        nc_pth = output_dir / f"{key}.nc"
        input_variables[key] = nc_pth

    if do_compare_output_files:
        nc_parent = tmp_path / simulation["name"]
        control.options["netcdf_output_dir"] = nc_parent

    runoff = PRMSRunoff(
        control=control,
        discretization=discretization,
        parameters=parameters,
        **input_variables,
        budget_type="error",
        calc_method=calc_method,
    )

    if do_compare_output_files:
        runoff.initialize_netcdf()
        # test that init netcdf twice raises a warning
        with pytest.warns(UserWarning):
            runoff.initialize_netcdf()

    if do_compare_in_memory:
        answers = {}
        for var in comparison_var_names:
            var_pth = output_dir / f"{var}.nc"
            answers[var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    for istep in range(control.n_times):
        control.advance()
        runoff.advance()
        runoff.calculate(1.0)
        runoff.output()

        for var in answers.values():
            var.advance()
        if do_compare_in_memory:
            compare_in_memory(
                runoff, answers, atol=atol, rtol=rtol, skip_missing_ans=True
            )

    runoff.finalize()

    if do_compare_output_files:
        compare_netcdfs(
            comparison_var_names,
            tmp_path / simulation["name"],
            output_dir,
            atol=atol,
            rtol=rtol,
        )

    return
