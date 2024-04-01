import pathlib as pl

import pytest
from utils_compare import compare_in_memory, compare_netcdfs

from pywatershed import Control, Parameters, PRMSGroundwater
from pywatershed.base.adapter import adapter_factory
from pywatershed.hydrology.prms_groundwater import has_prmsgroundwater_f
from pywatershed.parameters import PrmsParameters

# compare in memory (faster) or full output files?
do_compare_output_files = False
do_compare_in_memory = True
rtol = atol = 1.0e-13

calc_methods = ("numpy", "numba", "fortran")
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
        param_file = simulation["dir"] / "parameters_PRMSGroundwater.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(
    simulation, control, discretization, parameters, tmp_path, calc_method
):
    if not has_prmsgroundwater_f and calc_method == "fortran":
        pytest.skip(
            "PRMSGroundwater fortran code not available, skipping its test."
        )

    tmp_path = pl.Path(tmp_path)

    output_dir = simulation["output_dir"]
    input_variables = {}
    for key in PRMSGroundwater.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        # TODO: this is hacky for accommodating dprst_flag, improve the design
        # so people dont have to pass None for dead options.
        if not nc_path.exists():
            nc_path = None
        input_variables[key] = nc_path

    if do_compare_output_files:
        nc_parent = tmp_path / simulation["name"]
        control.options["netcdf_output_dir"] = nc_parent

    gw = PRMSGroundwater(
        control,
        discretization,
        parameters,
        **input_variables,
        budget_type="error",
        calc_method=calc_method,
    )

    if do_compare_output_files:
        gw.initialize_netcdf()

    if do_compare_in_memory:
        answers = {}
        for var in PRMSGroundwater.get_variables():
            var_pth = output_dir / f"{var}.nc"
            answers[var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    for istep in range(control.n_times):
        control.advance()
        gw.advance()
        gw.calculate(float(istep))
        gw.output()

        if do_compare_in_memory:
            compare_in_memory(gw, answers, atol=atol, rtol=rtol)

    gw.finalize()

    if do_compare_output_files:
        compare_netcdfs(
            PRMSGroundwater.get_variables(),
            tmp_path / simulation["name"],
            output_dir,
            atol=atol,
            rtol=rtol,
        )

    return
