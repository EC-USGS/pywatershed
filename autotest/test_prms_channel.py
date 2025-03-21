import pathlib as pl

import pytest
from utils_compare import compare_in_memory, compare_netcdfs

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.base.parameters import Parameters
from pywatershed.hydrology.prms_channel import PRMSChannel
from pywatershed.parameters import PrmsParameters

# compare in memory (faster) or full output files? or both!
do_compare_output_files = True
do_compare_in_memory = True
rtol = atol = 1.0e-7

calc_methods = ("numpy", "numba")
params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(simulation):
    if "obsin" in simulation["name"]:
        pytest.skip("Not testing passthrough flow graph for drb_2yr:nhm_obsin")

    ctl = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    if ctl.options["streamflow_module"] == "strmflow":
        pytest.skip(
            f"PRMSChannel not present in simulation {simulation['name']}"
        )
    del ctl.options["netcdf_output_dir"]
    del ctl.options["netcdf_output_var_names"]

    return ctl


@pytest.fixture(scope="function")
def discretization(simulation):
    dis_hru_file = simulation["dir"] / "parameters_dis_hru.nc"
    dis_seg_file = simulation["dir"] / "parameters_dis_seg.nc"
    dis = Parameters.merge(
        Parameters.from_netcdf(dis_hru_file, encoding=False),
        Parameters.from_netcdf(dis_seg_file, encoding=False),
    )
    return dis


@pytest.fixture(scope="function", params=params)
def parameters(simulation, control, request):
    if request.param == "params_one":
        param_file = simulation["dir"] / control.options["parameter_file"]
        params = PrmsParameters.load(param_file)
    else:
        param_file = simulation["dir"] / "parameters_PRMSChannel.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(
    simulation, control, discretization, parameters, tmp_path, calc_method
):
    tmp_path = pl.Path(tmp_path)
    output_dir = simulation["output_dir"]

    input_variables = {}
    for key in PRMSChannel.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = nc_path

    channel = PRMSChannel(
        control,
        discretization,
        parameters,
        **input_variables,
        budget_type="error",
        calc_method=calc_method,
    )

    if do_compare_output_files:
        nc_parent = tmp_path / simulation["name"].replace(":", "_")
        channel.initialize_netcdf(nc_parent)
        # test that init netcdf twice raises a warning
        with pytest.warns(UserWarning):
            channel.initialize_netcdf(nc_parent)

    if do_compare_in_memory:
        answers = {}
        for var in PRMSChannel.get_variables():
            var_pth = output_dir / f"{var}.nc"
            answers[var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    for istep in range(control.n_times):
        control.advance()
        channel.advance()
        channel.calculate(float(istep))
        channel.output()
        if do_compare_in_memory:
            for var in answers.values():
                var.advance()
            compare_in_memory(channel, answers, atol=atol, rtol=rtol)

    channel.finalize()

    if do_compare_output_files:
        compare_netcdfs(
            PRMSChannel.get_variables(),
            tmp_path / simulation["name"].replace(":", "_"),
            output_dir,
            atol=atol,
            rtol=rtol,
        )

    return
