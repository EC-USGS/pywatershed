import pathlib as pl

import pytest

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.base.parameters import Parameters
from pywatershed.hydrology.prms_channel import PRMSChannel, has_prmschannel_f
from pywatershed.parameters import PrmsParameters
from utils_compare import compare_in_memory, compare_netcdfs

# compare in memory (faster) or full output files?
do_compare_output_files = True
do_compare_in_memory = True
rtol = atol = 1.0e-7

fail_fast = False

calc_methods = ("numpy", "numba", "fortran")
params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(domain):
    return Control.load_prms(domain["control_file"], warn_unused_options=False)


@pytest.fixture(scope="function")
def discretization(domain):
    dis_hru_file = domain["dir"] / "parameters_dis_hru.nc"
    dis_seg_file = domain["dir"] / "parameters_dis_seg.nc"
    dis = Parameters.merge(
        Parameters.from_netcdf(dis_hru_file, encoding=False),
        Parameters.from_netcdf(dis_seg_file, encoding=False),
    )
    return dis


@pytest.fixture(scope="function", params=params)
def parameters(domain, request):
    if request.param == "params_one":
        params = PrmsParameters.load(domain["param_file"])
    else:
        param_file = domain["dir"] / "parameters_PRMSChannel.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(
    domain, control, discretization, parameters, tmp_path, calc_method
):
    if not has_prmschannel_f and calc_method == "fortran":
        pytest.skip(
            "PRMSChannel fortran code not available, skipping its test."
        )

    tmp_path = pl.Path(tmp_path)

    # load csv files into dataframes
    output_dir = domain["prms_output_dir"]
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
        nc_parent = tmp_path / domain["domain_name"]
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
            compare_in_memory(channel, answers, atol=atol, rtol=rtol)

    channel.finalize()

    if do_compare_output_files:
        compare_netcdfs(
            PRMSChannel.get_variables(),
            tmp_path / domain["domain_name"],
            output_dir,
            atol=atol,
            rtol=rtol,
        )

    return
