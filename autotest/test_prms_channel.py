import pathlib as pl

import pytest

from pywatershed.base.control import Control
from pywatershed.base.parameters import Parameters
from pywatershed.hydrology.prms_channel import PRMSChannel, has_prmschannel_f
from pywatershed.parameters import PrmsParameters
from pywatershed.utils.netcdf_utils import NetCdfCompare

fail_fast = False

calc_methods = ("numpy", "numba", "fortran")
params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


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
    nc_parent = tmp_path / domain["domain_name"]
    channel.initialize_netcdf(nc_parent)

    for istep in range(control.n_times):
        control.advance()

        channel.advance()

        channel.calculate(float(istep))

        channel.output()

    channel.finalize()

    output_compare = {}

    for key in PRMSChannel.get_variables():
        base_nc_path = output_dir / f"{key}.nc"
        compare_nc_path = tmp_path / domain["domain_name"] / f"{key}.nc"
        # PRMS does not output the storage change in the channel
        if not base_nc_path.exists():
            continue
        output_compare[key] = (base_nc_path, compare_nc_path)

    assert_error = False
    for key, (base, compare) in output_compare.items():
        print(f"\nbase_nc_path: {base}")
        print(f"compare_nc_path: {compare}")
        success, diff = NetCdfCompare(base, compare).compare()
        if not success:
            print(
                f"comparison for {key} failed: "
                + f"maximum error {diff[key][0]} "
                + f"(maximum allowed error {diff[key][1]}) "
                + f"in column {diff[key][2]}"
            )
            assert_error = True
            if fail_fast:
                assert False

        else:
            print(f"comparison for {key} passed")

    assert not assert_error, "comparison failed"

    return
