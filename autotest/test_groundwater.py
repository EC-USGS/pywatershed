import pathlib as pl

import pytest

from pywatershed.base.control import Control
from pywatershed.hydrology.PRMSGroundwater import (
    PRMSGroundwater,
    has_prmsgroundwater_f,
)
from pywatershed.parameters import PrmsParameters
from pywatershed.utils.netcdf_utils import NetCdfCompare

calc_methods = ("numpy", "numba")
if has_prmsgroundwater_f:
    calc_methods += ("fortran",)


@pytest.mark.parametrize("calc_method", calc_methods)
def test_run(domain, tmp_path, calc_method):
    tmp_path = pl.Path(tmp_path)
    params = PrmsParameters.load(domain["param_file"])

    # Set information from the control file
    control = Control.load(domain["control_file"], params=params)

    # load csv files into dataframes
    output_dir = domain["prms_output_dir"]
    input_variables = {}
    for key in PRMSGroundwater.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = nc_path

    gw = PRMSGroundwater(
        control,
        **input_variables,
        budget_type="error",
        calc_method=calc_method,
    )
    nc_parent = tmp_path / domain["domain_name"]
    gw.initialize_netcdf(nc_parent)

    output_compare = {}
    vars_compare = (
        "gwres_flow",
        "gwres_sink",
        "gwres_stor",
        "ssr_to_gw",
        "soil_to_gw",
    )
    for key in PRMSGroundwater.get_variables():
        if key not in vars_compare:
            continue
        base_nc_path = output_dir / f"{key}.nc"
        compare_nc_path = tmp_path / domain["domain_name"] / f"{key}.nc"
        output_compare[key] = (base_nc_path, compare_nc_path)

    print(f"base_nc_path: {base_nc_path}")
    print(f"compare_nc_path: {compare_nc_path}")

    for istep in range(control.n_times):
        control.advance()

        gw.advance()

        gw.calculate(float(istep))

        gw.output()

    gw.finalize()

    assert_error = False
    for key, (base, compare) in output_compare.items():
        success, diff = NetCdfCompare(base, compare).compare()
        if not success:
            print(
                f"comparison for {key} failed: "
                + f"maximum error {diff[key][0]} "
                + f"(maximum allowed error {diff[key][1]}) "
                + f"in column {diff[key][2]}"
            )
            assert_error = True
    assert not assert_error, "comparison failed"

    return
