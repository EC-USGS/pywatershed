import pathlib as pl

import pytest

from pywatershed.base.control import Control
from pywatershed.hydrology.SWBRootZone import SWBRootZone

from pywatershed import Parameters
from pywatershed.utils.netcdf_utils import NetCdfCompare


def test_init(domain, tmp_path):
    tmp_path = pl.Path(tmp_path)
    params = Parameters.from_yaml(domain["swb_param_file"])

    # TODO fix this. make it a data of DatasetDict.from_yaml()
    # nc_file = tmp_path / "swb_params_hru_1.nc"
    # params.to_netcdf(nc_file, use_xr=True)

    # Set information from the control file
    control = Control.load(domain["control_file"], params=params)

    # load csv files into dataframes
    output_dir = domain["prms_output_dir"]
    input_variables = {}
    for key in SWBRootZone.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = nc_path

    swb_rz = SWBRootZone(
        control,
        **input_variables,
        budget_type="warn",
    )

    nc_parent = tmp_path / domain["domain_name"]
    swb_rz.initialize_netcdf(nc_parent)

    # output_compare = {}
    # vars_compare = (
    #     "swb_rzres_flow",
    #     "swb_rzres_sink",
    #     "swb_rzres_stor",
    #     "ssr_to_swb_rz",
    #     "soil_to_swb_rz",
    # )
    # for key in SWBRootZone.get_variables():
    #     if key not in vars_compare:
    #         continue
    #     base_nc_path = output_dir / f"{key}.nc"
    #     compare_nc_path = tmp_path / domain["domain_name"] / f"{key}.nc"
    #     output_compare[key] = (base_nc_path, compare_nc_path)

    # print(f"base_nc_path: {base_nc_path}")
    # print(f"compare_nc_path: {compare_nc_path}")

    for istep in range(control.n_times):
        control.advance()
        swb_rz.advance()
        swb_rz.calculate(float(istep))
        swb_rz.output()

    swb_rz.finalize()

    asdf

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
