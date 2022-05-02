import pathlib as pl

import pytest

from pynhm.atmosphere.NHMSolarGeometry import NHMSolarGeometry
from pynhm.base.control import Control
from pynhm.hydrology.PRMSSnow import PRMSSnow
from pynhm.utils.netcdf_utils import NetCdfCompare
from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def solar_geom(domain):
    params = PrmsParameters.load(domain["param_file"])
    return NHMSolarGeometry(params)


class TestPRMSSnow:
    def test_init(self, domain, solar_geom, tmp_path):
        tmp_path = pl.Path(tmp_path)
        prms_params = PrmsParameters.load(domain["param_file"])

        # Set information from the control file
        control = Control.load(domain["control_file"])

        # load csv files into dataframes
        output_dir = domain["prms_output_dir"]
        input_variables = {}
        for key in PRMSSnow.get_inputs():
            if key == "soltab_horad_potsw":
                input_variables[key] = solar_geom["potential_sw_rad_flat"]
            else:
                nc_path = output_dir / f"{key}.nc"
                input_variables[key] = nc_path

        snow = PRMSSnow(control, prms_params, **input_variables)
        nc_parent = tmp_path / domain["domain_name"]
        snow.initialize_netcdf(nc_parent)

        output_compare = {}
        for key in PRMSSnow.get_variables():
            ans_path = output_dir / f"{key}.nc"
            result_path = tmp_path / domain["domain_name"] / f"{key}.nc"
            output_compare[key] = (result_path, ans_path)

        print(f"result_path: {result_path}")
        print(f"ans_path: {ans_path}")

        for istep in range(control.n_times):

            control.advance()
            snow.advance()
            snow.calculate(float(istep))

        #   gw.output()

        # gw.finalize()

        # assert_error = False
        # for key, (base, compare) in output_compare.items():
        #     success, diff = NetCdfCompare(base, compare).compare()
        #     if not success:
        #         print(f"comparison for {key} failed: maximum error {diff}")
        #         assert_error = True
        # assert not assert_error, "comparison failed"

        return
