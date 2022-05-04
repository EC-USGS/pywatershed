import pathlib as pl

import numpy as np
import pytest

from pynhm.atmosphere.NHMSolarGeometry import NHMSolarGeometry
from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.hydrology.PRMSSnow import PRMSSnow
from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def solar_geom(domain):
    params = PrmsParameters.load(domain["param_file"])
    return NHMSolarGeometry(params)


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


class TestPRMSSnow:
    def test_init(self, domain, control, params, solar_geom, tmp_path):
        tmp_path = pl.Path(tmp_path)
        output_dir = domain["prms_output_dir"]

        # get the answer data
        comparison_var_names = [
            # "ai",
            # "albedo",
            # "frac_swe",
            # "freeh2o",
            # "iasw",
            # "int_alb",
            # "iso",
            # "lso",
            # "lst",
            # "mso",
            # "pk_def",
            # "pk_den",
            # "pk_depth",
            # "pk_ice",
            # "pk_precip",
            # "pk_temp",
            # "pksv",
            "pkwater_ante",
            "pkwater_equiv",
            # "pptmix_nopack",
            # "pss",
            # "pst",
            # "salb",
            # "scrv",
            # "slst",
            # "snow_evap",
            # "snowcov_area",
            # "snowcov_areasv",
            # "snowmelt",
            # "snsv",
            "tcal",
        ]

        ans = {}
        for key in comparison_var_names:
            nc_pth = output_dir / f"{key}.nc"
            ans[key] = adapter_factory(nc_pth, variable_name=key)

        # setup the snow
        input_variables = {}
        for key in PRMSSnow.get_inputs():
            if key == "soltab_horad_potsw":
                input_variables[key] = solar_geom["potential_sw_rad_flat"]
            else:
                nc_path = output_dir / f"{key}.nc"
                input_variables[key] = nc_path

        snow = PRMSSnow(control, params, **input_variables)

        all_success = True
        for istep in range(control.n_times):

            control.advance()
            snow.advance()
            snow.calculate(float(istep))

            print("\n")
            print(control.current_time)

            # compare along the way
            atol = 1.0e-5
            for key, val in ans.items():
                val.advance()
            for key in ans.keys():
                a1 = ans[key].current
                a2 = snow[key]
                success = np.isclose(a1, a2, atol=atol).all()
                if not success:
                    all_success = False
                    diff = a1 - a2
                    diffmin = diff.min()
                    diffmax = diff.max()
                    print(f"time step {istep}")
                    print(f"output variable: {key}")
                    print(f"prms   {a1.min()}    {a1.max()}")
                    print(f"pynhm  {a2.min()}    {a2.max()}")
                    print(f"diff   {diffmin}  {diffmax}")

            if istep == 15:
                asdf

        cnp.finalize()

        if not all_success:
            raise Exception("pynhm results do not match prms results")

        asdf
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
