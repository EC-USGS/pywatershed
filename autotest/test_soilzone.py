import pathlib as pl

import numpy as np
import pytest

from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.hydrology.PRMSSoilzone import PRMSSoilzone
from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


class TestPRMSSoilzone:
    def test_init(self, domain, control, params, tmp_path):
        tmp_path = pl.Path(tmp_path)
        output_dir = domain["prms_output_dir"]

        # get the answer data
        comparison_var_names = []

        ans = {}
        for key in comparison_var_names:
            nc_pth = output_dir / f"{key}.nc"
            ans[key] = adapter_factory(nc_pth, variable_name=key)

        # setup the soilzone
        input_variables = {}
        for key in PRMSSoilzone.get_inputs():
            nc_path = output_dir / f"{key}.nc"
            input_variables[key] = nc_path

        soil = PRMSSoilzone(control, params, **input_variables)
        all_success = True
        # for istep in range(control.n_times):
        for istep in range(10):

            control.advance()
            soil.advance()
            soil.calculate(float(istep))

            print("\n")
            print(control.current_time)

            # compare along the way
            atol = 1.0e-5
            for key, val in ans.items():
                val.advance()
            for key in ans.keys():
                a1 = ans[key].current
                a2 = soil[key]
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

            # if istep == 15:
            #     asdf

        soil.finalize()

        if not all_success:
            raise Exception("pynhm results do not match prms results")

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
