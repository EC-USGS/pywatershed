import pathlib as pl

import numpy as np
import pytest

from pynhm.atmosphere.PRMSBoundaryLayer import PRMSBoundaryLayer
from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


class TestPRMSBoundaryLayer:
    def test_init(self, domain, control, tmp_path):

        output_dir = domain["prms_output_dir"]
        cbh_dir = domain["cbh_inputs"]["prcp"].parent.resolve()

        # get the answer data
        comparison_var_names = [
            "tmaxf",
            "tminf",
            "hru_ppt",
            "hru_rain",
            "hru_snow",
            "swrad",
            "potet",
        ]
        ans = {}
        for key in comparison_var_names:
            nc_pth = output_dir / f"{key}.nc"
            ans[key] = adapter_factory(
                nc_pth, variable_name=key, control=control
            )

        input_variables = {}
        for key in PRMSBoundaryLayer.get_inputs():
            nc_pth = cbh_dir / f"{key}.nc"
            input_variables[key] = nc_pth

        atm = PRMSBoundaryLayer(
            control=control,
            **input_variables,
            budget_type="strict",
            netcdf_output_dir=tmp_path,
        )

        all_success = True
        for istep in range(control.n_times):
            control.advance()
            atm.advance()
            atm.calculate(1.0)
            # print(atm.budget)

            # compare along the way
            atol = 1.0e-5
            for key, val in ans.items():
                val.advance()
            for key in ans.keys():
                a1 = ans[key].current
                a2 = atm[key]
                success = np.isclose(a2, a1, atol=atol).all()
                if not success:
                    all_success = False
                    diff = a1 - a2
                    diffmin = diff.min()
                    diffmax = diff.max()
                    print(f"time step {istep}")
                    print(f"output variable {key}")
                    print(f"prms   {a1.min()}    {a1.max()}")
                    print(f"pynhm  {a2.min()}    {a2.max()}")
                    print(f"diff   {diffmin}  {diffmax}")

        atm.finalize()

        if not all_success:
            raise Exception("pynhm results do not match prms results")
