import pathlib as pl
from warnings import warn

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
            "transp_on",
            "tmaxc",
            "tavgc",
            "tminc",
            "prmx",
            "pptmix",
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
            for key, val in ans.items():
                val.advance()
                
            for key in ans.keys():
                a1 = ans[key].current
                a2 = atm[key]

                tol = 5e-6
                if key == 'swrad':
                    tol = 5e-4
                    warn(f"using tol = {tol} for variable {key}")
                if key == 'tavgc':
                    tol = 1e-5
                    warn(f"using tol = {tol} for variable {key}")

                success_a = np.allclose(a2, a1, atol=tol, rtol=0.00)
                success_r = np.allclose(a2, a1, atol=0.00, rtol=tol)
                if (not success_a) and (not success_r):
                    diff = a2 - a1
                    diffratio = abs(diff / a2)
                    if (diffratio < 1e-6).all():
                        success = True
                        continue
                    all_success = False
                    diffmin = diff.min()
                    diffmax = diff.max()
                    abs_diff = abs(diff)
                    absdiffmax = abs_diff.max()
                    wh_absdiffmax = np.where(abs_diff)[0]
                    print(f"time step {istep}")
                    print(f"output variable {key}")
                    print(f"prms   {a1.min()}    {a1.max()}")
                    print(f"pynhm  {a2.min()}    {a2.max()}")
                    print(f"diff   {diffmin}  {diffmax}")
                    print(f"absdiffmax  {absdiffmax}")
                    print(f"wh_absdiffmax  {wh_absdiffmax}")
                    asdf

        atm.finalize()

        if not all_success:
            raise Exception("pynhm results do not match prms results")
