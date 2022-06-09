import pathlib as pl

import numpy as np
import pytest

from pynhm.base.adapter import Adapter, AdapterNetcdf, adapter_factory
from pynhm.base.control import Control
from pynhm.constants import one
from pynhm.hydrology.PRMSEt import PRMSEt
from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


class TestPRMSEt:
    def test_init(self, domain, control, params, tmp_path):

        tmp_path = pl.Path(tmp_path)
        output_dir = domain["prms_output_dir"]

        et_inputs = {}
        for key in PRMSEt.get_inputs():
            nc_path = output_dir / f"{key}.nc"
            et_inputs[key] = adapter_factory(nc_path, key)

        et = PRMSEt(control=control, params=params, **et_inputs)

        # ---------------------------------
        # get the answer data
        comparison_vars = ["hru_actet"]

        # Read PRMS output into ans for comparison with pynhm results
        ans = {}
        for key in comparison_vars:
            nc_pth = output_dir / f"{key}.nc"
            ans[key] = adapter_factory(nc_pth, variable_name=key)

        all_success = True
        for istep in range(control.n_times):

            print(istep)
            control.advance()
            et.advance()
            et.calculate(1.0)

            # advance the answer, which is being read from a netcdf file
            for key, val in ans.items():
                val.advance()

            # make a comparison check with answer
            check = True
            failfast = True
            detailed = True
            if check:
                atol = 1.0e-5
                success = self.check_timestep_results(
                    et, istep, ans, atol, detailed, failfast
                )
                if not success:
                    all_success = False

        # check at the end and error if one or more steps didn't pass
        if not all_success:
            raise Exception("pynhm results do not match prms results")

        return

    @staticmethod
    def check_timestep_results(
        storageunit,
        istep,
        ans,
        atol,
        detailed=False,
        failfast=False,
    ):
        all_success = True
        for key in ans.keys():
            a1 = ans[key].current
            a2 = storageunit[key]
            success = np.isclose(a1, a2, atol=atol).all()
            if not success:
                all_success = False
                diff = a1 - a2
                diffmin = diff.min()
                diffmax = diff.max()
                if True:
                    print(f"time step {istep}")
                    print(f"output variable {key}")
                    print(f"prms   {a1.min()}    {a1.max()}")
                    print(f"pynhm  {a2.min()}    {a2.max()}")
                    print(f"diff   {diffmin}  {diffmax}")
                    if detailed:
                        idx = np.where(np.abs(diff) > atol)[0]
                        for i in idx:
                            print(
                                f"hru {i} prms {a1[i]} pynhm {a2[i]} diff {diff[i]}"
                            )
                if failfast:
                    raise (ValueError)
        return all_success
