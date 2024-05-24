import pathlib as pl

import numpy as np
import pytest

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.hydrology.prms_et import PRMSEt
from pywatershed.parameters import PrmsParameters


@pytest.fixture(scope="function")
def control(simulation):
    if simulation["name"] != "nhm":
        pytest.skip("Only test for nhm configuration")

    return Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )


@pytest.fixture(scope="function")
def params(simulation, control):
    param_file = simulation["dir"] / control.options["parameter_file"]
    return PrmsParameters.load(param_file)


@pytest.mark.domain
class TestPRMSEt:
    def test_init(self, simulation, control, params, tmp_path):
        tmp_path = pl.Path(tmp_path)
        output_dir = simulation["output_dir"]

        et_inputs = {}
        for key in PRMSEt.get_inputs():
            nc_path = output_dir / f"{key}.nc"
            et_inputs[key] = adapter_factory(nc_path, key, control)

        et = PRMSEt(
            control=control,
            discretization=None,
            parameters=params,
            budget_type="strict",
            **et_inputs,
        )

        # ---------------------------------
        # get the answer data
        comparison_vars = [
            # "potet",
            # "hru_impervevap",
            # "hru_intcpevap",
            # "snow_evap",
            # "dprst_evap_hru",
            # "perv_actet",
            "hru_actet",
        ]

        # Read PRMS output into ans for comparison with pywatershed results
        ans = {}
        for key in comparison_vars:
            nc_pth = output_dir / f"{key}.nc"
            ans[key] = adapter_factory(
                nc_pth, variable_name=key, control=control
            )

        all_success = True
        for istep in range(control.n_times):
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
            raise Exception("pywatershed results do not match prms results")

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
            success = np.isclose(a2, a1, atol=atol).all()
            if not success:
                all_success = False
                diff = a1 - a2
                diffmin = diff.min()
                diffmax = diff.max()
                if True:
                    print(f"time step {istep}")
                    print(f"output variable {key}")
                    print(f"prms   {a1.min()}    {a1.max()}")
                    print(f"pywatershed  {a2.min()}    {a2.max()}")
                    print(f"diff   {diffmin}  {diffmax}")
                    if detailed:
                        idx = np.where(np.abs(diff) > atol)[0]
                        for i in idx:
                            print(
                                f"hru {i} prms {a1[i]} pywatershed {a2[i]} "
                                f"diff {diff[i]}"
                            )
                if failfast:
                    raise (ValueError)
        return all_success
