import pathlib as pl

import numpy as np
import pytest

from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.hydrology.PRMSCanopy import PRMSCanopy
from pynhm.hydrology.PRMSEt import PRMSEt
from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


class TestPRMSEt:
    def test_init(self, domain, control, params, tmp_path):

        tmp_path = pl.Path(tmp_path)
        output_dir = domain["prms_output_dir"]

        # Variable wiring overview
        # et:
        #   potet: TO canopy
        #   hru_intcpevap: FROM canopy
        #
        # canopy:
        #   potet: FROM et
        #   hru_intcpevap: TO et

        # et
        et_inputs = {}
        for key in PRMSEt.get_inputs():
            if key == "hru_intcpevap":
                et_inputs[key] = None
            else:
                nc_path = output_dir / f"{key}.nc"
                et_inputs[key] = adapter_factory(nc_path, key, control)

        et = PRMSEt(
            control=control, params=params, budget_type="strict", **et_inputs
        )

        # canopy
        canopy_inputs = {}
        for key in PRMSCanopy.get_inputs():
            if key == "potet":
                canopy_inputs[key] = None
            else:
                nc_pth = output_dir / f"{key}.nc"
                canopy_inputs[key] = adapter_factory(nc_pth, key, control)

        canopy = PRMSCanopy(
            control=control,
            params=params,
            budget_type="strict",
            **canopy_inputs,
        )

        # wire up shared variables
        # et
        et.set_input_to_adapter(
            "hru_intcpevap",
            adapter_factory(canopy.hru_intcpevap),
        )

        # canopy
        canopy.set_input_to_adapter("potet", adapter_factory(et.potet))

        # get the answer data
        comparison_vars_dict = {
            "et": [
                "hru_intcpevap",
                "hru_actet",
                "potet",
            ],
            "canopy": [
                "hru_intcpstor",
                "hru_intcpevap",
            ],
        }

        # Read PRMS output into ans for comparison with pynhm results
        ans = {key: {} for key in comparison_vars_dict.keys()}
        for unit_name, var_names in comparison_vars_dict.items():
            for vv in var_names:
                nc_pth = output_dir / f"{vv}.nc"
                ans[unit_name][vv] = adapter_factory(nc_pth, vv, control)

        all_success = True
        for istep in range(control.n_times):

            # advance
            control.advance()
            canopy.advance()
            et.advance()

            # calculate
            canopy.calculate(1.0)

            et.calculate(1.0)

            # check
            # advance the answer, which is being read from a netcdf file
            for unit_name, var_names in ans.items():
                for vv in var_names:
                    ans[unit_name][vv].advance()

            # make a comparison check with answer
            check = True
            failfast = True
            detailed = True
            if check:
                atol = 1.0e-5
                for unit_name in ans.keys():
                    success = self.check_timestep_results(
                        locals()[unit_name],
                        istep,
                        ans[unit_name],
                        atol,
                        detailed,
                        failfast,
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
                    print(f"pynhm  {a2.min()}    {a2.max()}")
                    print(f"diff   {diffmin}  {diffmax}")
                    if detailed:
                        idx = np.where(np.abs(diff) > atol)[0]
                        for i in idx:
                            print(
                                f"hru {i} prms {a1[i]} pynhm {a2[i]} "
                                f"diff {diff[i]}"
                            )
                if failfast:
                    raise (ValueError)
        return all_success
