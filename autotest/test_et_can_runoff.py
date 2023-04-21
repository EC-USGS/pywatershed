import pathlib as pl

import numpy as np
import pytest

from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.hydrology.PRMSCanopy import PRMSCanopy
from pynhm.hydrology.PRMSEt import PRMSEt
from pynhm.hydrology.PRMSRunoff import PRMSRunoff
from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


class TestPRMSCanopyRunoffDomain:
    def test_init(self, domain, control, tmp_path):
        tmp_path = pl.Path(tmp_path)
        output_dir = domain["prms_output_dir"]

        # Variable wiring overview
        # et:
        #   potet: TO canopy and runoff
        #   hru_intcpevap: FROM canopy
        #   hru_impervevap: FROM runoff
        #   dprst_evap_hru: FROM runoff
        #   snow_evap: FROM file, common with runnoff
        #
        # canopy:
        #   potet: FROM et
        #   hru_intcpevap: TO et and runoff
        #   net_rain: TO runoff
        #   net_snow: TO runoff
        #   net_ppt: TO runoff
        #
        # runoff:
        #   potet: FROM et
        #   hru_intcpevap: FROM canopy
        #   net_rain: FROM canopy
        #   net_snow: FROM canopy
        #   net_ppt: FROM canopy
        #   snow_evap: FROM file common with et
        #   dprst_evap_hru: TO et
        #
        # common inputs:
        #   snow_evap: runoff, et
        #

        # TODO: Eliminate potet and other variables from being used

        common_inputs = {}
        for key in ["snow_evap"]:
            nc_path = output_dir / f"{key}.nc"
            common_inputs[key] = adapter_factory(nc_path, key, control)

        # ---------------------------------
        # init ET
        et_inputs = {}
        for key in PRMSEt.get_inputs():
            if key in [
                "hru_intcpevap",
                "hru_impervevap",
                "dprst_evap_hru",
                "snow_evap",
            ]:
                et_inputs[key] = None
            else:
                nc_path = output_dir / f"{key}.nc"
                et_inputs[key] = adapter_factory(nc_path, key, control)

        et = PRMSEt(
            control=control,
            budget_type="error",
            **et_inputs,
        )

        # ---------------------------------
        # instantiate canopy
        canopy_inputs = {}
        for key in PRMSCanopy.get_inputs():
            if key in ["potet"]:
                canopy_inputs[key] = None
            else:
                nc_pth = output_dir / f"{key}.nc"
                canopy_inputs[key] = nc_pth

        canopy = PRMSCanopy(
            control=control,
            budget_type="error",
            **canopy_inputs,
        )

        # ---------------------------------
        # instantiate runoff
        runoff_inputs = {}
        for key in PRMSRunoff.get_inputs():
            # These are the variables from et and canopy
            if key in [
                "potet",
                "hru_intcpevap",
                "net_ppt",
                "net_rain",
                "net_snow",
                "snow_evap",
            ]:
                runoff_inputs[key] = None

            else:
                nc_pth = output_dir / f"{key}.nc"
                if "soil_moist" in str(nc_pth):
                    nc_pth = output_dir / "soil_moist_prev.nc"
                runoff_inputs[key] = nc_pth

        runoff = PRMSRunoff(
            control=control,
            **runoff_inputs,
            budget_type=None,
        )

        # ---------------------------------
        # wire up shared variables
        # et
        et.set_input_to_adapter(
            "hru_intcpevap",
            adapter_factory(canopy.hru_intcpevap, control=control),
        )
        et.set_input_to_adapter(
            "hru_impervevap",
            adapter_factory(runoff.hru_impervevap, control=control),
        )
        et.set_input_to_adapter(
            "dprst_evap_hru",
            adapter_factory(runoff.dprst_evap_hru, control=control),
        )
        et.set_input_to_adapter(
            "snow_evap", adapter_factory(common_inputs["snow_evap"])
        )

        # canopy
        canopy.set_input_to_adapter(
            "potet", adapter_factory(et.potet, control=control)
        )

        # runoff
        runoff.set_input_to_adapter(
            "potet", adapter_factory(et.potet, control=control)
        )
        runoff.set_input_to_adapter(
            "snow_evap", adapter_factory(common_inputs["snow_evap"])
        )
        runoff.set_input_to_adapter(
            "net_ppt", adapter_factory(canopy.net_ppt, control=control)
        )
        runoff.set_input_to_adapter(
            "net_rain", adapter_factory(canopy.net_rain, control=control)
        )
        runoff.set_input_to_adapter(
            "net_snow", adapter_factory(canopy.net_snow, control=control)
        )
        runoff.set_input_to_adapter(
            "hru_intcpevap",
            adapter_factory(canopy.hru_intcpevap, control=control),
        )

        # ---------------------------------
        # get the answer data
        comparison_vars_dict = {
            "canopy": [
                "net_rain",
                "net_snow",
                "net_ppt",
                "intcp_stor",
                "intcp_evap",
                "hru_intcpstor",
                "hru_intcpevap",
                "potet",
            ],
            "et": [
                "potet",
                "hru_impervevap",
                "hru_intcpevap",
                "snow_evap",
                "dprst_evap_hru",
                "perv_actet",
                "hru_actet",
            ],
            "runoff": [
                "infil",
                "dprst_stor_hru",
                "hru_impervstor",
                "sroff",
                "dprst_evap_hru",
            ],
        }

        # Read PRMS output into ans for comparison with pynhm results
        ans = {key: {} for key in comparison_vars_dict.keys()}
        for unit_name, var_names in comparison_vars_dict.items():
            for vv in var_names:
                nc_pth = output_dir / f"{vv}.nc"
                ans[unit_name][vv] = adapter_factory(
                    nc_pth, variable_name=vv, control=control
                )

        all_success = True
        for istep in range(control.n_times):
            # print(istep)
            control.advance()
            canopy.advance()
            runoff.advance()
            et.advance()

            canopy.calculate(1.0)
            runoff.calculate(1.0)
            et.calculate(1.0)

            # advance the answer, which is being read from a netcdf file
            for unit_name, var_names in ans.items():
                for vv in var_names:
                    ans[unit_name][vv].advance()

            # make a comparison check with answer
            check = True
            failfast = False
            detailed = True
            if check:
                atol = 1.0e-4
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
