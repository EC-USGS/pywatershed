import pathlib as pl

import numpy as np
import pytest

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.hydrology.prms_canopy import PRMSCanopy
from pywatershed.hydrology.prms_et import PRMSEt
from pywatershed.hydrology.prms_runoff import PRMSRunoff
from pywatershed.parameters import PrmsParameters


@pytest.fixture(scope="function")
def control(simulation):
    domain_config = simulation["name"].split(":")[1]
    if domain_config != "nhm":
        pytest.skip("Only test for nhm configuration")

    return Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )


@pytest.fixture(scope="function")
def params(simulation, control):
    param_file = simulation["dir"] / control.options["parameter_file"]
    return PrmsParameters.load(param_file)


def test_et_can_runoff(simulation, control, params, tmp_path):
    tmp_path = pl.Path(tmp_path)
    output_dir = simulation["output_dir"]

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
        discretization=None,
        parameters=params,
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
        discretization=None,
        parameters=params,
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
        discretization=None,
        parameters=params,
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

    # Read PRMS output into ans for comparison with pywatershed results
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
                success = check_timestep_results(
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
        raise Exception("pywatershed results do not match prms results")

    return


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
