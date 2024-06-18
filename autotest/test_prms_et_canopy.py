import pathlib as pl

import numpy as np
import pytest

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.hydrology.prms_canopy import PRMSCanopy
from pywatershed.hydrology.prms_et import PRMSEt
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


@pytest.mark.domain
def test_et(simulation, control, params, tmp_path):
    tmp_path = pl.Path(tmp_path)
    output_dir = simulation["output_dir"]

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
        control=control,
        discretization=None,
        parameters=params,
        budget_type="error",
        **et_inputs,
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
        discretization=None,
        parameters=params,
        budget_type="error",
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

    # Read PRMS output into ans for comparison with pywatershed results
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
