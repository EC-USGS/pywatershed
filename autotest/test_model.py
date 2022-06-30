import pathlib as pl
from copy import deepcopy
from pprint import pprint

import numpy as np
import pytest

from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.base.model import Model
from pynhm.atmosphere.PRMSSolarGeometry import PRMSSolarGeometry
from pynhm.hydrology.PRMSCanopy import PRMSCanopy
from pynhm.hydrology.PRMSEt import PRMSEt
from pynhm.hydrology.PRMSGroundwater import PRMSGroundwater
from pynhm.hydrology.PRMSRunoff import PRMSRunoff
from pynhm.hydrology.PRMSSnow import PRMSSnow

from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


test_models = {
    # "et": [PRMSEt],
    # "et_canopy": [PRMSEt, PRMSCanopy],
    # "et_canopy_runoff": [PRMSEt, PRMSCanopy, PRMSRunoff],
    # "et_canopy_runoff_soil": [PRMSEt, PRMSCanopy, PRMSRunoff, PRMSSoilzone],
    "solar_snow": [PRMSSolarGeometry, PRMSSnow],
    "canopy_snow_runoff": [
        PRMSCanopy,
        PRMSSnow,
        # PRMSRunoff,
        PRMSEt,
        # PRMSSoilzone,
        # PRMSGroundwater,
    ],
}


@pytest.mark.parametrize(
    test_models.values(),
    ids=test_models.keys(),
)
def test_model(domain, control, params, components, tmp_path):

    tmp_path = pl.Path(tmp_path)
    output_dir = domain["prms_output_dir"]

    # TODO: Eliminate potet and other variables from being used
    model = Model(
        *components, control=control, params=params, input_dir=output_dir
    )

    # ---------------------------------
    # get the answer data
    comparison_vars_dict_all = {
        "PRMSCanopy": [
            "net_rain",
            "net_snow",
            "net_ppt",
            "intcp_stor",
            "intcp_evap",
            "hru_intcpstor",
            "hru_intcpevap",
            "potet",
        ],
        "PRMSEt": [
            "potet",
            "hru_impervevap",
            "hru_intcpevap",
            "snow_evap",
            "dprst_evap_hru",
            "perv_actet",
            "hru_actet",
        ],
        "PRMSRunoff": [
            "infil",
            "dprst_stor_hru",
            "hru_impervstor",
            "sroff",
            "dprst_evap_hru",
        ],
        "PRMSSolarGeometry": [],
        "PRMSSnow": [],
    }

    comparison_vars_dict = {}
    for cls in components:
        key = cls.__name__
        comparison_vars_dict[key] = comparison_vars_dict_all[key]

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
        model.advance()
        model.calculate()

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
                    model.components[unit_name],
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
            diff = a2 - a1
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
