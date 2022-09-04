import pathlib as pl
import shutil

import numpy as np
import pytest

import pynhm
from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.base.model import Model
from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


compare_to_prms521 = False
n_time_steps = 101
budget_type = None
test_models = {
    "nhm": [
        pynhm.PRMSSolarGeometry,
        pynhm.PRMSAtmosphere,
        pynhm.PRMSCanopy,
        pynhm.PRMSSnow,
        pynhm.PRMSRunoff,
        pynhm.PRMSSoilzone,
        pynhm.PRMSGroundwater,
        pynhm.PRMSChannel,
    ],
}


@pytest.mark.parametrize(
    "processes",
    test_models.values(),
    ids=test_models.keys(),
)
def test_model(domain, control, processes, tmp_path):

    tmp_path = pl.Path(tmp_path)
    output_dir = domain["prms_output_dir"]

    # setup input_dir with symlinked prms inputs and outputs
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    for ff in output_dir.resolve().glob("*.nc"):
        shutil.copy(ff, input_dir / ff.name)
    for ff in output_dir.parent.resolve().glob("*.nc"):
        shutil.copy(ff, input_dir / ff.name)

    # TODO: Eliminate potet and other variables from being used
    model = Model(
        *processes,
        control=control,
        input_dir=input_dir,
        budget_type=budget_type,
    )

    # ---------------------------------
    # get the answer data against PRMS5.2.1
    comparison_vars_dict_all = {
        "PRMSSolarGeometry": [],
        "PRMSAtmosphere": [
            "tmaxf",
            "tminf",
            "hru_ppt",
            "hru_rain",
            "hru_snow",
            "swrad",
            "potet",
        ],
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
        "PRMSSnow": [
            "iso",
            "pkwater_equiv",
            "snow_evap",
            "tcal",
        ],
        "PRMSRunoff": [
            "infil",
            "dprst_stor_hru",
            "hru_impervstor",
            "sroff",
            "dprst_evap_hru",
            # "dprst_vol_open",
            "dprst_seep_hru",
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
        "PRMSSoilzone": [
            "hru_actet",
            "perv_actet",
            "potet_lower",
            "potet_rechr",
            "recharge",
            "slow_flow",
            "slow_stor",
            "soil_lower",
            "soil_moist",
            "soil_moist_tot",
            "infil",
            "soil_rechr",
            "soil_to_gw",
            "soil_to_ssr",
            "ssr_to_gw",
            "ssres_flow",
            "ssres_in",
            "ssres_stor",
        ],
        "PRMSGroundwater": [
            "soil_to_gw",  # input
            "ssr_to_gw",  # input
            "dprst_seep_hru",  # input
            "gwres_flow",
            "gwres_sink",
            "gwres_stor",
        ],
        "PRMSChannel": [
            "seg_lateral_inflow",
            "seg_upstream_inflow",
            "seg_outflow",
        ],
    }

    tol = {
        "PRMSSolarGeometry": 1.0e-5,
        "PRMSAtmosphere": 1.0e-4,
        "PRMSCanopy": 1.0e-5,
        "PRMSSnow": 5e-2,
        "PRMSRunoff": 1.0e-5,
        "PRMSSoilzone": 1.0e-4,  ##
        "PRMSGroundwater": 1.0e-2,
        "PRMSEt": 1.0e-5,
        "PRMSChannel": 1.0e-5,
    }

    comparison_vars_dict = {}
    for cls in processes:
        key = cls.__name__
        comparison_vars_dict[key] = comparison_vars_dict_all[key]

    # Read PRMS output into ans for comparison with pynhm results
    ans = {key: {} for key in comparison_vars_dict.keys()}
    for unit_name, var_names in comparison_vars_dict.items():
        for vv in var_names:
            if vv in ["tmin", "tmax", "prcp"]:
                nc_pth = input_dir.parent / f"{vv}.nc"
            else:
                nc_pth = input_dir / f"{vv}.nc"
            ans[unit_name][vv] = adapter_factory(
                nc_pth, variable_name=vv, control=control
            )

    # ---------------------------------
    # itimestep: process: variable: mean
    regression_ans = {
        9: {
            "PRMSChannel": {
                "seg_outflow": {
                    "drb_2yr": 1152.9343097732165,
                    "hru_1": 13.37194463307104,
                    "ucb_2yr": 1694.47282087809,
                },
            },
        },
        99: {
            "PRMSChannel": {
                "seg_outflow": {
                    "drb_2yr": 1785.102361378249,
                    "hru_1": 14.614890968139237,
                    "ucb_2yr": 710.4918051193282,
                },
            },
        },
    }

    all_success = True
    fail_prms_compare = False
    fail_regression = False
    # for istep in range(control.n_times):

    for istep in range(n_time_steps):

        # print(istep)
        model.advance()
        model.calculate()

        # PRMS5 answers
        # advance the answer, which is being read from a netcdf file
        for unit_name, var_names in ans.items():
            for vv in var_names:
                ans[unit_name][vv].advance()

        # make a comparison check with answer
        failfast = True
        detailed = True

        # compare_to_prms521 is a golbal variable
        if compare_to_prms521:
            for unit_name in ans.keys():
                success = check_timestep_results(
                    model.processes[unit_name],
                    istep,
                    ans[unit_name],
                    tol[unit_name],
                    detailed,
                    failfast,
                )
                if not success:
                    fail_prms_compare = True
                    all_success = False
                    if failfast:
                        assert False, "PRMS comparison failfast"

        # Regression checking
        if istep in regression_ans:
            for pp, var_ans in regression_ans[istep].items():
                for vv, aa in var_ans.items():
                    result = model.processes[pp][vv].mean()
                    reg_ans = aa[domain["domain_name"]]
                    if not reg_ans:
                        print(
                            f"\nreg_ans: [{istep}][{pp}][{vv}] mean: {result}"
                        )
                        success = False
                    else:
                        success = (abs(result - reg_ans) / abs(reg_ans)) < 1e-5

                        if not success:
                            fail_regression = True
                            all_success = False
                            if failfast:
                                assert False, "Regression failffast"

    # check at the end and error if one or more steps didn't pass
    if not all_success:
        if fail_prms_compare and fail_regression:
            msg = (
                "pynhm results both failed regression test and comparison "
                "with prms5.2.1"
            )
        elif fail_prms_compare:
            msg = "pynhm results failed comparison with prms5.2.1"
        elif fail_regression:
            msg = "pynhm results failed regression"
        else:
            assert False, "this should not be possible"

        raise Exception(msg)

    return


def check_timestep_results(
    storageunit,
    istep,
    ans,
    tol,
    detailed=False,
    failfast=False,
):
    # print(storageunit)
    all_success = True
    for key in ans.keys():
        # print(key)
        a1 = ans[key].current
        a2 = storageunit[key]
        success_a = np.isclose(a2, a1, atol=tol, rtol=0.0)
        success_r = np.isclose(a2, a1, atol=0.0, rtol=tol)
        success = success_a | success_r
        if not success.all():
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
                    idx = np.where(~success)
                    for i in idx:
                        print(
                            f"hru {i} prms {a1[i]} pynhm {a2[i]} "
                            f"diff {diff[i]}"
                        )
            if failfast:
                raise (ValueError)
    return all_success
