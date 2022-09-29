import pathlib as pl

import numpy as np
import pytest
import xarray as xr

import pynhm
from pynhm.base.control import Control
from pynhm.base.model import Model
from pynhm.utils.parameters import PrmsParameters

# test for a few timesteps a model with both unit/cell and global balance
# budgets


# This probably dosent need varied over domain
@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


n_time_steps = 10
model_procs = [pynhm.PRMSCanopy, pynhm.PRMSChannel]

# Things to parameterize
# optional variables to processes
# optional variables to budgets
budget_type = "error"  # vary this? shouldnt matter

check_vars = {
    "PRMSCanopy": [
        "hru_intcpstor",
        "hru_intcpstor_change",
        "net_snow",
    ],
    "PRMSChannel": [
        "channel_outflow_vol",
        "seg_stor_change",
        "seg_upstream_inflow",
    ],
}

budget_sum_vars_all = ["inputs_sum", "outputs_sum", "storage_changes_sum"]
check_budget_sum_vars_params = [False, True, "some"]


@pytest.mark.filterwarnings("ignore:Budget for")
@pytest.mark.parametrize(
    "budget_sum_param",
    check_budget_sum_vars_params,
    ids=[str(ii) for ii in check_budget_sum_vars_params],
)
def test_process_budgets(domain, control, tmp_path, budget_sum_param):

    tmp_dir = pl.Path(tmp_path)
    # print(tmp_dir)

    # Deal with parameter around what budget sum vars to write and check
    if budget_sum_param == "some":
        check_budget_sum_vars = budget_sum_vars_all[0:-1]
        budget_args = {"write_sum_vars": check_budget_sum_vars}
    elif budget_sum_param:
        check_budget_sum_vars = budget_sum_vars_all
        budget_args = {"write_sum_vars": budget_sum_param}
    elif not budget_sum_param:
        check_budget_sum_vars = []
        budget_args = {"write_sum_vars": budget_sum_param}
    else:
        raise ValueError("upexpected value")

    # dont need any PRMS inputs for the model specified, so this is sufficient
    input_dir = domain["prms_output_dir"]

    # TODO: Eliminate potet and other variables from being used
    model = Model(
        *model_procs,
        control=control,
        input_dir=input_dir,
        budget_type=budget_type,
    )

    check_dict = {proc: {} for proc in check_vars.keys()}

    model.initialize_netcdf(tmp_dir, budget_args=budget_args)

    for tt in range(n_time_steps):
        model.advance()
        model.calculate()
        model.output()

        for pp, pp_vars in check_vars.items():
            for vv in pp_vars:
                if tt == 0:
                    # use the output data to figure out the shape
                    check_dict[pp][vv] = np.zeros(
                        (n_time_steps, model.processes[pp][vv].shape[0])
                    )

                check_dict[pp][vv][tt, :] = model.processes[pp][vv]

            for bb in check_budget_sum_vars:
                if tt == 0:
                    # use the output data to figure out the shape
                    check_dict[pp][bb] = np.zeros(
                        (
                            n_time_steps,
                            model.processes[pp].budget[f"_{bb}"].shape[0],
                        )
                    )

                check_dict[pp][bb][tt, :] = model.processes[pp].budget[
                    f"_{bb}"
                ]

    model.finalize()

    # read the data back in
    for pp, pp_vars in check_vars.items():
        for vv in pp_vars:
            nc_data = xr.open_dataset(tmp_dir / f"{vv}.nc")[vv]
            assert np.allclose(check_dict[pp][vv], nc_data)

        for bb in check_budget_sum_vars:
            nc_data = xr.open_dataset(tmp_dir / f"{pp}_budget.nc")[bb]
            assert np.allclose(check_dict[pp][bb], nc_data)

        if budget_sum_param == "some":
            nc_data = xr.open_dataset(tmp_dir / f"{pp}_budget.nc")
            for nn in set(budget_sum_vars_all).difference(
                set(check_budget_sum_vars)
            ):
                assert nn not in nc_data.variables
        elif not budget_sum_param:
            assert not (tmp_dir / f"{pp}_budget.nc").exists()

    return
