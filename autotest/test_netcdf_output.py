import pathlib as pl
import shutil

import numpy as np
import pytest
import xarray as xr

import pywatershed
from pywatershed.base.control import Control
from pywatershed.base.model import Model
from pywatershed.parameters import PrmsParameters

# test for a few timesteps a model with both unit/cell and global balance
# budgets


n_time_steps = 10


# This probably dosent need varied over domain
@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain):
    control = Control.load(domain["control_file"])
    control.edit_n_time_steps(n_time_steps)
    return control


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
def test_process_budgets(domain, control, params, tmp_path, budget_sum_param):
    tmp_dir = pl.Path(tmp_path)
    # print(tmp_dir)
    model_procs = [pywatershed.PRMSCanopy, pywatershed.PRMSChannel]

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
        model_procs,
        control=control,
        discretization_dict=None,
        parameters=params,
        input_dir=input_dir,
        budget_type=budget_type,
    )

    check_dict = {proc: {} for proc in check_vars.keys()}

    # test outputting specific vars by only using check_vars
    output_vars = [
        item for sublist in list(check_vars.values()) for item in sublist
    ]
    output_vars = None

    model.initialize_netcdf(
        tmp_dir,
        budget_args=budget_args,
        output_vars=output_vars,
    )

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


@pytest.mark.parametrize(
    "separate",
    [False, True],
    ids=["grp_by_process", "separate"],
)
def test_separate_together(domain, control, params, tmp_path, separate):
    tmp_dir = pl.Path(tmp_path)

    model_procs = [
        pywatershed.PRMSSolarGeometry,
        pywatershed.PRMSAtmosphere,
        pywatershed.PRMSCanopy,
        pywatershed.PRMSChannel,
    ]

    # setup input_dir with symlinked prms inputs and outputs
    test_output_dir = tmp_dir / "test_results"
    domain_output_dir = domain["prms_output_dir"]
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    # Could limit this to just the variables in model_procs
    for ff in domain_output_dir.resolve().glob("*.nc"):
        shutil.copy(ff, input_dir / ff.name)
    for ff in domain_output_dir.parent.resolve().glob("*.nc"):
        shutil.copy(ff, input_dir / ff.name)

    model = Model(
        model_procs,
        control=control,
        discretization_dict=None,
        parameters=params,
        input_dir=input_dir,
        budget_type=budget_type,
    )

    model.initialize_netcdf(
        output_dir=test_output_dir,
        separate_files=separate,
    )

    for tt in range(n_time_steps):
        model.advance()
        model.calculate()
        model.output()

    model.finalize()

    if separate:
        for proc_key, proc in model.processes.items():
            for vv in proc.variables:
                nc_file = test_output_dir / f"{vv}.nc"
                assert nc_file.exists()
                ds = xr.open_dataset(nc_file, decode_timedelta=False)
                if isinstance(
                    proc[vv], pywatershed.base.timeseries.TimeseriesArray
                ):
                    assert (ds[vv].values == proc[vv].data).all()
                else:
                    assert (ds[vv][-1, :] == proc[vv]).all()

                del ds

    else:
        for proc_key, proc in model.processes.items():
            # non-budget
            nc_file = test_output_dir / f"{proc_key}.nc"
            assert nc_file.exists()
            ds = xr.open_dataset(nc_file, decode_timedelta=False)
            proc_vars = set(proc.get_variables())
            nc_vars = set(ds.data_vars)
            assert proc_vars == nc_vars
            for vv in proc.variables:
                if isinstance(
                    proc[vv], pywatershed.base.timeseries.TimeseriesArray
                ):
                    assert (ds[vv].values == proc[vv].data).all()

                else:
                    assert (ds[vv][-1, :] == proc[vv]).all()

            del ds

            # budget
            # no budgets for solar or atmosphere
            if proc_key in ["PRMSSolarGeometry", "PRMSAtmosphere"]:
                continue
            nc_file = test_output_dir / f"{proc_key}_budget.nc"
            ds = xr.open_dataset(nc_file)
            for ss in budget_sum_vars_all:
                assert (proc.budget[ss] == ds[ss][-1, :]).all()

            del ds
    return
