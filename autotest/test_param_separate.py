import pathlib as pl
import shutil

import numpy as np
import pytest

import pynhm
from pynhm.base.control import Control
from pynhm.base.timeseries import TimeseriesArray
from pynhm.utils import separate_domain_params_to_ncdf
from pynhm.parameters import PrmsParameters

n_time_steps = 10
budget_type = None
pynhm_processes = [
    # [pynhm.PRMSSolarGeometry],
    # [pynhm.PRMSAtmosphere],
    # [pynhm.PRMSCanopy],
    # [pynhm.PRMSSnow],
    # [pynhm.PRMSRunoff],
    # [pynhm.PRMSSoilzone],
    [pynhm.PRMSGroundwater],
    [pynhm.PRMSChannel],
    [pynhm.PRMSGroundwater, pynhm.PRMSChannel],
]
pynhm_process_id = [
    "_".join([proc.__name__ for proc in procs]) for procs in pynhm_processes
]


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


@pytest.mark.parametrize(
    "processes",
    pynhm_processes,
    ids=pynhm_process_id,
)
def test_param_sep(domain, control, processes, tmp_path):
    tmp_path = pl.Path(tmp_path)

    domain_name = domain["domain_name"]
    prms_param_file = domain["param_file"]

    out_dir = tmp_path / f"{domain_name}"
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    param_nc_files = separate_domain_params_to_ncdf(
        domain_name, prms_param_file, out_dir, process_list=processes
    )
    assert len(param_nc_files) == len(processes)
    params_sep = PrmsParameters.from_nc_files(param_nc_files)

    for proc_class in params_sep.parameters.keys():
        params_indiv = params_sep.parameters[proc_class]
        dims_indiv = params_sep.parameter_dimensions[proc_class]
        # same parameter names in both
        assert set(params_indiv.keys()) == set(proc_class.get_parameters())
        for param in proc_class.get_parameters():
            check_vals = (
                params_indiv[param] == control.params.parameters[param]
            )
            check_dims = (
                dims_indiv[param] == control.params.parameter_dimensions[param]
            )
            for check_it in [check_vals, check_dims]:
                if isinstance(check_it, np.ndarray):
                    assert check_it.all()
                else:
                    assert check_it

    # advanced model run-based checking code
    # probably NOT necessary but can prototype how to pass
    # parameter_dicts from netcdf files
    prms_output_dir = domain["prms_output_dir"]

    # setup input_dir with symlinked prms inputs and outputs
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    for ff in prms_output_dir.resolve().glob("*.nc"):
        shutil.copy(ff, input_dir / ff.name)
    for ff in prms_output_dir.parent.resolve().glob("*.nc"):
        shutil.copy(ff, input_dir / ff.name)

    model_prms_params = pynhm.Model(
        *processes,
        control=control,
        input_dir=input_dir,
        budget_type=budget_type,
    )
    model_prms_params.run(n_time_steps=n_time_steps)

    control_params_sep = Control.load(
        domain["control_file"], params=params_sep
    )
    model_sep_params = pynhm.Model(
        *processes,
        control=control_params_sep,
        input_dir=input_dir,
        budget_type=budget_type,
    )
    model_sep_params.run(n_time_steps=n_time_steps)

    for proc_key in model_sep_params.processes.keys():
        for var in model_sep_params.processes[proc_key].variables:
            prms_param_vals = model_prms_params.processes[proc_key][var]
            sep_param_vals = model_sep_params.processes[proc_key][var]
            if isinstance(prms_param_vals, TimeseriesArray):
                prms_param_vals = prms_param_vals.current
                sep_param_vals = sep_param_vals.current

            assert (prms_param_vals == sep_param_vals).all()
