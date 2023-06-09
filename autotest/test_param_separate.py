import pathlib as pl
import shutil

import numpy as np
import pytest

import pywatershed
from pywatershed.base.control import Control
from pywatershed.base.timeseries import TimeseriesArray
from pywatershed.parameters import PrmsParameters
from pywatershed.utils import separate_domain_params_to_ncdf

n_time_steps = 10
budget_type = None
pyws_processes = [
    [pywatershed.PRMSSolarGeometry],
    [pywatershed.PRMSAtmosphere],
    [pywatershed.PRMSCanopy],
    [pywatershed.PRMSSnow],
    [pywatershed.PRMSRunoff],
    [pywatershed.PRMSSoilzone],
    [pywatershed.PRMSGroundwater],
    [pywatershed.PRMSChannel],
    [pywatershed.PRMSGroundwater, pywatershed.PRMSChannel],
]
pyws_process_id = [
    "_".join([proc.__name__ for proc in procs]) for procs in pyws_processes
]


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


@pytest.mark.parametrize(
    "processes",
    pyws_processes,
    ids=pyws_process_id,
)
def test_param_sep(domain, control, processes, tmp_path):
    tmp_path = pl.Path(tmp_path)

    domain_name = domain["domain_name"]
    prms_param_file = domain["param_file"]

    out_dir = tmp_path / f"{domain_name}"
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    # separate to netcdf
    param_nc_files = separate_domain_params_to_ncdf(
        domain_name, prms_param_file, out_dir, process_list=processes
    )
    assert len(param_nc_files) == len(processes)
    # read back in
    params_sep = PrmsParameters.from_nc_files(param_nc_files)

    # check roundtrip
    for proc_class, proc_val in params_sep.items():
        assert set(proc_val.parameters.keys()) == set(proc_class.get_parameters())
        for param in proc_class.get_parameters():
            np.testing.assert_equal(
                proc_val.parameters[param], control.params.parameters[param]
            )
            np.testing.assert_equal(
                proc_val.metadata[param]["dims"],
                control.params.metadata[param]["dims"],
            )

    return
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

    model_prms_params = pywatershed.Model(
        *processes,
        control=control,
        input_dir=input_dir,
        budget_type=budget_type,
    )
    model_prms_params.run(n_time_steps=n_time_steps)

    # pulling out the first parameter value is a simple trick that wont work
    # if there are multiple processes in play.
    control_params_sep = Control.load(
        domain["control_file"], params=list(params_sep.values)[0]
    )
    model_sep_params = pywatershed.Model(
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
