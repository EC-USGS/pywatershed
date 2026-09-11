"""separate_domain_params_dis_to_ncdf with a cascade control.

Checks that the per-process parameter files written for the cascade
process classes are complete: every declared parameter is on file (or is
a discretization parameter), and every declared dimension is on file even
when no parameter uses it (nsegment). Runs only on cascade controls;
CI ignores this file in the broad steps.
"""

import pathlib as pl

import pytest

import pywatershed
from pywatershed.base.data_model import open_datasetdict
from pywatershed.parameters import PrmsParameters
from pywatershed.utils import separate_domain_params_dis_to_ncdf

cascade_processes = [
    pywatershed.PRMSRunoffCascadesNoDprst,
    pywatershed.PRMSSoilzoneCascadesNoDprst,
]


@pytest.fixture(scope="function")
def control(simulation):
    ctl = pywatershed.Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    if not ctl.options.get("cascade_flag", 0):
        pytest.skip("cascade_flag absent or 0")
    del ctl.options["netcdf_output_dir"]
    return ctl


@pytest.fixture(scope="function")
def params(simulation, control):
    param_file = simulation["dir"] / control.options["parameter_file"]
    return PrmsParameters.load(param_file)


def test_param_sep_cascades(simulation, control, params, tmp_path):
    prms_param_file = simulation["dir"] / control.options["parameter_file"]
    proc_nc_files = separate_domain_params_dis_to_ncdf(
        prms_param_file,
        simulation["name"],
        pl.Path(tmp_path),
        process_list=cascade_processes,
        control=control,
    )

    prms_names = set(params.variables.keys())
    for proc in cascade_processes:
        file_params = open_datasetdict(proc_nc_files[proc])
        file_names = set(file_params.variables.keys())
        # both loaders must see the same coordinates, including nhm_seg
        # which no parameter uses (only the global coordinates attribute
        # records it)
        file_params_nc4 = open_datasetdict(proc_nc_files[proc], use_xr=False)
        assert set(file_params_nc4.coords.keys()) == set(
            file_params.coords.keys()
        )
        assert "nhm_seg" in file_params.coords.keys()
        # every declared parameter is on file or comes from the PRMS file
        # via a discretization
        missing = set(proc.get_parameters()) - file_names - prms_names
        assert not missing, f"{proc.__name__} file lacks {missing}"
        # the cascade parameters are not in the PRMS file; they must be here
        assert "hru_route_order" in file_names
        assert "hru_down" in file_names
        # every declared dimension is on file, used by a parameter or not
        missing_dims = set(proc.get_dimensions()) - set(file_params.dims)
        assert not missing_dims, f"{proc.__name__} file lacks {missing_dims}"
