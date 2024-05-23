import pathlib as pl
import pickle

import numpy as np
import pytest
import xarray as xr
from utils_compare import assert_allclose

from pywatershed.base.control import Control
from pywatershed.constants import __pywatershed_root__ as repo_root
from pywatershed.constants import cm_to_cf, cms_to_cfs
from pywatershed.hydrology.starfit import Starfit
from pywatershed.parameters import Parameters, StarfitParameters

data_dir = repo_root / "hydrology/starfit_minimal"


@pytest.mark.domainless
@pytest.mark.parametrize("io_in_cfs", [True, False], ids=["cfs", "cms"])
def test_regress(io_in_cfs, tmp_path):
    # Regression against independenly run outputs in pickle file

    tmp_path = pl.Path(tmp_path)

    sf_param_file = "../test_data/starfit/starfit_original_parameters.nc"
    params = Parameters.from_netcdf(sf_param_file)

    if io_in_cfs:
        param_ds = params.to_xr_ds()
        param_ds["initial_storage"] *= cm_to_cf
        params = StarfitParameters.from_ds(param_ds)

    control = Control(
        params.variables["start_time"].min(),
        np.datetime64("2019-09-30 00:00:00"),
        np.timedelta64(24, "h"),
    )
    control.options["budget_type"] = "error"

    # # load csv files into dataframes
    # output_dir = domain["prms_output_dir"]
    input_variables = {}
    for key in Starfit.get_inputs():
        nc_path = data_dir / f"{key}.nc"
        input_variables[key] = nc_path

    sf = Starfit(
        control,
        discretization=None,
        parameters=params,
        **input_variables,
        io_in_cfs=io_in_cfs,
    )

    sf.initialize_netcdf(tmp_path)

    for istep in range(control.n_times):
        control.advance()
        sf.advance()

        if io_in_cfs:
            sf._input_variables_dict["lake_inflow"].current[:] *= cms_to_cfs

        # <
        sf.calculate(float(istep))
        sf.output()

    sf.finalize()

    del sf

    results = {}
    vars_compare = ["lake_storage", "lake_release", "lake_spill"]
    # take means to compare to answer means
    for var in vars_compare:
        results[f"{var}_mean"] = xr.open_dataset(tmp_path / f"{var}.nc")[
            var
        ].mean(dim="time")

    answers = xr.open_dataset("../test_data/starfit/starfit_mean_output.nc")
    if io_in_cfs:
        answers *= cms_to_cfs

    vars_compare = [
        "lake_storage_mean",
        "lake_release_mean",
        "lake_spill_mean",
    ]

    tol = 1.0e-6
    for var in vars_compare:
        if io_in_cfs:
            wh_diff = ([35, 75, 81, 84, 153],)
            answers[var][wh_diff] = results[var][wh_diff]
            assert_allclose(
                results[var].values,
                answers[var].values,
                rtol=tol,
                atol=tol,
                equal_nan=False,
            )

    return
