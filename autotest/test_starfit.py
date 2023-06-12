import pathlib as pl
import pickle

import numpy as np
import pytest
import xarray as xr

from pywatershed.base.control import Control
from pywatershed.hydrology.Starfit import Starfit
from pywatershed.parameters import StarfitParameters


# havent pared down the data yet to add it to the repo
@pytest.mark.xfail
def test_regress(tmp_path):
    # Regression against independenly run outputs in pickle file

    tmp_path = pl.Path(tmp_path)
    print(tmp_path)
    data_dir = pl.Path(
        "/Users/jmccreight/usgs/pywatershed/pywatershed/hydrology/starfit_minimal"
    )

    # TODO: make this work with the original source files
    param_files = {
        "grand_ids": data_dir / "domain_grand_ids.nc",
        "resops_domain": data_dir / "ResOpsUS.nc",
        "istarf_conus": data_dir / "starfit/ISTARF-CONUS.nc",
        "grand_dams": data_dir / "grand_dams.nc",
    }

    # passing the parameter names is optional
    starfit_parameters = Starfit.get_parameters()
    params = StarfitParameters.from_netcdf(
        **param_files,
        param_names=starfit_parameters,
    )

    control = Control(
        params.variables["start_time"].min(),
        np.datetime64("2019-09-30 00:00:00"),
        np.timedelta64(24, "h"),
        params=params,
    )

    # # load csv files into dataframes
    # output_dir = domain["prms_output_dir"]
    input_variables = {}
    for key in Starfit.get_inputs():
        nc_path = data_dir / f"{key}.nc"
        input_variables[key] = nc_path

    sf = Starfit(
        control,
        **input_variables,
        budget_type="warn",
    )

    sf.initialize_netcdf(tmp_path)

    for istep in range(control.n_times):
        control.advance()

        sf.advance()

        sf.calculate(float(istep))

        sf.output()

    sf.finalize()

    del sf

    results = {}
    vars_compare = ["lake_storage", "lake_release", "lake_spill"]
    for var in vars_compare:
        results[var] = xr.open_dataset(tmp_path / f"{var}.nc")[var]

    # Regression against independenly run outputs in pickle file
    with open(data_dir / "starfit_outputs.pickle", "rb") as handle:
        ans_dict = pickle.load(handle)

    answers = {}
    answers["lake_storage"] = ans_dict["Ssim_list"]
    answers["lake_release"] = ans_dict["Rsim_list"]
    answers["lake_spill"] = ans_dict["SPILLsim_list"]
    ans_grand_ids = ans_dict["grand_ids"]

    for ii, gi in enumerate(results["lake_storage"].grand_id):
        assert ans_grand_ids[ii] == gi  # extra dummy check
        for vv in vars_compare:
            res = results[vv].loc[dict(grand_id=gi)].load()
            ans = answers[vv][ii].to_xarray().rename(index="time")
            diff = res - ans
            # diff is only on common times
            # make sure that no diffs are nan
            assert not np.isnan(diff).any()
            # make sure no times simulated that were not requested
            assert len(np.where(~np.isnan(res))[0]) == len(ans)
            tol = 1.0e-6
            wh_abs_diff = np.where((abs(diff) > tol) & ((abs(diff) / ans) > tol))
            assert len(wh_abs_diff[0]) == 0

    return


@pytest.mark.xfail
@pytest.mark.parametrize("domain_tag", ["drb"])
def test_param_subset_write(tmp_path, domain_tag):
    tmp_path = pl.Path(tmp_path)
    print(tmp_path)
    data_dir = pl.Path(
        "/Users/jmccreight/usgs/pywatershed/pywatershed/hydrology/starfit_minimal"
    )

    # TODO: make this work with the original source files
    param_files = {
        "grand_ids": data_dir / f"grand_seg_ids_{domain_tag}.nc",
        "resops_domain": data_dir / "ResOpsUS.nc",
        "istarf_conus": data_dir / "starfit/ISTARF-CONUS.nc",
        "grand_dams": data_dir / "grand_dams.nc",
    }

    # passing the parameter names is optional
    starfit_parameters = Starfit.get_parameters()
    params = StarfitParameters.from_netcdf(
        **param_files,
        param_names=starfit_parameters,
    )

    params.to_netcdf(tmp_path / f"starfit_params_{domain_tag}.nc", use_xr=True)
