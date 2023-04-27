import pathlib as pl

import numpy as np
import pytest
import xarray as xr

from pywatershed.parameters import PrmsParameters
from pywatershed.utils.cbh_utils import cbh_files_to_df, cbh_files_to_netcdf
from utils import assert_or_print

var_cases = ["prcp", "rhavg", "tmax", "tmin"]


@pytest.fixture
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(params=["no_params", "params"])
def params_and_none(request, domain):
    if request.param == "params":
        return PrmsParameters.load(domain["param_file"])
    else:
        return None


answer_key = {
    "drb_2yr": {
        "files_to_df": {
            "prcp": 0.12495362248866715,
        },
        "df_concat": 40.7829059976932,
        "np_dict": {
            "prcp": 0.12495362248866718,
            "rhavg": 62.738591579267386,
            "tmax": 60.26121597238987,
            "tmin": 40.00686281662687,
            "time": 315532800.0,
        },
    },
    "hru_1": {
        "files_to_df": {
            "prcp": 0.13502103505843072,
        },
        "df_concat": 34.968545798553144,
        "np_dict": {
            "prcp": 0.13502103505843072,
            "tmax": 61.986048747913195,
            "tmin": 42.78456761268781,
            "time": 930873600.0,
        },
    },
    "ucb_2yr": {
        "files_to_df": {
            "prcp": 0.046485653521159784,
        },
        "df_concat": 34.05471072235577,
        "np_dict": {
            "prcp": 0.046485653521159805,
            "rhavg": 50.87569199962628,
            "tmax": 56.74147989347377,
            "tmin": 28.555185342801856,
            "time": 315532800.0,
        },
    },
}


@pytest.mark.parametrize("var", [var_cases[0]])
def test_cbh_files_to_df(domain, var, params):
    the_file = domain["cbh_inputs"][var]
    assert pl.Path(the_file).exists(), the_file
    df = cbh_files_to_df(the_file, params)
    results = {var: df.mean().mean()}
    answers = answer_key[domain["domain_name"]]["files_to_df"]
    assert_or_print(
        results, answers, "files_to_df", print_ans=domain["print_ans"]
    )
    return


def test_cbh_files_to_netcdf(domain, params, tmp_path):
    nc_file = tmp_path / "cbh_files_to_netcdf.nc"
    input_files_dict = domain["cbh_inputs"]
    _ = cbh_files_to_netcdf(input_files_dict, params, nc_file)

    answers = answer_key[domain["domain_name"]]["np_dict"]
    results_ds = xr.open_dataset(nc_file)
    results = {}
    for var in answers.keys():
        if var == "time":
            results[var] = (
                results_ds[var]
                .values.astype("datetime64[s]")
                .astype(np.int32)
                .mean()
            )
        else:
            results[var] = results_ds[var].mean().mean().values.tolist()

    assert_or_print(
        results,
        answers,
        "files_to_np_dict",
        print_ans=domain["print_ans"],
        close=True,
    )

    return
