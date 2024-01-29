import pathlib as pl

import numpy as np
import pytest
import xarray as xr

from pywatershed import Control
from pywatershed.parameters import PrmsParameters
from pywatershed.utils.cbh_utils import cbh_file_to_netcdf, cbh_files_to_df
from utils import assert_or_print

var_cases = ["prcp", "rhavg", "tmax", "tmin"]


@pytest.fixture(scope="function")
def control(simulation):
    ctl = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    del ctl.options["netcdf_output_dir"]
    return ctl


@pytest.fixture
def params(simulation, control):
    param_file = simulation["dir"] / control.options["parameter_file"]
    return PrmsParameters.load(param_file)


@pytest.fixture(params=["no_params", "params"])
def params_and_none(request, simulation, control):
    if request.param == "params":
        param_file = simulation["dir"] / control.options["parameter_file"]
        return PrmsParameters.load(param_file)
    else:
        return None


answer_key = {
    "drb_2yr:nhm": {
        "files_to_df": {
            "prcp": 0.12495362248866715,
        },
        "df_concat": 40.7829059976932,
        "file_to_netcdf": {
            "prcp": 0.12495362248866718,
            "rhavg": 62.738591579267386,
            "tmax": 60.26121597238987,
            "tmin": 40.00686281662687,
            "time": 315532800.0,
        },
    },
    "hru_1:nhm": {
        "files_to_df": {
            "prcp": 0.13502103505843072,
        },
        "df_concat": 34.968545798553144,
        "file_to_netcdf": {
            "prcp": 0.13502103505843072,
            "tmax": 61.986048747913195,
            "tmin": 42.78456761268781,
            "time": 930873600.0,
        },
    },
    "ucb_2yr:nhm": {
        "files_to_df": {
            "prcp": 0.046485653521159784,
        },
        "df_concat": 34.05471072235577,
        "file_to_netcdf": {
            "prcp": 0.046485653521159805,
            "rhavg": 50.87569199962628,
            "tmax": 56.74147989347377,
            "tmin": 28.555185342801856,
            "time": 315532800.0,
        },
    },
}


@pytest.mark.parametrize("var", [var_cases[0]])
def test_cbh_files_to_df(simulation, var, params):
    the_file = simulation["dir"] / f"{var}.cbh"
    assert pl.Path(the_file).exists(), the_file
    df = cbh_files_to_df(the_file, params)
    results = {var: df.mean().mean()}
    answers = answer_key[simulation["name"]]["files_to_df"]
    assert_or_print(
        results, answers, "files_to_df", print_ans=simulation["print_ans"]
    )
    return


def test_cbh_file_to_netcdf(simulation, params, tmp_path):
    input_files_dict = {
        ff.with_suffix("").name: ff for ff in simulation["dir"].glob("*.cbh")
    }

    results = {}
    for var, cbh_file in input_files_dict.items():
        nc_file = tmp_path / f"{var}.nc"
        _ = cbh_file_to_netcdf(cbh_file, params, nc_file)
        results[var] = xr.open_dataset(nc_file)[var]

    results["time"] = results[var].time

    answers = answer_key[simulation["name"]]["file_to_netcdf"]
    for var in answers.keys():
        if var == "time":
            results[var] = (
                results[var]
                .values.astype("datetime64[s]")
                .astype(np.int32)
                .mean()
            )
        else:
            results[var] = results[var].mean().mean().values.tolist()

    assert_or_print(
        results,
        answers,
        "file_to_netcdf",
        print_ans=simulation["print_ans"],
        close=True,
    )

    return
