import pathlib as pl
import warnings

import numpy as np
import pytest
import xarray as xr

from pynhm.preprocess.cbh import CBH
from pynhm.preprocess.cbh_utils import cbh_files_to_df
from pynhm.utils import PrmsParameters, timer
from pynhm.utils.prms5util import load_prms_statscsv
from utils import assert_or_print

var_cases = ["prcp", "rhavg", "tmax", "tmin"]


# Reduce IO/loading for the domain parameters
@pytest.fixture
def params(domain):
    return PrmsParameters.load(domain["param_file"])


# It appears difficult/ugly to reuse the above
@pytest.fixture(params=["no_params", "params"])
def params_and_none(request, domain):
    if request.param == "params":
        return PrmsParameters.load(domain["param_file"])
    else:
        return None


@timer
def cbh_files_to_df_timed(*args):
    df = cbh_files_to_df(*args)
    return df


# -------------------------------------------------------
# Test Invalid Instantiation
case_dict_init = {
    "invalid": ["a", "b"],
}


@pytest.mark.parametrize("case", list(case_dict_init.keys()))
def test_cbh_init_files_invalid(case):
    try:
        _ = CBH(case_dict_init[case])
        assert False  # This case is supposed to fail
    except ValueError:
        if case in ["invalid"]:
            assert True
        else:
            assert False
    return


# Valid instantiation happens in other tests.


# -------------------------------------------------------
# Test the file to df parser on a *single file*
# One case is likely sufficient
@pytest.mark.parametrize("var", [var_cases[0]])
@timer
def test_cbh_files_to_df(domain, var, params):
    the_file = domain["cbh_inputs"][var]
    assert pl.Path(the_file).exists(), the_file
    df = cbh_files_to_df(the_file, params)
    results = {var: df.mean().mean()}
    answers = domain["test_ans"]["preprocess_cbh"]["files_to_df"]
    assert_or_print(
        results, answers, "files_to_df", print_ans=domain["print_ans"]
    )
    return


# -------------------------------------------------------
# Test concatenation of individual files
@timer
def test_cbh_df_concat(domain):
    df = cbh_files_to_df(domain["input_files_dict"])
    results = {"df_concat": df.mean().mean()}
    answers = {"df_concat": domain["test_ans"]["preprocess_cbh"]["df_concat"]}
    assert_or_print(results, answers, print_ans=domain["print_ans"])
    return


# -------------------------------------------------------
# Test conversion to numpy dict
def test_cbh_np_dict(domain, params):
    # Could hide the above in the domain on read
    cbh = CBH(domain["input_files_dict"], parameters=params)
    dont_average = ["datetime", "hru_ind", "nhm_id"]
    results = {
        key: val.mean()
        for key, val in cbh.state.items()
        if key not in dont_average
    }
    results["datetime"] = cbh.state["datetime"].astype(int).mean()
    answers = domain["test_ans"]["preprocess_cbh"]["np_dict"]
    assert_or_print(results, answers, "np_dict", print_ans=domain["print_ans"])
    return


# -------------------------------------------------------
# Test forcing adjustment
def test_cbh_adj(domain, params_and_none):
    cbh = CBH(
        domain["input_files_dict"], parameters=params_and_none, adjust=True
    )
    if params_and_none is not None:
        results = {
            key: val.mean() for key, val in cbh.state.items() if "_adj" in key
        }
        answers = domain["test_ans"]["preprocess_cbh"]["adj"]["params"]
        assert_or_print(
            results, answers, "adj:params", print_ans=domain["print_ans"]
        )
    else:
        results = {
            key: val.mean()
            for key, val in cbh.state.items()
            if key in var_cases
        }
        answers = domain["test_ans"]["preprocess_cbh"]["adj"]["no_params"]
        assert_or_print(
            results, answers, "adj:params", print_ans=domain["print_ans"]
        )

    return


# -------------------------------------------------------
# Test forcing adj, compare to prms output
def test_cbh_adj_prms_output(domain, params):
    cbh = CBH(domain["input_files_dict"], params, adjust=True)
    for var, var_file in domain["prms_outputs"].items():
        if var in ["soltab"]:
            continue
        if var in ["swrad", "potet"]:
            msg = (
                f"Skipping {var} as it is not currently preprocessed, "
                f"this skip should be removed when it is"
            )
            warnings.warn(msg)
            continue
        prms_output = load_prms_statscsv(var_file)
        p_dates = prms_output.index.values
        p_array = prms_output.to_numpy()
        wh_dates = np.where(np.isin(cbh.state["datetime"], p_dates))
        result = np.isclose(
            p_array,
            cbh.state[var][wh_dates, :],
            rtol=1e-121,
            atol=1e-04,  # Only the atol matters here, if atol < 1e-4 fails
        )
        assert result.all()

    return


# -------------------------------------------------------
# # Test netcdf write
# @pytest.mark.parametrize("use_params", [False, True])
def test_cbh_to_netcdf(domain, tmp_path, params_and_none):
    cbh = CBH(domain["input_files_dict"], params_and_none, adjust=True)
    tmp_file = tmp_path / "test_cbh_to_netcdf.nc"
    print(tmp_file)
    global_atts = {"domain_name": domain["domain_name"]}
    _ = cbh.to_netcdf(tmp_file, global_atts=global_atts)
    ds = xr.open_dataset(tmp_file)
    # compare the calculated means in memory and after reading the data from disk
    for vv in ds.variables:
        if vv in ["hru_ind", "nhm_id"]:
            assert np.isclose(
                ds[vv].astype(float).values.mean(),
                cbh.state[vv].astype(float).mean(),
            )
        elif vv == "datetime":
            assert np.isclose(
                ds[vv].astype(float).values.mean(),
                cbh.state["datetime"]
                .astype("datetime64[ns]")
                .astype(float)
                .mean(),
            )
        else:
            assert np.isclose(
                ds[vv].mean().values.tolist(),
                cbh.state[vv].mean(),
            )

    ds.close()
    return
