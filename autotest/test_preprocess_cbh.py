import pathlib as pl

import numpy as np
import pytest
import xarray as xr

from pynhm.preprocess.cbh import CBH
from pynhm.preprocess.cbh_utils import cbh_files_to_df
from pynhm.utils import PrmsParameters, timer
from pynhm.utils.prms5util import load_prms_statscsv
from utils import assert_or_print


var_cases = ["prcp", "rhavg", "tmax", "tmin"]


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
def test_cbh_files_to_df(domain, var):
    the_file = domain["cbh_inputs"][var]
    assert pl.Path(the_file).exists(), the_file
    df = cbh_files_to_df_timed(the_file)
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
    answers = domain["test_ans"]["preprocess_cbh"]
    assert_or_print(results, answers, print_ans=domain["print_ans"])
    return


# -------------------------------------------------------
# Test conversion to numpy dict
def test_cbh_np_dict(domain):
    # Could hide the above in the domain on read
    cbh = CBH(domain["input_files_dict"])
    results = {
        key: val.mean() for key, val in cbh.state.items() if key != "datetime"
    }
    results["datetime"] = cbh.state["datetime"].astype(int).mean()
    answers = domain["test_ans"]["preprocess_cbh"]["np_dict"]
    assert_or_print(results, answers, "np_dict", print_ans=domain["print_ans"])
    return


# -------------------------------------------------------
# Test forcing adjustment
def test_cbh_adj(domain):
    params = PrmsParameters(domain["param_file"])
    cbh = CBH(domain["input_files_dict"])  # JLM also test passing on init?
    _ = cbh.adjust(params)
    results = {
        key: val.mean() for key, val in cbh.state.items() if "_adj" in key
    }
    answers = domain["test_ans"]["preprocess_cbh"]["adj"]
    assert_or_print(results, answers, "adj", print_ans=domain["print_ans"])
    return


# -------------------------------------------------------
# Test forcing adj, compare to prms output
case_list_adj_prms_output = ["drb_2yr"]


def test_cbh_adj_prms_output(domain):
    params = PrmsParameters(domain["param_file"])
    cbh = CBH(domain["input_files_dict"])

    @timer
    def adj_timed(params):
        return cbh.adjust(params)

    _ = adj_timed(params)

    for var, var_file in domain["prms_outputs"].items():
        prms_output = load_prms_statscsv(var_file)
        p_dates = prms_output.index.values
        p_array = prms_output.to_numpy()
        wh_dates = np.where(np.isin(cbh.state["datetime"], p_dates))
        # print(var)
        result = np.isclose(
            p_array,
            cbh.state[var][wh_dates, :],
            rtol=1e-121,
            atol=1e-04,  # Only the atol matters here, if atol < 1e-4 fails
        )
        # print(result.sum()/np.prod(result.shape))
        assert result.all()

    return


# -------------------------------------------------------
# # Test netcdf write
def test_cbh_to_netcdf(domain, tmp_path):
    cbh = CBH(domain["input_files_dict"])
    tmp_file = tmp_path / "test_cbh_to_netcdf.nc"
    _ = cbh.to_netcdf(tmp_file)
    ds = xr.open_dataset(tmp_file)
    # compare the calculated means in memory and after reading the data from disk
    for vv in ds.variables:
        if vv == "hruid":
            # for now this is just an index
            assert ds[vv].mean().values.tolist() == ((len(ds[vv]) - 1) / 2)
        elif vv == "datetime":
            assert np.isclose(
                ds[vv].mean().values.tolist(),
                cbh.state["datetime"]
                .astype("datetime64[ns]")
                .astype(int)
                .mean(),
            )
        else:
            assert np.isclose(
                ds[vv].mean().values.tolist(),
                cbh.state[vv].mean(),
            )

    ds.close()
    return
