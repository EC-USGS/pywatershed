import numpy as np
import pathlib as pl
from pynhm.preprocess.cbh import CBH
from pynhm.utils import timer
from pynhm.utils.prms5util import load_prms_statscsv
from pynhm.utils import PrmsParameters
from pynhm.preprocess.cbh_utils import cbh_files_to_df
import pytest
import xarray as xr

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


# Valid instantiation
@pytest.mark.parametrize("case", var_cases)
def test_cbh_init_files_valid(domain, case):
    _ = CBH(domain["dir"] / domain["cbh_inputs"][case])
    return


# -------------------------------------------------------
# Test the file to df parser on a *single file*
@pytest.mark.parametrize("var", var_cases)
def test_cbh_files_to_df(domain, var):
    the_file = domain["dir"] / domain["cbh_inputs"][var]
    assert pl.Path(the_file).exists(), the_file  # rm pathlib?
    df = cbh_files_to_df_timed(the_file)
    # print(repr(f"'{var}': ('{repr(df.iloc[-1, :])}'),"))
    assert (
        repr(df.iloc[-1, :])
        == domain["test_ans"]["preprocess_cbh"]["files_to_df"][var]
    )
    return


# -------------------------------------------------------
# Test concatenation of individual files
def test_cbh_concat(domain):
    df = cbh_files_to_df(domain["input_files_dict"])
    # print(repr(f"'{case}': ('{repr(df.iloc[-1, :])}'),"))
    assert (
        repr(df.iloc[-1, :])
        == domain["test_ans"]["preprocess_cbh"]["df_concat"]
    )
    return


# -------------------------------------------------------
# Test conversion to numpy dict
def test_cbh_np_dict(domain):
    # Could hide the above in the domain on read
    cbh = CBH(domain["input_files_dict"])
    with np.printoptions(threshold=8):
        for key, val in cbh.state.items():
            # print(repr(f'"{key}": ("{repr(val)}"),'))
            assert (
                repr(val)
                == domain["test_ans"]["preprocess_cbh"]["np_dict"][key]
            )
    return


# -------------------------------------------------------
# Test forcing adjustment
def test_cbh_adj(domain):
    params = PrmsParameters(domain["dir"] / domain["param_file"])
    cbh = CBH(domain["input_files_dict"])  # JLM also test passing on init?
    _ = cbh.adjust(params)
    with np.printoptions(threshold=2):
        for key, val in cbh.state.items():
            if not "_adj" in key:
                continue
            # Could just use the final row or col
            # print(repr(f'"{key}": ("{repr(val)}"),'))
            assert (
                repr(val) == domain["test_ans"]["preprocess_cbh"]["adj"][key]
            )
    return


# -------------------------------------------------------
# Test forcing adj, compare to prms output
case_list_adj_prms_output = ["drb_2yr"]


def test_cbh_adj_prms_output(domain):
    params = PrmsParameters(domain["dir"] / domain["param_file"])
    cbh = CBH(domain["input_files_dict"])

    @timer
    def adj_timed(params):
        return cbh.adjust(params)

    _ = adj_timed(params)

    for var, var_file in domain["prms_outputs"].items():
        prms_output = load_prms_statscsv(domain["dir"] / var_file)
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
    # print(tmp_file)
    cbh.to_netcdf(tmp_file)
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
    return
