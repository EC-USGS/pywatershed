import numpy as np
import os
import pathlib as pl
import sys

import pytest

cwd = os.getcwd()
if cwd.endswith("autotest"):
    sys.path.append("..")
    rel_path = pl.Path("..")
elif cwd.endswith("pynhm"):
    sys.path.append(".")
    rel_path = pl.Path(".")

from pynhm.preprocess.cbh import CBH
from pynhm.utils import timer
from pynhm.prms5util import load_prms_statscsv
from pynhm import PrmsParameters
from pynhm.preprocess.cbh_utils import (
    cbh_files_to_df, cbh_files_to_np_dict)

# establish some input data


# helper
def mk_path(arg):
    return base_path / arg


@timer
def cbh_files_to_df_timed(*args):
    df = cbh_files_to_df(*args)
    return df


base_path = rel_path / ("test_data")

cbh_input_dict = {
    'drb_2yr': {
        'prcp': mk_path("drb_2yr/prcp.cbh"),
        'rhavg': mk_path("drb_2yr/rhavg.cbh"),
        'tmax': mk_path("drb_2yr/tmax.cbh"),
        'tmin': mk_path("drb_2yr/tmin.cbh"),},}

param_file_dict = {
    'drb_2yr': mk_path("drb_2yr/myparam.param"), }

prms_output_dict = {
    'drb_2yr': {
        'prcp_adj': mk_path("drb_2yr/output/nhru_hru_ppt.csv"),
        'rainfall_adj': mk_path("drb_2yr/output/nhru_hru_rain.csv"),
        'snowfall_adj': mk_path("drb_2yr/output/nhru_hru_snow.csv"),
        'tmax_adj': mk_path("drb_2yr/output/nhru_tmaxf.csv"),
        'tmin_adj': mk_path("drb_2yr/output/nhru_tminf.csv"), }, }

# -------------------------------------------------------
# Test Instantiation from a string vs dict vs other(list)
case_dict_init = {
    'drb_2yr': cbh_input_dict['drb_2yr'],
    'invalid': ['a', 'b'], }


@pytest.mark.parametrize("case", list(case_dict_init.keys()))
def test_cbh_init_files(case):
    try:
        _ = CBH(case_dict_init[case])
        assert True
    except ValueError:
        if case in ['invalid']:
            assert True
        else:
            assert False
    return


# -------------------------------------------------------
# Test the file to df parser on a *single file*
# Check the repr of the final row
ans_dict_files_to_df = {
    'drb_2yr': {
        'prcp': (
            'prcp000    0.0\nprcp001    0.0\nprcp002    0.0\nprcp003    0.0\nprcp004    0.0'
            '\n          ... \nprcp760    0.0\nprcp761    0.0\nprcp762    0.0\nprcp763    0.0'
            '\nprcp764    0.0\n'
            'Name: 1980-12-31 00:00:00, Length: 765, dtype: float64'),
        'rhavg': (
            'rhavg000    49.11\nrhavg001    48.35\nrhavg002    48.75\nrhavg003    50.05\n'
            'rhavg004    48.90\n            ...  \nrhavg760    52.08\nrhavg761    54.73\n'
            'rhavg762    53.93\nrhavg763    56.40\nrhavg764    50.92\n'
            'Name: 1980-12-31 00:00:00, Length: 765, dtype: float64'),
        'tmax': (
            'tmax000    40.95\ntmax001    40.55\ntmax002    40.93\ntmax003    38.87\n'
            'tmax004    39.57\n           ...  \ntmax760    18.15\ntmax761    18.92\n'
            'tmax762    19.65\ntmax763    19.93\ntmax764    16.13\n'
            'Name: 1980-12-31 00:00:00, Length: 765, dtype: float64'),
        'tmin': (
            'tmin000    26.88\ntmin001    26.51\ntmin002    26.78\ntmin003    25.53\n'
            'tmin004    25.59\n           ...  \ntmin760     1.47\ntmin761     1.49\n'
            'tmin762     2.91\ntmin763     4.56\ntmin764     1.79\n'
            'Name: 1980-12-31 00:00:00, Length: 765, dtype: float64'), }, }

case_list_file_to_df = [
    ('drb_2yr', 'prcp'),
    ('drb_2yr', 'rhavg'),
    ('drb_2yr', 'tmax'),
    ('drb_2yr', 'tmin'), ]


@pytest.mark.parametrize("domain, var", case_list_file_to_df)
def test_cbh_files_to_df(domain, var):
    the_file = cbh_input_dict[domain][var]
    assert pl.Path(the_file).exists(), the_file  # rm pathlib?
    df = cbh_files_to_df_timed(the_file)
    #print(repr(f"'{var}': ('{repr(df.iloc[-1, :])}'),"))
    assert repr(df.iloc[-1, :]) == ans_dict_files_to_df[domain][var]
    return


# -------------------------------------------------------
# Test concatenation of individual files
case_list_df_concat = ['drb_2yr']

ans_dict_df_concat = {
    'drb_2yr': (
        'prcp000    0.00\nprcp001    0.00\nprcp002    0.00\nprcp003    0.00\nprcp004    0.00'
        '\n           ... \ntmin760    1.47\ntmin761    1.49\ntmin762    2.91\ntmin763    4.56\n'
        'tmin764    1.79\n'
        'Name: 1980-12-31 00:00:00, Length: 3060, dtype: float64'), }

@pytest.mark.parametrize("case", case_list_df_concat)
def test_cbh_concat(case):
    df = cbh_files_to_df(cbh_input_dict[case])
    # print(repr(f"'{case}': ('{repr(df.iloc[-1, :])}'),"))
    assert repr(df.iloc[-1, :]) == ans_dict_df_concat[case]
    return


# -------------------------------------------------------
# Test conversion to numpy dict
case_list_np_dict = ['drb_2yr']

ans_dict_np_dict = {
    "drb_2yr": {
        "datetime": (
            "array([\'1979-01-01\', \'1979-01-02\', \'1979-01-03\', ..., \'1980-12-29\',\n"
            "       \'1980-12-30\', \'1980-12-31\'], "
            "dtype=\'datetime64[D]\')"),
        "prcp": (
            "array([[0.48, 0.48, 0.48, ..., 0.88, 0.79, 0.97],\n       [1.01, 1.01, 1.01, ..., 1.15, 1.04, 1.45],\n"
            "       [0.  , 0.  , 0.  , ..., 0.06, 0.05, 0.09],\n       ...,\n       [0.12, 0.05, 0.09, ..., 0.  , 0.  , 0.  ],\n"
            "       [0.  , 0.  , 0.  , ..., 0.  , 0.  , 0.  ],\n       [0.  , 0.  , 0.  , ..., 0.  , 0.  , 0.  ]])"),
        "rhavg": (
            "array([[82.46, 82.6 , 83.03, ..., 66.93, 68.67, 64.31],\n       [81.99, 82.35, 82.43, ..., 68.38, 69.92, 66.48],\n"
            "       [68.49, 70.05, 68.94, ..., 61.18, 62.23, 58.88],\n       ...,\n"
            "       [76.27, 71.1 , 75.26, ..., 66.49, 69.15, 63.22],\n       [65.54, 64.7 , 65.11, ..., 63.96, 66.34, 61.66],\n"
            "       [49.11, 48.35, 48.75, ..., 53.93, 56.4 , 50.92]])"),
        "tmax": (
            "array([[56.23, 55.13, 55.26, ..., 51.23, 49.8 , 47.29],\n       [54.33, 53.69, 53.82, ..., 50.04, 48.84, 46.2 ],\n"
            "       [41.57, 41.27, 41.39, ..., 26.76, 28.  , 25.96],\n       ...,\n"
            "       [47.08, 46.67, 46.8 , ..., 41.94, 41.  , 38.15],\n       [45.45, 44.87, 45.  , ..., 35.14, 35.2 , 31.42],\n"
            "       [40.95, 40.55, 40.93, ..., 19.65, 19.93, 16.13]])"),
        "tmin": (
            "array([[44.53, 44.33, 44.33, ..., 40.21, 38.67, 36.  ],\n       [36.17, 36.05, 35.93, ..., 26.49, 27.66, 25.15],\n"
            "       [14.97, 15.17, 15.05, ...,  2.75,  3.92,  0.88],\n       ...,\n"
            "       [35.3 , 35.33, 35.33, ..., 27.93, 28.07, 24.7 ],\n       [25.49, 25.43, 25.43, ...,  1.  ,  2.94, -0.68],\n"
            "       [26.88, 26.51, 26.78, ...,  2.91,  4.56,  1.79]])"), }, }

@pytest.mark.parametrize("case", case_list_np_dict)
def test_cbh_np_dict(case):
    cbh = CBH(cbh_input_dict[case])
    with np.printoptions(threshold=8):
        for key, val in cbh.state.items():
            # print(repr(f'"{key}": ("{repr(val)}"),'))
            assert repr(val) == ans_dict_np_dict[case][key]
    return


# -------------------------------------------------------
# Test forcing adjustment
case_list_adj = ['drb_2yr']

ans_dict_adj = {
    'drb_2yr': {
        "tmax_adj": (
            "array([[57.41188 , 56.47271 , 57.6613  , ..., 51.909633, 51.564279,\n        50.29    ],\n"
            "       [55.51188 , 55.03271 , 56.2213  , ..., 50.719633, 50.604279,\n        49.2     ],\n"
            "       [42.75188 , 42.61271 , 43.7913  , ..., 27.439633, 29.764279,\n        28.96    ],\n"
            "       ...,\n       [48.26188 , 48.01271 , 49.2013  , ..., 42.619633, 42.764279,\n        41.15    ],\n"
            "       [46.63188 , 46.21271 , 47.4013  , ..., 35.819633, 36.964279,\n        34.42    ],\n"
            "       [42.13188 , 41.89271 , 43.3313  , ..., 20.329633, 21.694279,\n        19.13    ]])"),
        "tmin_adj": (
            "array([[46.12091 , 45.76805 , 44.46652 , ..., 38.069587, 37.981065,\n        36.876292],\n"
            "       [37.76091 , 37.48805 , 36.06652 , ..., 24.349587, 26.971065,\n        26.026292],\n"
            "       [16.56091 , 16.60805 , 15.18652 , ...,  0.609587,  3.231065,\n         1.756292],\n"
            "       ...,\n       [36.89091 , 36.821522, 35.46652 , ..., 25.789587, 27.381065,\n        25.586652],\n"
            "       [27.08091 , 26.921522, 25.56652 , ..., -1.140413,  2.251065,\n         0.206652],\n"
            "       [28.47091 , 28.001522, 26.91652 , ...,  0.769587,  3.871065,\n         2.676652]])"),
        "prcp_adj": (
            "array([[0., 0., 0., ..., 0., 0., 0.],\n       [0., 0., 0., ..., 0., 0., 0.],\n"
            "       [0., 0., 0., ..., 0., 0., 0.],\n       ...,\n       [0., 0., 0., ..., 0., 0., 0.],\n"
            "       [0., 0., 0., ..., 0., 0., 0.],\n       [0., 0., 0., ..., 0., 0., 0.]])"),
        "rainfall_adj": (
            "array([[0.3139296 , 0.2478048 , 0.248112  , ..., 0.48323264, 0.57291748,\n        0.49815223],\n"
            "       [0.6605602 , 0.5214226 , 0.522069  , ..., 0.6314972 , 0.75422048,\n        0.74466055],\n"
            "       [0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ],\n"
            "       ...,\n       [0.0784824 , 0.025813  , 0.046521  , ..., 0.        , 0.        ,\n        0.        ],\n"
            "       [0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ],\n"
            "       [0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ]])"),
        "snowfall_adj": (
            "array([[0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ],\n"
            "       [0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ],\n"
            "       [0.        , 0.        , 0.        , ..., 0.03856368, 0.025     ,\n        0.04753728],\n"
            "       ...,\n       [0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ],\n"
            "       [0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ],\n"
            "       [0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ]])"),
        "precip_adj": (
            "array([[0.3139296 , 0.2478048 , 0.248112  , ..., 0.48323264, 0.57291748,\n        0.49815223],\n"
            "       [0.6605602 , 0.5214226 , 0.522069  , ..., 0.6314972 , 0.75422048,\n        0.74466055],\n"
            "       [0.        , 0.        , 0.        , ..., 0.03856368, 0.025     ,\n        0.04753728],\n"
            "       ...,\n       [0.0784824 , 0.025813  , 0.046521  , ..., 0.        , 0.        ,\n        0.        ],\n"
            "       [0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ],\n"
            "       [0.        , 0.        , 0.        , ..., 0.        , 0.        ,\n        0.        ]])"), }, }


@pytest.mark.parametrize("case", case_list_adj)
def test_cbh_adj(case):
    params = PrmsParameters(param_file_dict[case])
    cbh = CBH(cbh_input_dict[case])
    _ = cbh.adjust(params)
    with np.printoptions(threshold=2):
        for key, val in cbh.state.items():
            if not '_adj' in key:
                continue
            # Could just use the final row or col
            # print(repr(f'"{key}": ("{repr(val)}"),'))
            assert repr(val) == ans_dict_adj[case][key]
    return


# -------------------------------------------------------
# Test forcing adj, compare to prms output
case_list_adj_prms_output = ['drb_2yr']


@pytest.mark.parametrize("case", case_list_adj_prms_output)
def test_cbh_adj_prms_output(case):
    params = PrmsParameters(param_file_dict[case])
    cbh = CBH(cbh_input_dict[case])
    @timer
    def adj_timed(params):
        return cbh.adjust(params)
    _ = adj_timed(params)

    for var, var_file in prms_output_dict[case].items():
        prms_output = load_prms_statscsv(var_file)
        p_dates = prms_output.index.values
        p_array = prms_output.to_numpy()
        wh_dates = np.where(np.isin(cbh.state['datetime'], p_dates))
        print(var)
        result = np.isclose(
            p_array,
            cbh.state[var][wh_dates,:],
            rtol=1e-121, atol=1e-04  # Only the atol matters here, if atol < 1e-4 fails
        )
        print(result.sum()/np.prod(result.shape))
        print(result.all())

    return


# -------------------------------------------------------
# # Test netcdf write
# case_dict_netcdf = {
#     # 0: cbh_input_dict['merced']['all'],
#     1: {key: val for key, val in cbh_input_dict['acf'].items()
#         if key in ['prcp', 'tmax', 'tmin']},
# }


# @pytest.mark.parametrize("case", list(case_dict.keys()))
# def test_cbh_to_netcdf(case):
#     cbh = CBH(case_dict[case])
#     file = 'dummy.nc'
#     cbh.to_netcdf(file)
#     return
