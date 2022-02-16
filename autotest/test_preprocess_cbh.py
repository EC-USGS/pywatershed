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
ans_dict_cbh_files_to_df = {
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
    assert repr(df.iloc[-1, :]) == ans_dict_cbh_files_to_df[domain][var]
    return


# -------------------------------------------------------
# Test concatenation of individual files
case_list_df_concat = ['drb_2yr']

ans_dict_cbh_df_concat = {
    'drb_2yr': (
        'prcp000    0.00\nprcp001    0.00\nprcp002    0.00\nprcp003    0.00\nprcp004    0.00'
        '\n           ... \ntmin760    1.47\ntmin761    1.49\ntmin762    2.91\ntmin763    4.56\n'
        'tmin764    1.79\n'
        'Name: 1980-12-31 00:00:00, Length: 3060, dtype: float64'), }

@pytest.mark.parametrize("case", case_list_df_concat)
def test_cbh_concat(case):
    df = cbh_files_to_df(cbh_input_dict[case])
    # print(repr(f"'{case}': ('{repr(df.iloc[-1, :])}'),"))
    assert repr(df.iloc[-1, :]) == ans_dict_cbh_df_concat[case]
    return


# -------------------------------------------------------
# Test conversion to numpy dict
case_list_np_dict = ['drb_2yr']

ans_dict_np_dict = {
    "drb_2yr": {
        "datetime": (
            "array([\'1979-01-01T00:00:00.000000000\', \'1979-01-02T00:00:00.000000000\',\n"
            "       \'1979-01-03T00:00:00.000000000\', ...,\n       \'1980-12-29T00:00:00.000000000\'"
            ", \'1980-12-30T00:00:00.000000000\',\n       \'1980-12-31T00:00:00.000000000\'], "
            "dtype=\'datetime64[ns]\')"),
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
# could reuse some of these
case_dict_adj = {
    'drb_2yr': cbh_input_dict['drb_2yr'], }

# @pytest.mark.parametrize("case", list(case_dict_adj.keys()))
# def test_cbh_adj(case):
#     params = PrmsParameters(param_file_dict[case])
#     cbh = CBH(case_dict_adj[case])
#     cbh.adjust(params)
#     return


# -------------------------------------------------------
# # JLM move
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
