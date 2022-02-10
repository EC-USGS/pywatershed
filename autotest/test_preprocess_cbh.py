import os
import pathlib as pl  # rm pathlib?
import sys

import pytest

cwd = os.getcwd()
if cwd.endswith("autotest"):
    sys.path.append("..")
    rel_path = ".."
elif cwd.endswith("prmsNHMpy"):
    sys.path.append(".")
    rel_path = "."
from preprocess import cbh_to_df
from pynhm import timer


@timer
def cbh_to_df_timed(*args):
    return cbh_to_df(*args)


def mk_path(*args):
    return os.path.join(base_path, *args)


# establish some input data
base_path = os.path.join(rel_path, "..", "prmsNHMpy", "test_models")
cbh_dict = {
    'merced': {
        'all': mk_path("merced", "input", "mercd.data"), },
    'acf': {
        'orad': mk_path("acf", "input", "ACF_current_Orad.data"),
        'prcp': mk_path("acf", "input", "ACF_current_Prcp.data"),
        'ptet': mk_path("acf", "input", "ACF_current_Ptet.data"),
        'tmax': mk_path("acf", "input", "ACF_current_Tmax.data"),
        'tmin': mk_path("acf", "input", "ACF_current_Tmin.data"), },
    'acfb_water_use': {
        'orad': mk_path("acfb_water_use", "input", "Orad_1950-99.prms"),
        'prcp': mk_path("acfb_water_use", "input", "Prcp_1950-99.prms"),
        'ptet': mk_path("acfb_water_use", "input", "Ptet_1950-99.prms"),
        'tmax': mk_path("acfb_water_use", "input", "Tmax_1950-99.prms"),
        'tmin': mk_path("acfb_water_use", "input", "Tmin_1950-99.prms"), }, }

cbh_ans_dict = {
    'merced': {
        'all': (
            'tmax00       75.0\ntmax01       61.0\ntmax02     -999.0\ntmax03       69.0\n'
            'tmax04     -999.0\n            ...  \nprecip22   -999.0\nprecip23   -999.0\n'
            'precip24   -999.0\nprecip25   -999.0\nrunoff0      95.0\n'
            'Name: 1995-09-30 00:00:00, Length: 79, dtype: float64'), },
    'acf': {
        'orad': (
            'orad000    161.06\norad001    180.02\norad002    138.34\norad003    145.76\n'
            'orad004    152.46\n            ...  \norad253    152.13\norad254    151.77\n'
            'orad255    152.22\norad256    152.01\norad257    151.88\n'
            'Name: 1990-09-30 00:00:00, Length: 258, dtype: float64'),
        'prcp': (
            'prcp000    0.00\nprcp001    0.01\nprcp002    0.01\nprcp003    0.00\n'
            'prcp004    0.00\n           ... \nprcp253    0.43\nprcp254    0.41\n'
            'prcp255    0.22\nprcp256    0.30\nprcp257    0.42\n'
            'Name: 1990-09-30 00:00:00, Length: 258, dtype: float64'),
        'ptet': (
            'ptet000    0.04\nptet001    0.04\nptet002    0.03\nptet003    0.04\n'
            'ptet004    0.04\n           ... \nptet253    0.04\nptet254    0.04\n'
            'ptet255    0.04\nptet256    0.04\nptet257    0.04\n'
            'Name: 1990-09-30 00:00:00, Length: 258, dtype: float64'),
        'tmax': (
            'tmax000    74.92\ntmax001    75.17\ntmax002    76.99\ntmax003    76.16\n'
            'tmax004    76.35\n           ...  \ntmax253    77.33\ntmax254    80.32\n'
            'tmax255    79.07\ntmax256    77.28\ntmax257    80.16\n'
            'Name: 2000-09-30 00:00:00, Length: 258, dtype: float64'),
        'tmin': (
            'tmin000    52.12\ntmin001    52.48\ntmin002    55.40\ntmin003    54.53\n'
            'tmin004    54.86\n           ...  \ntmin253    67.18\ntmin254    68.12\n'
            'tmin255    67.78\ntmin256    67.15\ntmin257    62.98\n'
            'Name: 2000-09-30 00:00:00, Length: 258, dtype: float64'), },
    'acfb_water_use': {
        'orad': (
            'orad00    358.36\norad01    437.09\norad02    224.55\norad03    280.73\n'
            'orad04    120.42\norad05    346.94\norad06    137.04\norad07    303.29\n'
            'orad08    284.72\norad09    314.25\norad10    120.26\norad11    288.21\n'
            'orad12    269.50\norad13    211.57\norad14    255.46\norad15    233.39\n'
            'orad16    262.19\norad17    292.65\norad18    309.17\norad19    292.32\n'
            'orad20    310.13\norad21    275.37\norad22    254.89\norad23    308.33\n'
            'orad24    323.41\norad25    210.54\norad26    250.98\n'
            'Name: 1999-12-31 00:00:00, dtype: float64'),
        'prcp': (
            'precip00    0.00\nprecip01    0.00\nprecip02    0.00\nprecip03    0.00\n'
            'precip04    0.01\nprecip05    0.00\nprecip06    0.01\nprecip07    0.00\n'
            'precip08    0.00\nprecip09    0.00\nprecip10    0.01\nprecip11    0.00\n'
            'precip12    0.00\nprecip13    0.00\nprecip14    0.00\nprecip15    0.00\n'
            'precip16    0.00\nprecip17    0.00\nprecip18    0.00\nprecip19    0.00\n'
            'precip20    0.00\nprecip21    0.00\nprecip22    0.00\nprecip23    0.00\n'
            'precip24    0.00\nprecip25    0.00\nprecip26    0.00\n'
            'Name: 1999-12-31 00:00:00, dtype: float64'),
        'ptet': (
            'ptet00    0.06\nptet01    0.07\nptet02    0.03\nptet03    0.05\n'
            'ptet04    0.02\nptet05    0.05\nptet06    0.02\nptet07    0.05\n'
            'ptet08    0.05\nptet09    0.05\nptet10    0.02\nptet11    0.04\n'
            'ptet12    0.05\nptet13    0.03\nptet14    0.04\nptet15    0.04\n'
            'ptet16    0.04\nptet17    0.05\nptet18    0.05\nptet19    0.05\n'
            'ptet20    0.05\nptet21    0.05\nptet22    0.04\nptet23    0.06\n'
            'ptet24    0.06\nptet25    0.03\nptet26    0.05\n'
            'Name: 1999-12-31 00:00:00, dtype: float64'),
        'tmax': (
            'tmax00    62.04\ntmax01    60.72\ntmax02    60.02\ntmax03    62.83\n'
            'tmax04    54.69\ntmax05    58.70\ntmax06    53.39\ntmax07    60.80\n'
            'tmax08    60.61\ntmax09    58.51\ntmax10    55.11\ntmax11    57.19\n'
            'tmax12    61.05\ntmax13    59.13\ntmax14    57.07\ntmax15    59.85\n'
            'tmax16    59.59\ntmax17    58.19\ntmax18    59.41\ntmax19    59.86\n'
            'tmax20    60.00\ntmax21    59.36\ntmax22    60.32\ntmax23    62.14\n'
            'tmax24    59.59\ntmax25    59.60\ntmax26    62.80\n'
            'Name: 1999-12-31 00:00:00, dtype: float64'),
        'tmin': (
            'tmin00    31.23\ntmin01    28.21\ntmin02    27.56\ntmin03    31.97\n'
            'tmin04    23.81\ntmin05    27.95\ntmin06    23.30\ntmin07    31.78\n'
            'tmin08    32.42\ntmin09    30.23\ntmin10    25.87\ntmin11    29.08\n'
            'tmin12    32.69\ntmin13    32.16\ntmin14    29.84\ntmin15    33.47\n'
            'tmin16    32.25\ntmin17    31.85\ntmin18    32.61\ntmin19    33.10\n'
            'tmin20    31.99\ntmin21    31.92\ntmin22    32.14\ntmin23    34.52\n'
            'tmin24    31.73\ntmin25    31.74\ntmin26    35.21\n'
            'Name: 1999-12-31 00:00:00, dtype: float64')}, }

cbh_case_list = []
for domain, file_dict in cbh_dict.items():
    for var, file in file_dict.items():
        cbh_case_list += [(domain, var)]


@pytest.mark.parametrize("domain, var", cbh_case_list)
def test_preprocess_cbh(domain, var):
    assert pl.Path(file).exists(), file  # rm pathlib?
    the_file = cbh_dict[domain][var]
    df = cbh_to_df_timed(the_file)
    # print(repr(f"'{var}': ('{repr(df.iloc[-1, :])}'),"))
    assert repr(df.iloc[-1, :]) == cbh_ans_dict[domain][var]
    return


if __name__ == "__main__":

    for idx, case in cbh_case_list:
        domain, var = case
        test_preprocess_cbh(domain, var)
