import os
import pathlib as pl

import netCDF4 as nc
import numpy as np
import pandas as pd

from pynhm.preprocess import CsvFile

csv_test_vars = ["hru_ppt", "intcp_stor", "potet", "gwres_stor"]


def compare_netcdf(csv, nc_name):
    ds = nc.Dataset(nc_name)

    np_data = csv.data
    hrus = csv.nhm_id
    _ = csv.nhm_seg
    variables = csv.variable_names

    for variable in variables:
        nc_data = ds[variable][:]
        csv_data = np.zeros(nc_data.shape, dtype=float)
        for idx, hru in enumerate(hrus):
            key = f"{variable}_{hru}"
            csv_data[:, idx] = np_data[key][:]
        assert np.allclose(
            nc_data, csv_data
        ), f"comparision for {variable} is not close"

    # close the netcdf file and return
    ds.close()
    return


def test_single_csv(domain):
    var = "gwres_stor"
    csv = CsvFile(path=domain["prms_output_dir"] / f"{var}.csv")

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


def test_single_csv_to_netcdf(domain):
    var = "gwres_stor"
    path = domain["prms_output_dir"] / f"{var}.csv"
    csv = CsvFile(path=path)

    basedir = pl.Path(path.parent)
    nc_file = basedir / "single_variable.nc"
    csv.to_netcdf(nc_file)

    compare_netcdf(csv, nc_file)

    # clean up netcdf file
    os.remove(nc_file)

    return


def test_multiple_csv(domain):
    csv = CsvFile()
    imax = 0
    for var in csv_test_vars:
        csv.add_path(domain["prms_output_dir"] / f"{var}.csv")
    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


def test_multiple_csv_to_netcdf(domain, tmp_path):
    csv = CsvFile()
    for var in csv_test_vars:
        csv.add_path(domain["prms_output_dir"] / f"{var}.csv")
    nc_file = tmp_path / "multiple_variables.nc"
    csv.to_netcdf(nc_file)
    compare_netcdf(csv, nc_file)
    return
