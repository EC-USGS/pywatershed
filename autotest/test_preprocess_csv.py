import os
import pathlib as pl

import netCDF4 as nc
import numpy as np
import pandas as pd

from pynhm.preprocess import CsvFile

max_csv_files = 4


def compare_netcdf(csv, nc_name):
    ds = nc.Dataset(nc_name)

    np_data = csv.data
    hrus = csv.hru_ids
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
    csv = CsvFile(name=list(domain["prms_outputs"].values())[-1])

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


def test_single_csv_to_netcdf(domain):
    path = list(domain["prms_outputs"].values())[-1]
    csv = CsvFile(name=path)

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
    for name, path in domain["prms_outputs"].items():
        if path.suffix in (".csv",):
            csv.add_path(path)
            if imax >= max_csv_files:
                break
            imax += 1

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


def test_multiple_csv_to_netcdf(domain):
    csv = CsvFile()
    basedir = None
    imax = 0
    for name, path in domain["prms_outputs"].items():
        if path.suffix in (".csv",):
            basedir = pl.Path(path.parent)
            csv.add_path(path)
            if imax >= max_csv_files:
                break
            imax += 1

    nc_file = basedir / "multiple_variables.nc"
    csv.to_netcdf(nc_file)

    compare_netcdf(csv, nc_file)

    # clean up netcdf file
    os.remove(nc_file)

    return
