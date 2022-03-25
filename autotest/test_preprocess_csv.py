import os
import pathlib as pl

import netCDF4 as nc
import numpy as np
import pandas as pd

from pynhm.preprocess import CsvFile


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
    csv = CsvFile(name=list(domain["prms_outputs"].values())[0])

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


def test_single_csv_to_netcdf(domain):
    files = list(domain["prms_outputs"].values())[0]
    csv = CsvFile(name=files)

    nc_file = pl.Path(files.with_suffix(".nc"))
    csv.to_netcdf(nc_file)

    compare_netcdf(csv, nc_file)

    # clean up netcdf file
    os.remove(nc_file)

    return


def test_multiple_csv(domain):
    csv = CsvFile()
    for name, path in domain["prms_outputs"].items():
        if path.suffix in (".csv",):
            csv.add_path(path)

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


def test_multiple_csv_to_netcdf(domain):
    csv = CsvFile()
    basedir = None
    for name, path in domain["prms_outputs"].items():
        if path.suffix in (".csv",):
            basedir = pl.Path(path.parent)
            csv.add_path(path)

    nc_file = basedir / "multiple_variables.nc"
    csv.to_netcdf(nc_file)

    compare_netcdf(csv, nc_file)

    # clean up netcdf file
    os.remove(nc_file)

    return
