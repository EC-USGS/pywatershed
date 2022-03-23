import os
import pathlib as pl

import pandas as pd

from pynhm.preprocess import CsvFile

# JLM: Need to add/configure box01 domain yaml


def test_single_csv(domain):
    csv = CsvFile(name=list(domain["prms_outputs"].values())[0])

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


def test_single_csv_to_netcdf(domain):
    files = list(domain["prms_outputs"].values())[0]
    csv = CsvFile(name=files)

    csv_file = pl.Path(files.with_suffix(".nc"))
    csv.to_netcdf(csv_file)

    return


def test_multiple_csv(domain):
    csv = CsvFile()
    for name, path in domain["prms_outputs"].items():
        if path.suffix in (".csv",):
            csv.add_path(path)

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return
