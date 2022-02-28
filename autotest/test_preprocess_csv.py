import os
import pathlib as pl
import sys

import pandas as pd
import pytest

from pynhm.preprocess import CsvFile

# JLM: Need to add/configure box01 domain yaml


def test_single_csv(domain):
    csv = CsvFile(name=list(domain["prms_outputs"].values())[0])

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return


def test_multiple_csv(domain):
    csv = CsvFile()
    for name, path in domain["prms_outputs"].items():
        csv.add_path(path)

    df = csv.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    return
