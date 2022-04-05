import os
import pathlib as pl

import netCDF4 as nc
import numpy as np

from pynhm.utils.netcdf_utils import NetCdfRead

def test_netcdf(domain):
    output_files = list(domain["prms_outputs"].values())
    output_path = output_files[0]
    nc_pth = output_path.with_suffix(".nc")
    variable = nc_pth.stem

    nc_data = NetCdfRead(nc_pth)

    for idx in range(nc_data.ntimes):
        arr = nc_data.advance(variable)