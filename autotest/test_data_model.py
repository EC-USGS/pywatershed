import pathlib as pl

import netCDF4 as nc4
import numpy as np
import pandas as pd
import xarray as xr

from pynhm.base import data_model as dm
from pynhm.base.data_model import DatasetDict


nc_file = pl.Path("../test_data/drb_2yr/prcp.nc")


def mk_ds0():
    # a simple dataset taken from xarray Dataset documentation
    np.random.seed(0)
    temperature = 15 + 8 * np.random.randn(2, 2, 3)
    precipitation = 10 * np.random.rand(2, 2, 3)
    lon = [[-99.83, -99.32], [-99.79, -99.23]]
    lat = [[42.25, 42.21], [42.63, 42.59]]
    time = pd.date_range("2014-09-06", periods=3)
    reference_time = pd.Timestamp("2014-09-05")

    ds = xr.Dataset(
        data_vars=dict(
            temperature=(["x", "y", "time"], temperature),
            precipitation=(["x", "y", "time"], precipitation),
        ),
        coords=dict(
            lon=(["x", "y"], lon),
            lat=(["x", "y"], lat),
            time=time,
            reference_time=reference_time,
        ),
        attrs=dict(description="Weather related data."),
    )

    return ds


ds0 = mk_ds0()
xr_dd0 = ds0.to_dict(numpy_data=True, encoding=True)
dd0 = DatasetDict.from_ds(ds0)


def assert_encodings_identical(ds1, ds2):
    assert sorted(ds1.variables) == sorted(ds2.variables)
    for vv in ds1.variables:
        np.testing.assert_equal(ds1[vv].encoding, ds2[vv].encoding)


# helpers and data above
# -------------------------------------------------------
# tests below


def test_to_dd_identical():
    dd_xr = dm.xr_ds_to_dd(nc_file)  # pass the file
    dd_nc4 = dm.nc4_ds_to_dd(nc_file)

    enc_xr = dd_xr.pop("encoding")
    enc_nc4 = dd_nc4.pop("encoding")

    np.testing.assert_equal(dd_xr, dd_nc4)


def test_xr_dd_xr(tmp_path):
    tmp_path = pl.Path(tmp_path)

    # start with ds loaded: ds -> dd -> ds
    ds1 = xr.open_dataset(nc_file)
    ds2 = dm.dd_to_xr_ds(dm.xr_ds_to_dd(ds1))  # round trip from ds
    xr.testing.assert_identical(ds1, ds2)
    assert_encodings_identical(ds1, ds2)

    # start with ds as file: file -> ds -> dd -> ds
    ds3 = dm.dd_to_xr_ds(dm.xr_ds_to_dd(nc_file))  # round trip from file
    xr.testing.assert_identical(ds1, ds3)
    assert_encodings_identical(ds1, ds3)

    # start with a dd
    # round trip through a xr.dict: dd -> xr.dict -> dd
    dd1 = dm.xr_dd_to_dd(ds1.to_dict(numpy_data=True, encoding=True))
    dd2 = dm.xr_dd_to_dd(dm.dd_to_xr_dd(dd1))
    np.testing.assert_equal(dd1, dd2)

    # roundtrip a dd though xr.ds: dd -> ds -> dd
    dd3 = dm.xr_ds_to_dd(dm.dd_to_xr_ds(dd1))
    np.testing.assert_equal(dd1, dd3)

    # roundtrip a dd though netcdf file: dd -> ds -> netcdf file -> ds -> dd
    ff = tmp_path / "ds.nc"
    _ = dm.dd_to_xr_ds(dd1).to_netcdf(ff)
    dd4 = dm.xr_ds_to_dd(xr.open_dataset(ff))
    # The source file is obviously different
    for kk in dd1["encoding"].keys():
        del dd1["encoding"][kk]["source"]
        del dd4["encoding"][kk]["source"]

    return


def test_nc4_dd_nc4(tmp_path):
    # pass the file and round trip it to file, compare via xarry
    out_file = pl.Path(tmp_path / "rt_1.nc")
    _ = dm.dd_to_nc4_ds(dm.nc4_ds_to_dd(nc_file), out_file)
    ds1 = xr.open_dataset(nc_file)
    ds2 = xr.open_dataset(out_file)
    xr.testing.assert_identical(ds1, ds2)

    # the source file was not the same
    # coordinates encoded on the coordinate variables in netcdf
    for kk in ["source", "coordinates"]:
        for vv in sorted(ds1.variables):
            for dd in [ds1, ds2]:
                if kk in dd[vv].encoding.keys():
                    del dd[vv].encoding[kk]

    assert_encodings_identical(ds1, ds2)

    # test passing a nc4.Dataset, cant not use_xr_encodings
    # round trip via file and compare in xarray
    out_file = pl.Path(tmp_path / "rt_2.nc")
    nc_file_ds = nc4.Dataset(nc_file)
    _ = dm.dd_to_nc4_ds(
        dm.nc4_ds_to_dd(nc_file_ds, use_xr_enc=False), out_file
    )
    ds1 = xr.open_dataset(nc_file)
    ds2 = xr.open_dataset(out_file)
    xr.testing.assert_identical(ds1, ds2)

    # nc4 encoding defaults dont match the source since we didnt read them
    # with xarray
    for kk in [
        "source",
        "zlib",
        "complevel",
        "shuffle",
        "contiguous",
        "chunksizes",
        "coordinates",
    ]:
        for vv in sorted(ds1.variables):
            for dd in [ds1, ds2]:
                if kk in dd[vv].encoding.keys():
                    del dd[vv].encoding[kk]

    assert_encodings_identical(ds1, ds2)

    return


def dd_basics():
    # init
    # properties and accessors
    #  - return dicts are still editable, which is dangerous?
    # repr
    # str
    pass


def dd_validate():
    # top-level keys
    # coord and var dims
    pass


def dd_spatial_coord_names():
    pass


def dd_netcdf(tmp_path):
    tmp_path = pl.Path(tmp_path)
    dd0 = DatasetDict.from_ds(ds0)  # from xarray ds
    dd0_file = tmp_path / "dd0.nc"
    dd0.to_netcdf(dd0_file)  # to file via xr (default)
    ds1 = xr.open_dataset(dd0_file)
    xr.testing.assert_identical(ds1, ds0)  # ds -> file(via xr) -> ds
    dd1 = DatasetDict.from_ds(ds1)
    del ds1
    dd0.to_netcdf(dd0_file, use_xr=False)  # file via nc4
    ds1 = xr.open_dataset(dd0_file)
    xr.testing.assert_identical(ds1, ds0)  # ds -> file(via nc4) -> ds
    dd1_nc4 = DatasetDict.from_ds(ds1)
    np.testing.assert_equal(dd1.data, dd1_nc4.data)
    del ds1

    # read using nc4
    # ds = nc4.
    np.testing.assert_equal(dd1.data, dd1_nc4.data)

    # from
    # to
    # from
    # roundtrip
    # same as round trips above?
    return


def dd_subset():
    pass


def dd_get_parameters():
    pass


def dd_merge():
    dd0 = {
        "dims": {},
        "coords": {},
    }
    pass


# method to get dd as a dict
# methods to build dds from data and coord_name or coord_data
