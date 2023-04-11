import pathlib as pl

import numpy as np

from pynhm.utils import data_model as dm

nc_file = pl.Path("../test_data/drb_2yr/prcp.nc")


def test_to_dd_identical():
    dd_xr = dm.xr_ds_to_dd(nc_file)  # pass the file
    dd_nc4 = dm.nc4_ds_to_dd(nc_file)

    enc_xr = dd_xr.pop("encoding")
    enc_nc4 = dd_nc4.pop("encoding")

    np.testing.assert_equal(dd_xr, dd_nc4)


def test_xr_dd_xr():
    import xarray as xr

    # pass the ds
    ds1 = xr.open_dataset(nc_file)
    dd1 = dm.xr_ds_to_dd(ds1)
    ds2 = dm.dd_to_xr_ds(dd1)
    xr.testing.assert_identical(ds1, ds2)

    # pass the file
    ds1 = xr.open_dataset(nc_file)
    dd1 = dm.xr_ds_to_dd(nc_file)
    ds2 = dm.dd_to_xr_ds(dd1)
    xr.testing.assert_identical(ds1, ds2)


# def test_nc4_dd_nc4(tmp_path):
#     # pass the file
#     out_file = pl.Path(tmp_path / "test_nc4_dd_nc4.nc")
#     dm.dd_to_nc4_ds(dm.nc4_ds_to_dd(nc_file), out_file)
#     asdf
#     import xarray as xr

#     ds1 = xr.open_dataset(nc_file)
#     ds2 = xr.open_dataset(out_file)

#     xr.testing.assert_identical(ds1, ds2)
#     # looks like the time variable needs more data.
#     # pass the ds
#     # ds1 = nc4.Dataset(nc_file)


# #     dd1 = dm.xr_ds_to_dd(ds1)
# #     ds2 = dm.dd_to_xr_ds(dd1)
# #     xr.testing.assert_identical(ds1, ds2)
