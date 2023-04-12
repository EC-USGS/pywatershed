import pathlib as pl

import numpy as np
import xarray as xr

from pynhm.utils import data_model as dm

nc_file = pl.Path("../test_data/drb_2yr/prcp.nc")


def test_to_dd_identical():
    dd_xr = dm.xr_ds_to_dd(nc_file)  # pass the file
    dd_nc4 = dm.nc4_ds_to_dd(nc_file)

    enc_xr = dd_xr.pop("encoding")
    enc_nc4 = dd_nc4.pop("encoding")

    np.testing.assert_equal(dd_xr, dd_nc4)


def test_xr_dd_xr(tmp_path):
    tmp_path = pl.Path(tmp_path)

    def assert_encodings_identical(ds1, ds2):
        assert sorted(ds1.variables) == sorted(ds2.variables)
        for vv in ds1.variables:
            np.testing.assert_equal(ds1[vv].encoding, ds2[vv].encoding)

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
    dd1 = dm.xr_dict_to_dd(ds1.to_dict(numpy_data=True, encoding=True))
    dd2 = dm.xr_dict_to_dd(dm.dd_to_xr_dict(dd1))
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
