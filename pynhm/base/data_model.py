from copy import deepcopy
import itertools
import pprint

import netCDF4 as nc4
import numpy as np
import xarray as xr

from .accessor import Accessor
from ..constants import listish

# This file defines the data model for pywatershed. It is called a
# "dataset_dict" and has a invertible mapping with non-hierarchical netcdf
# or xarray datasets.

# Show a basic example of how an xr_dd is changed to a dd.
# What is the difference


# TODO: what about hierarchical/groups in netcdf files?


# This is what a dataset_dict looks like. Metadata for coord and data_vars
# is found in metadata.
template_dd = {
    "dims": {},
    "coords": {},
    "data_vars": {},
    "metadata": {},
    "encoding": {},
}

template_xr_dd = {
    "attrs": {},
    "dims": {},
    "coords": {},
    "data_vars": {},
    "encoding": {},
}


class DatasetDict(Accessor):
    def __init__(
        self,
        dims: dict = {},
        coords: dict = {},
        data_vars: dict = {},
        metadata: dict = {},
        encoding: dict = {},
    ) -> "DatasetDict":
        self._data_vars = data_vars
        self._dims = dims
        self._coords = coords
        self._metadata = metadata
        self._encoding = encoding

        return

    @property
    def dims(self) -> dict:
        """Return the dimensions"""
        return self._dims

    @property
    def coords(self) -> dict:
        """Return the coordinates"""
        return self._coords

    @property
    def data_vars(self) -> dict:
        """Return the data_vars."""
        return self._data_vars

    @property
    def variables(self) -> dict:
        """Return coords and data_vars together"""
        return {**self._coords, **self._data_vars}

    @property
    def metadata(self) -> dict:
        """Return the metadata"""
        return self._metadata

    @property
    def encoding(self) -> dict:
        """Return the encoding"""
        return self._encoding

    @property
    def data(self) -> dict:
        return {
            "dims": self._dims,
            "coords": self._coords,
            "data_vars": self._data_vars,
            "metadata": self._metadata,
            "encoding": self._encoding,
        }

    def _keys(self) -> list:
        return ["dims", "coords", "data_vars", "metadata", "encoding"]

    # are __repr__ and __str__ better than default?
    # def __repr__(self):
    #     return pprint.pformat(
    #         {
    #             "dims": self.dims,
    #             "coords": self.coords,
    #             "data_vars": self.data_vars,
    #             "metadata": self.metadata,
    #             "encoding": self.encoding,
    #         }
    #     )

    # def __str__(self):
    #     return

    @staticmethod
    def from_dict(dict_in):
        return DatasetDict(**dict_in)

    @property
    def spatial_coord_names(self) -> dict:
        """Return the spatial coordinate names."""
        attrs = self._metadata["global"]
        return {kk: vv for kk, vv in attrs.items() if "spatial" in kk}

    @staticmethod
    def from_ds(ds):
        # detect typ as xr or nc4
        if isinstance(ds, xr.Dataset):
            return DatasetDict(**xr_ds_to_dd(ds))
        elif isinstance(ds, nc4.Dataset):
            return DatasetDict(**nc4_ds_to_dd(ds))
        else:
            raise ValueError("Passed dataset neither from xarray nor netCDF4")

    @staticmethod
    def from_netcdf(nc_file_list) -> "DatasetDict":
        """Load parameters object from a netcdf file(s?)"""
        # Provide in base class
        # handle more than one file? # see prms ?
        raise NotImplementedError

    def to_xr_ds(self) -> xr.Dataset:
        return dd_to_xr_ds(self.data)

    def to_xr_dd(self) -> dict:
        return dd_to_xr_dd(self.data)

    def to_nc4_ds(self, filename) -> None:
        return dd_to_nc4_ds(self.data, filename)

    def to_netcdf(self, filename, use_xr=False) -> None:
        """Write parameters to a netcdf file"""
        if use_xr:
            self.to_xr_ds.to_netcdf(filename)
        else:
            self.to_nc4_ds(filename)
        return

    def subset(
        self,
        keys: listish = None,
        # process: str = None,
    ) -> "DatasetDict":
        """Returns a Parameters object with a subset of the data."""
        subset = {}
        for dd in self._data_vars.keys():
            if dd not in keys:
                continue
            subset[dd] = self._data_vars[kk]
            asdf

    def get_parameters(
        self, keys: listish, process: str = None, dims: bool = False
    ) -> dict:
        """Returns Parameter data for requested keys."""
        # Provide in base class
        raise NotImplementedError

    def validate(self):
        # TODO:
        #    * check for all required keys
        #    * check metadata keys against data
        #    * check all dims, coords
        raise NotImplementedError

    @staticmethod
    def merge(*dd_list):
        # check each dd's data_vars for duplicates across the list
        all_keys = [list(dd["data_vars"].keys()) for dd in dd_list]
        all_keys = np.array(list(itertools.chain(*all_keys)))

        # def merge_dicts(dict_list):
        #     result = {}
        #     for d in dict_list:
        #         for key, value in d.items():
        #             if key not in result:
        #                 result[key] = value
        #             elif isinstance(value, dict) and isinstance(result[key], dict):
        #                 result[key] = merge_dicts([result[key], value])
        #             elif value == result[key]:
        #                 pass
        #             else:
        #                 raise ValueError(f"Duplicate key '{key}' with non-identical data")
        #     return result

        # if there are duplicate vars, check for their data equality
        raise NotImplementedError

    def _merge_dicts(*dict_list):
        merged_dict = {}

        for d in dict_list:
            for key in d:
                if key in merged_dict:
                    raise ValueError(
                        f"Duplicate key '{key}' found in dictionaries"
                    )
                merged_dict[key] = d[key]

        return merged_dict


# module scope function


def xr_ds_to_dd(file_or_ds, schema_only=False) -> dict:
    """Xarray dataset to a pywatershed dataset dict

    The pyws data model moves metadata off the variables to a separate
    var_metadata dictionary with the same keys found in the union of
    the keys of coords and data_vars.
    """
    if not isinstance(file_or_ds, xr.Dataset):
        xr_ds = xr.open_dataset(file_or_ds)
    else:
        xr_ds = file_or_ds

    dd = xr_ds.to_dict(numpy_data=True, data=(not schema_only), encoding=True)

    dd = xr_dd_to_dd(dd)

    return dd


def xr_dd_to_dd(xr_dd: dict) -> dict:
    dd = deepcopy(xr_dd)

    # Move the global encoding to a global key of itself
    dd["encoding"] = {"global": dd.get("encoding", {})}

    # rename data to metadata for var and coord
    var_metadata = dd.pop("data_vars")
    coord_metadata = dd.pop("coords")

    # create empty data dicts and move the data out of the metadata
    # and move the encoding to encoding[var]
    dd["data_vars"] = {}
    for key, val in var_metadata.items():
        dd["data_vars"][key] = val.pop("data")
        dd["encoding"][key] = val.pop("encoding")

    dd["coords"] = {}
    for key, val in coord_metadata.items():
        dd["coords"][key] = val.pop("data")
        dd["encoding"][key] = val.pop("encoding")

    dd["metadata"] = {**coord_metadata, **var_metadata}
    dd["metadata"]["global"] = dd.pop("attrs")

    return dd


def dd_to_xr_dd(dd: dict) -> dict:
    dd = deepcopy(dd)

    # remove metadata and encoding from the dict/model
    meta = dd.pop("metadata")
    encoding = dd.pop("encoding")

    # loop over meta data, putting it back on the data_vars and/or coords
    # and moving the data to "data"
    dd["attrs"] = meta.pop("global", {})
    dd["encoding"] = encoding.pop("global", {})

    for key, val in meta.items():
        cv = None
        if key in dd["data_vars"].keys():
            cv = "data_vars"
        elif key in dd["coords"].keys():
            cv = "coords"

        dd[cv][key] = {
            **val,
            "data": dd[cv][key],
            "encoding": encoding[key],
        }

    return dd


def dd_to_xr_ds(dd: dict) -> xr.Dataset:
    """pywatershed dataset dict to xarray dataset

    The pyws data model moves metadata off the variables to a separate
    var_metadata dictionary with the same keys found in the union of
    the keys of coords and data_vars. This maps the metadata back to the
    variables.
    """
    return xr.Dataset.from_dict(dd_to_xr_dd(dd))


def _nc4_var_to_datetime64(var, attrs, encoding):
    """netCDF4 conversion of time to numpy.datetime64 based on metadata."""
    if not (hasattr(var, "units") and "since" in var.units):
        return var[:], attrs, encoding

    # Check if the variable has a calendar attribute
    if hasattr(var, "calendar"):
        time_data = nc4.num2date(
            var[:], var.units, var.calendar, only_use_cftime_datetimes=False
        )
    else:
        time_data = nc4.num2date(var[:], var.units)

    time_data = time_data.filled().astype("datetime64[us]")

    for aa in ["calendar", "units"]:
        if aa in attrs.keys():
            encoding[aa] = attrs.pop(aa)

    return time_data, attrs, encoding


def _datetime64_to_nc4_var(var):
    """vectorized conversion of numpy.datetime64 to a netcdf4 variable."""
    # Based on what xarray does
    # https://github.com/pydata/xarray/blob/
    #   a1f5245a48146bd8fc5bdb07ef8ae6077d6e511c/xarray/coding/times.py#L687
    from xarray.coding.times import encode_cf_datetime, decode_cf_datetime

    # This can take encoding info but we are using defaults.
    (data, units, calendar) = decode_cf_datetime(var)

    return {"data": data, "units": units, "calendar": calendar}


def nc4_ds_to_xr_dd(file_or_ds, xr_enc: dict = None) -> dict:
    """Convert a netCDF4 dataset to and xarray dataset dictionary"""

    if not isinstance(file_or_ds, nc4.Dataset):
        ds = nc4.Dataset(file_or_ds, "r")
    else:
        ds = file_or_ds

    # An empty xr_dd dictionary to hold the data
    xr_dd = deepcopy(template_xr_dd)

    # xr_dd["attrs"] = nc_file.__dict__  # ugly
    for attrname in ds.ncattrs():
        xr_dd["attrs"][attrname] = ds.getncattr(attrname)

    for dimname, dim in ds.dimensions.items():
        xr_dd["dims"][dimname] = len(dim)

    for varname, var in ds.variables.items():
        # _Encoding is used for string encoding in nc4
        var_encoding = {}
        if "_Encoding" in var.__dict__.keys():
            var_encoding["_Encoding"] = var._Encoding

        data_dict = {"dims": var.dimensions}

        var_attrs = {}
        for attrname in var.ncattrs():
            var_attrs[attrname] = var.getncattr(attrname)

        var_data, var_attrs, var_encoding = _nc4_var_to_datetime64(
            var,
            var_attrs,
            var_encoding,
        )

        if isinstance(var_data, np.ma.core.MaskedArray):
            var_data = var_data.data

        for aa in ["_FillValue"]:
            if aa in var_attrs.keys():
                var_encoding[aa] = var_attrs.pop(aa)

        data_dict["data"] = var_data
        data_dict["attrs"] = var_attrs
        data_dict["encoding"] = var_encoding

        if varname in ds.dimensions:
            xr_dd["coords"][varname] = data_dict
        else:
            xr_dd["data_vars"][varname] = data_dict

    ds.close()

    if xr_enc:
        xr_dd["encoding"] = {**xr_dd["encoding"], **xr_enc.pop("global")}
        for cc in xr_dd["coords"].keys():
            if cc in xr_enc.keys():
                xr_dd["coords"][cc]["encoding"] = {
                    **xr_dd["coords"][cc]["encoding"],
                    **xr_enc.pop(cc),
                }
        for vv in xr_dd["data_vars"].keys():
            if vv in xr_enc.keys():
                xr_dd["data_vars"][vv]["encoding"] = {
                    **xr_dd["data_vars"][vv]["encoding"],
                    **xr_enc.pop(vv),
                }

    return xr_dd


def _get_xr_encoding(nc_file) -> dict:
    ds = xr.open_dataset(nc_file)
    encoding = {}
    encoding["global"] = ds.encoding
    for vv in ds.variables:
        encoding[vv] = ds[vv].encoding

    ds.close()
    return encoding


def nc4_ds_to_dd(
    nc4_file_ds, subset: np.ndarray = None, use_xr_enc=True
) -> dict:
    """netCDF4 dataset to a pywatershed dataset dict."""
    xr_enc = None
    if not isinstance(nc4_file_ds, nc4.Dataset):
        if use_xr_enc:
            xr_enc = _get_xr_encoding(nc4_file_ds)
        nc4_file_ds = nc4.Dataset(nc4_file_ds)
    else:
        if use_xr_enc:
            raise ValueError(
                "Must pass a file and not an nc4.Dataset to use_xr_enc argument"
            )

    xr_dd = nc4_ds_to_xr_dd(nc4_file_ds, xr_enc=xr_enc)
    dd = xr_dd_to_dd(xr_dd)
    return dd


def dd_to_nc4_ds(dd, nc_file):
    """nc4_ds is a bit of a misnomer since it's on disk, dd_to_nc4"""
    # import cftime
    from xarray.coding.times import encode_cf_datetime

    dd = deepcopy(dd)

    # work from xrarray's dict representation of a dataset
    xr_dd = dd_to_xr_dd(dd)
    del dd

    # create a new netCDF4 file
    with nc4.Dataset(nc_file, "w") as ds:
        ds.set_fill_on()

        for key, value in xr_dd["attrs"].items():
            setattr(ds, key, value)

        for dim, size in xr_dd["dims"].items():
            ds.createDimension(dim, size)

        for coord_name, values in xr_dd["coords"].items():
            enc = values["encoding"]
            # handle time encoding
            if np.issubdtype(values["data"].dtype, np.datetime64):
                dates_enc, units, calendar = encode_cf_datetime(
                    values["data"].astype("datetime64[ns]"),
                    enc.get("units", None),
                    enc.get("calendar", None),
                )

                # TODO: what if .astype("datetime64[ns]") fails?
                # this would be a manual way that's looped and that
                # dosent guess at units and calendar
                # cf_dates = cftime.date2num(
                #     dates,
                #     units=enc.get("units"),
                #     calendar=enc.get("calendar"),
                # )

                values["data"] = dates_enc
                values["attrs"]["units"] = units
                values["attrs"]["calendar"] = calendar
                for kk in ["units", "calendar"]:
                    if kk in enc.keys():
                        del enc[kk]

            var_type = values["attrs"].get("type", values["data"].dtype)
            var = ds.createVariable(
                coord_name,
                var_type,
                dimensions=xr_dd["coords"][coord_name]["dims"],
                fill_value=enc.get("_FillValue", None),
                # This is not complete. Defaults from nc4
                zlib=enc.get("zlib", False),
                complevel=enc.get("complevel", 4),
                shuffle=enc.get("shuffle", True),
                contiguous=enc.get("contiguous", False),
                fletcher32=enc.get("fletcher32", False),
                chunksizes=enc.get("chunksizes", None),
            )
            var[:] = values["data"]
            var.setncatts(values["attrs"])
            var.coordinates = coord_name

        for var_name, values in xr_dd["data_vars"].items():
            enc = values["encoding"]
            var_type = values["attrs"].get("type", values["data"].dtype)
            var = ds.createVariable(
                var_name,
                var_type,
                dimensions=values["dims"],
                fill_value=enc.get("_FillValue", None),
                # This is not complete. Defaults from nc4
                zlib=enc.get("zlib", False),
                complevel=enc.get("complevel", 4),
                shuffle=enc.get("shuffle", True),
                contiguous=enc.get("contiguous", False),
                fletcher32=enc.get("fletcher32", False),
                chunksizes=enc.get("chunksizes", None),
            )

            var[:] = values["data"]
            var.setncatts(values["attrs"])

    return
