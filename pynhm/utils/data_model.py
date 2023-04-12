from copy import deepcopy

import netCDF4 as nc4
import numpy as np
import xarray as xr

# This file defines the data model for pywatershed. It is called a
# "dataset_dict" and has a invertible mapping with non-hierarchical netcdf
# or xarray datasets.

# TODO: dataset_dict class?
#       why?: methods for consistency?

# TODO: what about hierarchical/groups in netcdf files?

# TODO: nc4_ds_to_dd inverse

# TODO: def validate dataset_dict
#         check for all required keys
#         check metadata keys against data
#         check all dims, coords

# This is what a dataset_dict looks like. Metadata for coord and data_vars
# is found in metadata.
template_dataset_dict = {
    "dims": {},
    "coords": {},
    "data_vars": {},
    "metadata": {},
    "encoding": {},
}

# Show a basic example of how an xr_dict is changed to a dd.


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

    dd = xr_dict_to_dd(dd)

    return dd


def xr_dict_to_dd(xr_dict: dict) -> dict:
    dd = deepcopy(xr_dict)

    # Move the global encoding to a global key of itself
    dd["encoding"] = {"global": dd["encoding"]}

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


def dd_to_xr_dict(dd: dict) -> dict:
    dd = deepcopy(dd)

    # remove metadata and encoding from the dict/model
    meta = dd.pop("metadata")
    encoding = dd.pop("encoding")

    # loop over meta data, putting it back on the data_vars and/or coords
    # and moving the data to "data"
    dd["attrs"] = meta.pop("global")
    dd["encoding"] = encoding.pop("global")

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
    return xr.Dataset.from_dict(dd_to_xr_dict(dd))


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
    from xarray.coding.times import encode_cf_datetime

    # This can take encoding info but we are using defaults.
    (data, units, calendar) = encode_cf_datetime(var)

    return {"data": data, "units": units, "calendar": calendar}


def nc4_ds_to_dd(nc_file, subset: np.ndarray = None):
    """netCDF4 dataset to a pywatershed dataset dict."""
    nc = nc4.Dataset(nc_file)

    # Create an empty dictionary to hold the data
    dataset = {
        "attrs": {},
        "dims": {},
        "coords": {},
        "data_vars": {},
        "encoding": {},
    }

    # dataset["attrs"] = nc_file.__dict__  # ugly
    for attrname in nc.ncattrs():
        dataset["attrs"][attrname] = nc.getncattr(attrname)

    for dimname, dim in nc.dimensions.items():
        dataset["dims"][dimname] = len(dim)

    for varname, var in nc.variables.items():
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

        if varname in nc.dimensions:
            dataset["coords"][varname] = data_dict
        else:
            dataset["data_vars"][varname] = data_dict

    nc.close()

    dataset["encoding"] = {}

    # convert to pyws dataset_dict
    data_dict = xr_dict_to_dd(dataset)

    return data_dict


# def dd_to_nc4_ds(dd, nc_file):
#     import cftime
#     import datetime

#     dd = dd_to_xr_dict(deepcopy(dd))

#     # create a new netCDF4 file
#     with nc4.Dataset(nc_file, "w") as ds:
#         # add global attributes
#         for key, value in dd["attrs"].items():
#             setattr(ds, key, value)

#         # add dimensions
#         for dim, size in dd["dims"].items():
#             ds.createDimension(dim, size)

#         # add coordinates
#         for coord_name, values in dd["coords"].items():
#             var = ds.createVariable(
#                 coord_name, values["attrs"]["type"], dimensions=(coord_name,)
#             )

#             if not isinstance(var, np.datetime64):
#                 # dates = values["data"].astype(datetime.datetime).tolist()
#                 # cf_dates = cftime.num2date(dates, calendar=')

#                 encoding_dict = _datetime64_to_nc4_var(dates)
#                 asdf

#                 values["data"] = encoding_dict.pop("data")
#                 values["attrs"] = {**values["attrs"], **encoding_dict}

#             var[:] = values["data"]
#             var.setncatts(values["attrs"])

#         # add data variables
#         for var_name, values in dd["data_vars"].items():
#             var = ds.createVariable(
#                 var_name, values["attrs"]["type"], dimensions=values["dims"]
#             )
#             var[:] = values["data"]
#             var.setncatts(values["attrs"])
