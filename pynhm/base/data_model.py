from copy import deepcopy
import itertools
from pprint import pprint

import cftime
import netCDF4 as nc4
import numpy as np
import xarray as xr

from .accessor import Accessor
from ..constants import listish, fill_values_dict, np_type_to_netcdf_type_dict

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
        validate: bool = True,
    ) -> "DatasetDict":
        self._data_vars = data_vars
        self._dims = dims
        self._coords = coords
        self._metadata = metadata
        self._encoding = encoding

        if validate:
            self.validate()

        return

    # TODO: add copy to properties? assuming that's necessary

    @property
    def dims(self, copy=False) -> dict:
        """Return the dimensions"""
        if copy:
            return deepcopy(self._dims)
        return self._dims

    @property
    def coords(self, copy=False) -> dict:
        """Return the coordinates"""
        if copy:
            return deepcopy(self._coords)
        return self._coords

    @property
    def data_vars(self, copy=False) -> dict:
        """Return the data_vars."""
        if copy:
            return deepcopy(self._data_vars)
        return self._data_vars

    @property
    def variables(self, copy=False) -> dict:
        """Return coords and data_vars together"""
        vars = {**self._coords, **self._data_vars}
        if copy:
            return deepcopy(vars)
        return vars

    @property
    def metadata(self, copy=False) -> dict:
        """Return the metadata"""
        if copy:
            return deepcopy(self._metadata)
        return self._metadata

    @property
    def encoding(self, copy=False) -> dict:
        """Return the encoding"""
        if copy:
            return deepcopy(self._encoding)
        return self._encoding

    @property
    def data(self, copy=False) -> dict:
        data = {
            "dims": self._dims,
            "coords": self._coords,
            "data_vars": self._data_vars,
            "metadata": self._metadata,
            "encoding": self._encoding,
        }
        if copy:
            return deepcopy(data)
        return data

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
    def from_dict(dict_in, copy=False):
        if copy:
            return DatasetDict(**deepcopy(dict_in))
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
    def from_netcdf(nc_file, use_xr=False) -> "DatasetDict":
        """Load from a netcdf file"""
        # handle more than one file?
        if use_xr:
            return DatasetDict(**xr_ds_to_dd(nc_file))
        else:
            return DatasetDict(**nc4_ds_to_dd(nc_file))

    def to_xr_ds(self) -> xr.Dataset:
        return dd_to_xr_ds(self.data)

    def to_xr_dd(self) -> dict:
        return dd_to_xr_dd(self.data)

    def to_nc4_ds(self, filename) -> None:
        return dd_to_nc4_ds(self.data, filename)

    def to_netcdf(self, filename, use_xr=False) -> None:
        """Write parameters to a netcdf file"""
        if use_xr:
            self.to_xr_ds().to_netcdf(filename)
        else:
            self.to_nc4_ds(filename)
        return

    def _get_var_dims(self, var_name, data=False, copy=False):
        dim_names = self._metadata[var_name]["dims"]
        if not data:
            return dim_names  # tuple, no need to copy
        else:
            dim_data = {
                kk: vv for kk, vv in self._dims.items() if kk in dim_names
            }
            if copy:
                return deepcopy(dim_data)
            return dim_data

    def _get_dim_coords(self, dim_list, data=False, copy=False):
        """Given a set of dimensions, get the corresponding coords"""
        # if all of a coords dims are in supplied dims, take the coordinate
        coords_out = []
        coord_dims = {cc: self._get_var_dims(cc) for cc in self._coords.keys()}
        for c_name, c_dims in coord_dims.items():
            # if c_dims is empty (a scalar coord), it automatically goes
            if len(set(c_dims) - set(dim_list)) == 0:
                coords_out += [c_name]
        if not data:
            return coords_out

        coord_data = {
            ck: cv for ck, cv in self._coords.items() if ck in coords_out
        }
        if copy:
            return deepcopy(coord_data)
        return coord_data

    def subset(
        self,
        keys: listish = None,
        # process: str = None,
        copy=True,
        keep_global_metadata=False,
        keep_global_encoding=False,
    ) -> "DatasetDict":
        """Subset a DatasetDict to keys in data_vars or coordinates."""
        # Instantiate the DatasetDict at end as deepcopy will be used
        # on the constructed subset dict (if requested)
        subset = deepcopy(template_dd)

        subset["metadata"]["global"] = {}
        subset["encoding"]["global"] = {}
        if keep_global_metadata:
            subset["metadata"]["global"] = self.metadata["global"]
        if keep_global_encoding:
            subset["encoding"]["global"] = self.encoding["global"]

        for vv in self.variables.keys():
            if vv not in keys:
                continue

            is_coord = vv in self._coords.keys()
            if is_coord:
                subset["coords"][vv] = self._coords[vv]
            else:
                subset["data_vars"][vv] = self._data_vars[vv]

            subset["metadata"][vv] = self._metadata[vv]
            subset["encoding"][vv] = self._encoding[vv]
            var_dim_data = self._get_var_dims(vv, copy=True, data=True)
            # dims
            for dd in var_dim_data:
                if dd not in subset["dims"].keys():  # faster?
                    subset["dims"][dd] = var_dim_data[dd]

            if not is_coord:
                # build coords from variables
                var_coord_data = self._get_dim_coords(
                    list(var_dim_data.keys()), data=True
                )
                for ck, cv in var_coord_data.items():
                    if ck not in subset["coords"].keys():
                        subset["coords"][ck] = cv

        # build metadata and encoding from coords
        for cc in subset["coords"].keys():
            for aa in ["metadata", "encoding"]:
                if cc in subset[aa].keys():
                    continue
                subset[aa][cc] = self[aa][cc]

        result = DatasetDict.from_dict(subset, copy=copy)
        return result

    def get_parameters(
        self, keys: listish, process: str = None, dims: bool = False
    ) -> dict:
        """Returns Parameter data for requested keys."""
        # Provide in base class
        raise NotImplementedError

    def validate(self):
        # required keys
        assert sorted(self.data.keys()) == sorted(
            ["dims", "coords", "data_vars", "metadata", "encoding"]
        )

        # can not have same names in coords and data_vars
        common_keys = set(self.coords.keys()).intersection(
            set(self.data_vars.keys())
        )
        assert len(common_keys) == 0

        # metadata and encoding keys against variable keys
        meta_keys = set(self.metadata.keys())
        enc_keys = set(self.encoding.keys())
        var_keys = set(self.variables.keys())
        assert meta_keys == enc_keys
        assert meta_keys == var_keys.union(set(["global"]))

        # check all vars dims exist
        var_dims = [self.metadata[kk]["dims"] for kk in self.variables.keys()]
        dims = set([dim for dims in var_dims for dim in dims])
        for dd in dims:
            assert dd in self.dims.keys()

        return

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

    if isinstance(time_data, cftime.real_datetime):
        time_data = np.datetime64(time_data)
    else:
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
    from xarray.coding.times import decode_cf_datetime

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

    # have to promote coord variables here?
    all_coords = [
        vv["attrs"].pop("coordinates")
        for vv in xr_dd["data_vars"].values()
        if "coordinates" in vv["attrs"].keys()
    ]
    all_coords = sorted(set(" ".join(all_coords).split(" ")))
    for cc in all_coords:
        if cc == "":
            continue
        xr_dd["coords"][cc] = xr_dd["data_vars"].pop(cc)

    # handle bools
    for meta in [xr_dd["coords"], xr_dd["data_vars"]]:
        for kk, vv in meta.items():
            if "dtype" in vv["attrs"].keys():
                dtype = vv["attrs"]["dtype"]
                if dtype == "bool":
                    vv["data"] = vv["data"].astype("bool")
                    _ = vv["attrs"].pop("dtype")

    # bring in the encoding information using xarray (cheating?)
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

        for var_name, values in xr_dd["data_vars"].items():
            enc = values["encoding"]
            is_bool = values["data"].dtype == "bool"
            np_type = values["data"].dtype
            nc_type = np_type_to_netcdf_type_dict[np_type]
            var_type = values["attrs"].get("type", nc_type)
            if is_bool:
                var_type = "i1"
            default_fill = fill_values_dict[np_type]
            var = ds.createVariable(
                var_name,
                var_type,
                dimensions=values["dims"],
                fill_value=enc.get("_FillValue", default_fill),
                # This is not complete. Defaults from nc4
                zlib=enc.get("zlib", False),
                complevel=enc.get("complevel", 4),
                shuffle=enc.get("shuffle", True),
                contiguous=enc.get("contiguous", False),
                fletcher32=enc.get("fletcher32", False),
                chunksizes=enc.get("chunksizes", None),
            )

            if is_bool:
                var[:] = values["data"].astype("int8")
            else:
                var[:] = values["data"]
            var.setncatts(values["attrs"])
            # solve coords for this var: any coord which has any of its dims
            coords = []
            # these rules are a bit opaque
            coords_just_cuz = ["reference_time"]
            coords_no_way = ["time"]
            for c_name, c_val in xr_dd["coords"].items():
                if c_name in coords_no_way:
                    continue
                common = set(values["dims"]).intersection(set(c_val["dims"]))
                if len(common) or c_name in coords_just_cuz:
                    coords += [c_name]
            if len(coords):
                var.coordinates = " ".join(sorted(coords))
            if is_bool:
                var.setncattr("dtype", "bool")

    return
