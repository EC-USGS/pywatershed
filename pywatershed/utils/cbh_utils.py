import math
import pathlib as pl
from datetime import datetime
from typing import Union

import netCDF4 as nc4
import numpy as np
import pandas as pd

from ..base import meta
from ..parameters import PrmsParameters

zero = np.zeros((1))[0]
one = np.ones((1))[0]

# Compound types

nc4_to_np_types = {
    "i4": "int32",
    "i8": "int64",
    "f4": "float32",
    "f8": "float64",
}
np_to_nc4_types = {val: key for key, val in nc4_to_np_types.items()}

created_line = "Created "
written_line = "Written "
hash_line = "##############"
hash_line_official = "########################################"


def _cbh_file_to_df(
    the_file: Union[str, pl.Path], params: PrmsParameters = None
) -> pd.DataFrame:
    # Only take cbh files that contain single variables
    # JLM: this can be substantially simplified as
    # we are only reading NHM CBH files now.
    # Revisit with hru1.yaml

    meta_lines = []
    with open(the_file, "r") as file_open:
        wh_hash_line = -1
        the_line = ""
        while hash_line not in the_line:
            wh_hash_line += 1
            the_line = file_open.readline()
            if (
                ("//" not in the_line)
                and (created_line not in the_line)
                and (written_line not in the_line)
                and (hash_line not in the_line)
            ):
                meta_lines += [the_line.strip()]

    col_names = []
    var_count_dict = {}
    line = meta_lines[-1]
    line_split = line.split(" ")
    key = line_split[0]
    count = int(line_split[-1])
    zs = math.ceil(math.log(count, 10))
    if (params is None) or ("nhm_id" not in params.parameters):
        col_names += [f"{key}{str(ii).zfill(zs)}" for ii in range(count)]
    else:
        col_names = np.char.add(
            np.array([key]), params.parameters["nhm_id"].astype(str)
        ).tolist()
    var_count_dict[key] = count

    if len(var_count_dict) > 1:
        msg = (
            "cbh input files should contain only one variable each: "
            f"{the_file}"
        )
        raise ValueError(msg)

    dtypes = (["str"] * 6) + (["float64"] * len(col_names))
    date_cols = ["Y", "m", "d", "H", "M", "S"]
    col_names = date_cols + col_names
    assert len(dtypes) == len(col_names)
    dtype_dict = dict(zip(col_names, dtypes))
    assert len(dtype_dict) == len(col_names)

    # Some files erroneously specify the number of columns
    # catch that here
    data = pd.read_csv(
        the_file,
        nrows=1,
        # names=col_names,  # dont specify names
        header=None,  # dont use default header from first line
        index_col=False,
        skiprows=wh_hash_line + 1,
        delim_whitespace=True,
        dtype=dtype_dict,
    )
    msg = (
        "Number of actual data columns does not match metadata info: "
        f"{meta_lines}"
    )
    assert len(data.columns) == len(col_names), msg
    # JLM: is the above sufficient?

    data = pd.read_csv(
        the_file,
        names=col_names,
        index_col=False,
        skiprows=wh_hash_line + 1,
        delim_whitespace=True,
        dtype=dtype_dict,
    )

    data["date"] = pd.to_datetime(
        data.Y + "-" + data.m.str.zfill(2) + "-" + data.d.str.zfill(2)
    )
    # JLM TODO: Set datetime resolution to hours? or mins. Could do days but
    # might look forward a bit.
    data = data.drop(columns=set(date_cols))
    data = data.set_index("date")

    return data


def _cbh_files_to_df(
    file_dict: dict, params: PrmsParameters = None
) -> pd.DataFrame:
    if isinstance(file_dict, dict):
        dfs = [_cbh_file_to_df(val, params) for val in file_dict.values()]
    else:
        dfs = [_cbh_file_to_df(val, params) for val in file_dict]
    return pd.concat(dfs, axis=1)


def cbh_files_to_df(
    files: Union[str, pl.Path], params: PrmsParameters = None
) -> pd.DataFrame:
    if isinstance(files, (str, pl.Path)):
        df = _cbh_file_to_df(files, params)
    elif isinstance(files, (dict, list)):
        df = _cbh_files_to_df(files, params)
    else:
        raise ValueError(
            f'"files" argument of type {type(files)} not accepted.'
        )
    return df


def _col_name_split(string: str) -> tuple:
    char_list = list(string)
    wh_digit = [ii for ii in range(len(char_list)) if char_list[ii].isdigit()]
    return string[: wh_digit[0]], string[wh_digit[0] :]


def cbh_df_to_np_dict(df: pd.DataFrame) -> dict:
    # Convert to a multi-index? For getting to a dict of np arrays
    # This might go elsewhere
    col_name_tuples = [_col_name_split(kk) for kk, vv in df.items()]
    df.columns = pd.MultiIndex.from_tuples(col_name_tuples)
    var_names = df.columns.unique(level=0)
    np_dict = {}
    np_dict["time"] = df.index.to_numpy(copy=True).astype("datetime64[s]")
    spatial_ids = (
        df.loc[:, var_names[0]].columns.values.astype(float).astype(int)
    )
    if spatial_ids[0] == "000":
        np_dict["hru_ind"] = spatial_ids
    else:
        np_dict["nhm_id"] = spatial_ids
    for vv in var_names:
        np_dict[vv] = df[vv].to_numpy(copy=True)
    return np_dict


def cbh_n_hru(np_dict: dict) -> int:
    odd_shapes = ["time", "hru_ind", "nhm_id"]
    shapes = [
        var.shape for key, var in np_dict.items() if key not in odd_shapes
    ]
    for ss in shapes:
        assert shapes[0] == ss
    return shapes[0][1]


def cbh_n_time(np_dict: dict) -> int:
    return np_dict["time"].shape[0]


def cbh_files_to_np_dict(
    files: Union[str, pl.Path], params: PrmsParameters
) -> dict:
    np_dict = cbh_df_to_np_dict(cbh_files_to_df(files, params))
    return np_dict


def cbh_file_to_netcdf(
    input_file: Union[str, pl.Path],
    parameters: PrmsParameters,
    nc_file: Union[str, pl.Path],
    clobber: bool = True,
    output_vars: list = None,
    zlib: bool = True,
    complevel: int = 4,
    global_atts: dict = None,
    rename_vars: dict = None,
    chunk_sizes: dict = None,
    time_units="days since 1979-01-01 00:00:00",
    time_calendar="standard",
) -> None:
    """Convert PRMS native CBH files to NetCDF format for pywatershed

    Args:
        input_file: the CBH file to read
        parameters: the Parameters object of PRMS parameters for this domain
        nc_file: the NetCDF output file
        clobber: Overwrite an existing NetCDF file?
        output_vars: Subset of variables in the CBH to write?
        zlib: use zlib compression?
        complevel: The level of compression.
        global_atts: A dictionary of meta data for the variables.
        rename_vars: a dictionary for mapping CBH variable names
        chunk_sizes: along each dimension, default {"time": 30, "hru": 0},
        time_units: default "days since 1979-01-01 00:00:00",
        time_calendar: default"standard",
    """

    if rename_vars is None:
        rename_vars = {}
    if global_atts is None:
        global_atts = {}
    if chunk_sizes is None:
        chunk_sizes = {"time": 30, "hru": 0}

    np_dict = cbh_files_to_np_dict(input_file, parameters)

    # Default time chunk is for a read pattern of ~monthly at a time.

    ds = nc4.Dataset(nc_file, "w", clobber=clobber)
    ds.setncattr("Description", "Climate by HRU")
    for key, val in global_atts.items():
        ds.setncattr(key, val)

    # Dimensions
    # None for the len argument gives an unlimited dim
    ds.createDimension("time", None)  # cbh_n_time(np_dict))
    ds.createDimension("hru", cbh_n_hru(np_dict))

    # Dim Variables
    time = ds.createVariable("time", "f4", ("time",))

    time[:] = nc4.date2num(np_dict["time"].astype(datetime), time_units)
    time.units = time_units
    # time.calendar = time_calendar

    hru_name = "hru_ind" if "hru_ind" in np_dict.keys() else "nhm_id"
    hruid = ds.createVariable(
        hru_name, meta.get_types(hru_name)[hru_name], ("hru")
    )
    hru_meta_dict = {
        "hru_ind": {
            "type": "i4",
            "desc": "Hydrologic Response Unit (HRU) index",
            "cf_role": "timeseries_id",
        },
        "nhm_id": {
            "type": "i4",
            "desc": "NHM Hydrologic Response Unit (HRU) ID",
            "cf_role": "timeseries_id",
        },
    }

    for att, val in hru_meta_dict[hru_name].items():
        hruid.setncattr(att, val)
    hruid[:] = np_dict[hru_name]

    # Variables
    var_list = set(np_dict.keys()).difference({"time", "hru_ind", "nhm_id"})
    if output_vars is not None:
        var_list = [var for var in var_list if var in output_vars]

    for vv in var_list:
        vv_meta = meta.get_vars(vv)[vv]
        vv_type = vv_meta["type"]
        var_name_out = vv
        if vv in rename_vars.keys():
            var_name_out = rename_vars[vv]
        var = ds.createVariable(
            var_name_out,
            vv_type,
            ("time", "hru"),
            fill_value=nc4.default_fillvals[np_to_nc4_types[vv_type]],
            zlib=zlib,
            complevel=complevel,
            chunksizes=tuple(chunk_sizes.values()),
        )
        for att, val in vv_meta.items():
            if att in ["_FillValue", "type", "dimensions"]:
                continue
            var.setncattr(att, val)
        ds.variables[var_name_out][:, :] = np_dict[vv]

    ds.close()
    print(f"Wrote netcdf file: {nc_file}")
    return
