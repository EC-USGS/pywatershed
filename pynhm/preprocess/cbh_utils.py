import math
import pathlib as pl
from datetime import datetime
from typing import Union

import netCDF4 as nc4
import numpy as np
import pandas as pd

from ..utils.parameters import PrmsParameters
from .cbh_metadata import cbh_metadata

zero = np.zeros((1))[0]
one = np.ones((1))[0]

# Compound types
file_type = Union[str, pl.Path]
fileish = Union[str, pl.Path, dict]


cbh_units = {
    key: val["units"]
    for key, val in cbh_metadata.items()
    if "units" in list(val.keys())
}

created_line = "Created "
written_line = "Written "
hash_line = "##############"
hash_line_official = "########################################"


def _cbh_file_to_df(
    the_file: file_type, params: PrmsParameters = None
) -> pd.DataFrame:
    # Only take cbh files that contain single variables

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
    for posn in range(len(meta_lines)):
        key, count = meta_lines[posn].split(" ")
        count = int(count)
        zs = math.ceil(math.log(count, 10))
        if params is None:
            col_names += [f"{key}{str(ii).zfill(zs)}" for ii in range(count)]
        else:
            col_names = np.char.add(
                np.array([key]), params._parameter_data["nhm_id"].astype(str)
            ).tolist()
        var_count_dict[key] = count

    if len(var_count_dict) > 1:
        msg = f"cbh input files should contain only one variable each: {the_file}"
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
    msg = f"Number of actual data columns does not match metadata info: {meta_lines}"
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
    # JLM TODO: Set datetime resolution to hours? or mins. Could do days but might look forward a bit.
    data = data.drop(columns=set(date_cols))
    data = data.set_index("date")

    return data


def _cbh_files_to_df(
    file_dict: dict, params: PrmsParameters = None
) -> pd.DataFrame:
    dfs = [_cbh_file_to_df(val, params) for val in file_dict.values()]
    return pd.concat(dfs, axis=1)


def cbh_files_to_df(
    files: fileish, params: PrmsParameters = None
) -> pd.DataFrame:
    if isinstance(files, (str, pl.Path)):
        df = _cbh_file_to_df(files, params)
    elif isinstance(files, (dict)):
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
    col_name_tuples = [_col_name_split(kk) for kk, vv in df.iteritems()]
    df.columns = pd.MultiIndex.from_tuples(col_name_tuples)
    var_names = df.columns.unique(level=0)
    np_dict = {}
    np_dict["datetime"] = df.index.to_numpy(copy=True).astype("datetime64[s]")
    spatial_ids = df.loc[:, var_names[0]].columns.values
    if spatial_ids[0] == "000":
        np_dict["hru_ind"] = spatial_ids
    else:
        np_dict["nhm_id"] = spatial_ids
    for vv in var_names:
        np_dict[vv] = df[vv].to_numpy(copy=True)
    return np_dict


def cbh_n_hru(np_dict: dict) -> int:
    odd_shapes = ["datetime", "hru_ind", "nhm_id"]
    shapes = [
        var.shape for key, var in np_dict.items() if key not in odd_shapes
    ]
    for ss in shapes:
        assert shapes[0] == ss
    return shapes[0][1]


def cbh_n_time(np_dict: dict) -> int:
    return np_dict["datetime"].shape[0]


def cbh_files_to_np_dict(files: fileish, params: PrmsParameters) -> dict:
    np_dict = cbh_df_to_np_dict(cbh_files_to_df(files, params))
    return np_dict


def cbh_adjust(cbh_dict: dict, params: PrmsParameters) -> dict:
    # Param object has no defined interface at this time.
    if params is None:
        raise ValueError("Parameters have not been supplied for adjustment.")
    param_data = params._parameter_data
    nhru = params._dimensions["nhru"]

    # I dislike using pd for something that seems like it should exist in np
    month_ind_12 = pd.to_datetime(cbh_dict["datetime"]).month - 1  # (time)
    month_ind_1 = np.zeros(cbh_dict["datetime"].shape, dtype=int)  # (time)

    # adjust temp ---------------------------------
    # Temperature bias corrections are calibrated by PRMS/HNM
    # https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/physics/sm_temperature_hru.f90#L91
    tmax_param = param_data["tmax_cbh_adj"]  # (12 or 1, space)
    tmin_param = param_data["tmin_cbh_adj"]

    # i suppose thse could have different shapes.
    # throw an error if that happens
    if tmax_param.shape != tmin_param.shape:
        msg = "Not implemented: tmin/tmax cbh adj parameters with different shapes"
        raise NotImplementedError(msg)

    if tmax_param.shape[0] == 12:
        month_ind = month_ind_12
    elif tmax_param.shape[0] == 1:
        month_ind = month_ind_1
    else:
        msg = (
            "Unexpected month dimension for cbh temperature adjustment params"
        )
        raise ValueError(msg)

    cbh_dict["tmax_adj"] = np.zeros(
        cbh_dict["tmax"].shape, dtype=cbh_dict["tmax"].dtype
    )  # (time, space)
    cbh_dict["tmin_adj"] = np.zeros(
        cbh_dict["tmin"].shape, dtype=cbh_dict["tmin"].dtype
    )
    cbh_dict["tmax_adj"] = cbh_dict["tmax"] + tmax_param[month_ind]
    cbh_dict["tmin_adj"] = cbh_dict["tmin"] + tmin_param[month_ind]

    # To check the resulting shape used above
    # tmax_time_params = tmax_param[month_ind]
    # for tt in range(len(month_ind)):
    #     assert (tmax_time_params[tt,:] == tmax_param[month_ind[tt]]).all()

    # adjust precip/snowfall ---------------------------------
    # Snow/rain partitioning of total precip is (adjusted) temperature dependent in addition
    # to being dependent on the following parameters (not sure if all are calibrated).
    # https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/physics/sm_precipitation.f90#L227
    tmax_allsnow_param = param_data["tmax_allsnow"]
    tmax_allrain_offset_param = param_data["tmax_allrain_offset"]
    snow_cbh_adj_param = param_data["snow_cbh_adj"]
    rain_cbh_adj_param = param_data["rain_cbh_adj"]
    adjmix_rain_param = param_data["adjmix_rain"]

    # i suppose thse could have different shapes.
    # throw an error if that happens
    shape_list = np.array(
        [
            tmax_allsnow_param.shape[0],
            tmax_allrain_offset_param.shape[0],
            snow_cbh_adj_param.shape[0],
            rain_cbh_adj_param.shape[0],
            adjmix_rain_param.shape[0],
        ]
    )
    if not (shape_list == 12).all():
        msg = "Not implemented: tmin/tmax cbh adj parameters with different shapes"
        raise NotImplementedError(msg)

    tmax_allrain_param = tmax_allsnow_param + tmax_allrain_offset_param

    if tmax_allsnow_param.shape[0] == 12:
        month_ind = month_ind_12
    elif tmax_allsnow_param.shape[0] == 1:
        month_ind = month_ind_1
    else:
        msg = "Unexpected month dimension for cbh precip adjustment params"
        raise ValueError(msg)

    cbh_dict["prcp_adj"] = np.zeros(
        cbh_dict["prcp"].shape, dtype=cbh_dict["prcp"].dtype
    )  # (time, space)
    cbh_dict["rainfall_adj"] = np.zeros(
        cbh_dict["prcp"].shape, dtype=cbh_dict["prcp"].dtype
    )
    cbh_dict["snowfall_adj"] = np.zeros(
        cbh_dict["prcp"].shape, dtype=cbh_dict["prcp"].dtype
    )

    # JLM: in my test vectorization was 45x faster for drb_2yr: 3.728s:.083s for loop:vectorized
    prmx = np.zeros(cbh_dict["prcp"].shape, dtype=cbh_dict["prcp"].dtype)

    # Order MATTERS in calculating the prmx mask
    # The logic in PRMS is if(all_snow),elif(all_rain),else(mixed)
    # so we set the mask in the reverse order

    # Calculate the mix everywhere, then set the precip/rain/snow amounts from the conditions.
    tdiff = cbh_dict["tmax_adj"] - cbh_dict["tmin_adj"]
    prmx = (
        (cbh_dict["tmax_adj"] - tmax_allsnow_param[month_ind]) / tdiff
    ) * adjmix_rain_param[month_ind]
    del tdiff

    wh_all_snow = np.where(
        cbh_dict["tmax_adj"] <= tmax_allsnow_param[month_ind]
    )
    wh_all_rain = np.where(
        np.logical_or(
            cbh_dict["tmin_adj"] > tmax_allsnow_param[month_ind],
            cbh_dict["tmax_adj"] >= tmax_allrain_param[month_ind],
        )
    )
    prmx[wh_all_rain] = one
    prmx[wh_all_snow] = zero

    # Recalculate/redefine these now based on prmx instead of the temperature logic
    wh_all_snow = np.where(prmx <= zero)
    wh_all_rain = np.where(prmx >= one)

    # Mixed case (everywhere, to be overwritten by the all-snow/rain-fall cases)
    cbh_dict["prcp_adj"] = cbh_dict["prcp"] * snow_cbh_adj_param[month_ind]
    cbh_dict["rainfall_adj"] = prmx * cbh_dict["prcp_adj"]
    cbh_dict["snowfall_adj"] = cbh_dict["prcp_adj"] - cbh_dict["rainfall_adj"]
    del prmx

    # All precip is snow case
    # The condition to be used later:
    cbh_dict["prcp_adj"][wh_all_snow] = (
        cbh_dict["prcp"] * snow_cbh_adj_param[month_ind]
    )[wh_all_snow]
    cbh_dict["snowfall_adj"][wh_all_snow] = cbh_dict["prcp_adj"][wh_all_snow]
    cbh_dict["rainfall_adj"][wh_all_snow] = zero

    # All precip is rain case
    # The condition to be used later:
    cbh_dict["prcp_adj"][wh_all_rain] = (
        cbh_dict["prcp"] * rain_cbh_adj_param[month_ind]
    )[wh_all_rain]
    cbh_dict["rainfall_adj"][wh_all_rain] = cbh_dict["prcp_adj"][wh_all_rain]
    cbh_dict["snowfall_adj"][wh_all_rain] = zero

    return None


def cbh_convert_units():
    pass


def cbh_check(cbh_dict: dict, verbosity: int = 0) -> None:

    assert (~np.isnat(cbh_dict["datetime"])).all()

    for adj in ["", "_adj"]:

        # assume one variable represents if there are any adjustments
        if f"tmax{adj}" not in cbh_dict.keys():
            continue

        tmaxvar = f"tmax{adj}"
        tminvar = f"tmin{adj}"
        if (tminvar in cbh_dict.keys()) and (tmaxvar in cbh_dict.keys()):
            if not (cbh_dict[tmaxvar] >= cbh_dict[tminvar]).all():
                msg = f"{tmaxvar} < {tminvar}"
                raise ValueError(msg)

        # assert (cbh_dict['tmax'] >= zero).all()  # use units for checking max/minimums?

        prcpvar = f"prcp{adj}"
        if prcpvar in cbh_dict.keys():
            if not (cbh_dict[prcpvar] >= zero).all():
                msg = f"{prcpvar} contains negative values"
                raise ValueError(msg)

    if verbosity >= 1:
        print("cbh state check successful")
    return


def cbh_to_netcdf(
    np_dict: dict,
    filename: fileish,
    clobber: bool = True,
    output_vars: list = None,
    zlib: bool = True,
    complevel: int = 4,
    global_atts: dict = {},
    chunk_sizes={"time": 30, "hru": 0},
) -> None:
    # Default time chunk is for a read pattern of ~monthly at a time.

    ds = nc4.Dataset(filename, "w", clobber=clobber)
    ds.setncattr("Description", "Climate by HRU")
    for key, val in global_atts.items():
        ds.setncattr(key, val)

    # Dimensions
    # None for the len argument gives an unlimited dim
    ds.createDimension("time", None)  # cbh_n_time(np_dict))
    ds.createDimension("hru", cbh_n_hru(np_dict))

    # Dim Variables
    time = ds.createVariable(
        "datetime", cbh_metadata["datetime"]["type"], ("time")
    )
    for att, val in cbh_metadata["datetime"].items():
        time.setncattr(att, val)
    time[:] = nc4.date2num(
        np_dict["datetime"].astype(datetime),
        units=cbh_metadata["datetime"]["units"],
        calendar=cbh_metadata["datetime"]["calendar"],
    )

    hru_name = "hru_ind" if "hru_ind" in np_dict.keys() else "nhm_id"
    hruid = ds.createVariable(
        hru_name, cbh_metadata[hru_name]["type"], ("hru")
    )
    for att, val in cbh_metadata[hru_name].items():
        hruid.setncattr(att, val)
    hruid[:] = np_dict[hru_name]

    # Variables
    var_list = set(np_dict.keys()).difference(
        {"datetime", "hru_ind", "nhm_id"}
    )
    if output_vars is not None:
        var_list = [var for var in var_list if var in output_vars]

    for vv in var_list:
        vvtype = cbh_metadata[vv]["type"]
        var = ds.createVariable(
            vv,
            vvtype,
            ("time", "hru"),
            fill_value=nc4.default_fillvals[vvtype],  # JLM: sus
            zlib=zlib,
            complevel=complevel,
            chunksizes=tuple(chunk_sizes.values()),
        )
        for att, val in cbh_metadata[vv].items():
            if att in ["_FillValue"]:
                continue
            var.setncattr(att, val)
        ds.variables[vv][:, :] = np_dict[vv]

    ds.close()
    print(f"Wrote netcdf file: {filename}")
    return
