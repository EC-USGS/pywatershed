import datetime
import functools
import os
from io import StringIO

import numpy as np
import pandas as pd

inch_to_meter = 0.0254
acre_to_meter_squared = 4046.8564224

conversions = {
    # forcings
    "precipitation": inch_to_meter,
    "temp_min": 1.0,
    "temp_max": 1.0,
    # stats
    "basin_potet": inch_to_meter,
    # wbal files
    "soilzone_last_sm": inch_to_meter,
    # parameters
    "intcp_stor_max": inch_to_meter,
    # output
    "net_ppt": inch_to_meter,
    "soil_moist_ante": inch_to_meter,
    "hru_sroffi": inch_to_meter,
    "hru_sroffp": inch_to_meter,
}


def unit_conversion(data, verbose=False):
    if isinstance(data, (dict,)):
        if verbose:
            print("dictionary conversion...")
        keys = list(data.keys())
    elif isinstance(data, (pd.DataFrame,)):
        if verbose:
            print("dataframe conversion...")
        keys = data.columns.to_list()
    elif isinstance(data, (np.ndarray,)):
        if verbose:
            print("numpy array conversion...")
        try:
            keys = list(data.dtype.names)
        except:
            keys = []
    else:
        keys = []

    for key in keys:
        if key in conversions.keys():
            if verbose:
                print(f"converting {key}...using {conversions[key]} multipler")
            data[key] *= conversions[key]

    return data


def load_prms_input(
    input_data_path, datanames, filenames, convert=True, verbose=False
):
    # load prms input
    templist = []
    for dataname, cbhname in zip(datanames, filenames):
        fpath = os.path.join(input_data_path, cbhname)
        print(f"Loading {fpath}")
        filelist = f"date,{dataname}\n"
        with open(fpath) as f:
            for i, line in enumerate(f):
                if i > 2:
                    ll = line.strip().split()
                    yr = int(ll[0])
                    mo = int(ll[1])
                    da = int(ll[2])
                    data = float(ll[-1])
                    # dt = datetime.datetime(yr, mo, da)
                    filelist += f"{da:02d}/{mo:02d}/{yr:04d}, {data}\n"
        tdf = pd.read_csv(
            StringIO(filelist),
            parse_dates=["date"],
            index_col=["date"],
            dtype=float,
            float_precision="high",
        )
        templist.append(tdf)

    # concatenate individual dataframes
    df = pd.concat([v for v in templist], axis=1)

    # unit conversion
    if convert:
        unit_conversion(df, verbose=verbose)

    return df


class PrmsParameters:
    """PRMS parameter class

    parameter_file: str
        path to PRMS parameter file
    """

    def __init__(self, parameter_file):
        (
            self._dimensions,
            self._parameter_data,
            self._parameter_dimensions,
            self._parameter_types,
        ) = _load_prms_parameters(parameter_file)

    @property
    def get_dimensions(self):
        return self._dimensions

    @property
    def get_parameter_data(self):
        return self._parameter_data

    @property
    def get_parameter_dimensions(self):
        return self._parameter_dimensions

    @property
    def get_parameter_types(self):
        return self._parameter_types


def _load_prms_parameters(parameter_file):
    """Read a PRMS parameter file

    :param parameter_file:
    :return:
    """
    line_num = 0
    vals = {}
    dims = {}
    param_dims = {}
    param_type = {}

    with open(parameter_file) as f:
        reading_dims = False
        for line in f:
            try:
                line = line.rstrip()  # remove '\n' at end of line
                line_num += 1
                if line == "** Dimensions **":
                    reading_dims = True
                    line = f.readline().rstrip()
                    line_num += 1

                if line == "** Parameters **":
                    reading_dims = False
                    break

                if reading_dims:
                    line = f.readline().rstrip()
                    line_num += 1
                    dim_name = line

                    line = f.readline().rstrip()
                    line_num += 1
                    size = line

                    if dim_name in dims.keys():
                        pass
                    else:
                        dims[dim_name] = int(size)
            except:
                msg = (
                    f"read parameters exception line = {line}\n"
                    + f"read parameters exception line_num = {str(line_num)}\n"
                )
                raise ValueError(msg)

        #        read params
        for line in f:
            try:
                line = line.rstrip()  # remove '\n' at end of line
                line_num += 1

                if line == "####":
                    line = f.readline().rstrip()
                    line = line.split(" ", 1)[
                        0
                    ]  # old format parameter files have a blank (' ') and then a width format value. Strip this off.
                    param_name = line
                    line_num += 1

                    line = f.readline().rstrip()
                    line_num += 1
                    num_dims = int(line)
                    pd = [None] * num_dims
                    for ii in range(num_dims):
                        line = f.readline().rstrip()
                        pd[ii] = line
                        line_num += 1

                    param_dims[param_name] = pd

                    line = f.readline().rstrip()
                    line_num += 1
                    num_vals = int(line)
                    line = f.readline().rstrip()
                    line_num += 1
                    tp = int(line)
                    param_type[param_name] = tp

                    if tp == 2:
                        vs = np.zeros(num_vals, dtype=float)
                        for jj in range(num_vals):
                            line = f.readline().rstrip()
                            line_num += 1
                            vs[jj] = float(line)

                    elif tp == 1:
                        vs = np.zeros(num_vals, dtype=int)
                        for jj in range(num_vals):
                            line = f.readline().rstrip()
                            line_num += 1
                            vs[jj] = int(line)

                    else:
                        vs = np.zeros(num_vals, dtype=np.chararray)
                        for jj in range(num_vals):
                            line = f.readline().rstrip()
                            line_num += 1
                            vs[jj] = line

                    if num_dims == 2:
                        vs.shape = (dims[pd[1]], dims[pd[0]])

                    if param_name in vals.keys():
                        print(
                            "parameter ",
                            param_name,
                            " is already in ",
                            parameter_file,
                        )
                    else:
                        vals[param_name] = vs

            except:
                raise ValueError(
                    f"read parameters exception line_num = {line_num}"
                )

    return dims, vals, param_dims, param_type


def load_prms_output(output_data_path, csvfiles, convert=True, verbose=False):
    templist = []
    for csvname in csvfiles:
        fpath = os.path.join(output_data_path, csvname)
        tdf = pd.read_csv(
            fpath,
            parse_dates=["Date"],
            index_col=["Date"],
            dtype=float,
            float_precision="high",
        )
        colname = os.path.splitext(csvname)[0]
        tdf.columns = [colname]
        tdf.index.names = ["date"]
        templist.append(tdf)

    # concatenate individual dataframes
    df = pd.concat([v for v in templist], axis=1)

    # unit conversion
    if convert:
        unit_conversion(df, verbose=verbose)

    return df


def load_prms_statscsv(fname, convert=True, verbose=False):
    # read stats.csv
    with open(fname) as f:
        line = f.readline()
        colnames = line.strip().split(",")
        line = f.readline()
        units = line.strip().split(",")
    df = pd.read_csv(
        fname,
        skiprows=[1],
        parse_dates=["Date"],
        index_col=["Date"],
        dtype=float,
        float_precision="high",
    )
    df.index.names = ["date"]

    # unit conversion
    if convert:
        unit_conversion(df, verbose=verbose)

    return df


def load_wbl_output(output_data_path, convert=True, verbose=False):
    wbl_outflows = {
        "soilzone.wbal": (
            "perv ET",
            "sz2gw",
            "interflow",
            "soil2gw",
            "downflow",
            "pref flow",
            "pfr dunn",
            "dunnian gvr",
            "lake evap",
        ),
        "intcp.wbal": (
            "Netppt",
            "Intcpevap",
        ),
        "srunoff_smidx.wbal": (
            "Sroff",
            "Infil",
            "Impervevap",
            "Dprst_evap",
            "Dprst_seep",
            "Perv Sro",
            "Imperv Sro",
            "Dprst Sro",
            "CFGI Sro",
        ),
        "snowcomp.wbal": (
            "Snowmelt",
            "Snowevap",
        ),
        "gwflow.wbal": (
            "GW flow",
            "GW sink",
            "downflow",
        ),
    }
    wbl_files = [
        os.path.join(output_data_path, file)
        for file in os.listdir(output_data_path)
        if file.endswith(".wbal")
    ]
    df_dict = {}
    for fpth in wbl_files:
        key = os.path.basename(fpth)
        wbal_type = key.replace(".wbal", "")
        df = pd.read_fwf(fpth, parse_dates=["Date"], index_col=["Date"])

        # change sign of outflows
        for col in wbl_outflows[key]:
            df[col] *= -1.0

        # rename columns
        column_dict = {}
        for column in df.columns:
            column_out = "_".join(column.lower().split())
            column_dict[column] = f"{wbal_type}_{column_out}"
        df.rename(columns=column_dict, inplace=True)

        # unit conversion
        if convert:
            unit_conversion(df, verbose=verbose)

        # add dataframe to the dataframe dictionary
        df_dict[key] = df

    return pd.concat([v for k, v in df_dict.items()], axis=1)
