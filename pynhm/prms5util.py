import os

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
    "hru_ppt": inch_to_meter,
    # "net_ppt": inch_to_meter,
    # "soil_moist_ante": inch_to_meter,
    # "hru_sroffi": inch_to_meter,
    # "hru_sroffp": inch_to_meter,
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
