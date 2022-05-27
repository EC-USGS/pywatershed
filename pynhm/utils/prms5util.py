import os
import pathlib as pl
from typing import Tuple

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


# JLM: what is the fate of this? I deprecated in my preprocessing
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
                    min = int(ll[3])
                    sec = int(ll[4])
                    data = ll[5:]
                    # dt = datetime.datetime(yr, mo, da)
                    filelist += (
                        f"{da:02d}/{mo:02d}/{yr:04d},{','.join(data)}\n"
                    )
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
    # JLM: for prms_summary.csv?
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


def load_nhru_output_csv(fname, convert=True, verbose=False):
    # read stats.csv
    df = pd.read_csv(
        fname,
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


def load_soltab_debug(file_path: pl.Path) -> Tuple[np.ndarray, np.ndarray]:
    """Load the PRMS soltab_debug output file.

    With `print_debug` set to 5 in the control file, PRMS 5.2.1 prings the
    `soltab_debug` file in the run directory. This function parses it.

    Args:
        file_path: a pathlib.Path object representing the file

    Returns:
        Tuple of length 2 of np.ndarrays of shape [366, n_hru], the first
        for potential solar radiation and the second the total number of
        hours of sun.
    """
    # the data are 8 characters wide, 13 per line
    # for 28 full lines and a final line of 2 (366 total)
    width = 17
    slices = [slice(ii * width, (ii + 1) * width) for ii in range(13)]
    not_data = ["", "\n"]

    with file_path.open("r") as file_open:
        lines = file_open.readlines()
        data_list = []
        for ll in lines:
            if ll.startswith(" HRU:"):
                hru_ind = int((ll.strip().replace(" ", "").split(":"))[1]) - 1
                data_list += [{}]
            elif ll.startswith(" ***Soltab_sunhrs***"):
                key = "sunhrs"
                data_list[hru_ind][key] = []
                accum_count = 0
            elif ll.startswith(" ***Soltab_potsw***"):
                key = "potsw"
                data_list[hru_ind][key] = []
                accum_count = 0
            elif ll.startswith(" ***Soltab_horad_potsw***"):
                key = "potsw_flat"
                data_list[hru_ind][key] = []
                accum_count = 0
            else:
                # There are two lines at the end of the file that are
                # not like the others. accum_count will skip those
                accum_count += 1
                if accum_count > 29:
                    continue
                str_vals = [ll[slice] for slice in slices]
                data_list[hru_ind][key] += [
                    float(str) for str in str_vals if str not in not_data
                ]

    n_hru = len(data_list)

    # check number of hrus (not known apriori)
    # assert n_hru == 765  # drb

    # check lengths of each
    for v0 in data_list:
        for k1, v1 in v0.items():
            assert len(v1) == 366

    potential_sw_rad = np.zeros((366, n_hru))
    potential_sw_rad_flat = np.zeros((366, n_hru))
    sun_hrs = np.zeros((366, n_hru))

    for hh in range(n_hru):
        potential_sw_rad[:, hh] = np.array(data_list[hh]["potsw"])
        potential_sw_rad_flat[:, hh] = np.array(data_list[hh]["potsw_flat"])
        sun_hrs[:, hh] = np.array(data_list[hh]["sunhrs"])

    return potential_sw_rad, potential_sw_rad_flat, sun_hrs
