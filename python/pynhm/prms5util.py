import os
import functools
import numpy as np
import pandas as pd
import datetime


def load_prms_input(input_data_path, datanames, filenames):
    # load prms input
    templist = []
    for dataname, cbhname in zip(datanames, filenames):
        fpath = os.path.join(input_data_path, cbhname)
        print(f"Loading {fpath}")
        filelist = []
        with open(fpath) as f:
            for i, line in enumerate(f):
                if i > 2:
                    ll = line.strip().split()
                    yr = int(ll[0])
                    mo = int(ll[1])
                    da = int(ll[2])
                    data = float(ll[-1])
                    dt = datetime.datetime(yr, mo, da)
                    filelist.append([dt, data])
        df = pd.DataFrame(filelist, columns=["date", dataname])
        templist.append(df)
    input_df = functools.reduce(
        lambda left, right: pd.merge(left, right, on="date"), templist
    )
    input_df.set_index("date")
    return input_df


def load_prms_output(output_data_path, csvfiles):
    templist = []
    for csvname in csvfiles:
        fpath = os.path.join(output_data_path, csvname)
        df = pd.read_csv(fpath)
        colname = os.path.splitext(csvname)[0]
        df.columns = ["date", colname]
        templist.append(df)
    output_df = functools.reduce(
        lambda left, right: pd.merge(left, right, on="date"), templist
    )
    output_df["date"] = pd.to_datetime(output_df["date"])
    return output_df


def load_prms_statscsv(fname):
    # read stats.csv
    with open(fname) as f:
        line = f.readline()
        colnames = line.strip().split(",")
        line = f.readline()
        units = line.strip().split(",")
    prms_stats_df = pd.read_csv(fname, skiprows=2, header=None)
    prms_stats_df.columns = colnames
    prms_stats_df.rename(columns={"Date": "date"}, inplace=True)
    prms_stats_df.set_index("date")
    prms_stats_df["date"] = pd.to_datetime(prms_stats_df["date"])
    return prms_stats_df


def load_wbl_output(output_data_path):
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
        for col in wbl_outflows[key]:
            df[col] *= -1.0
        column_dict = {}
        for column in df.columns:
            column_out = "_".join(column.lower().split())
            column_dict[column] = f"{wbal_type}_{column_out}"
        df.rename(columns=column_dict, inplace=True)
        df_dict[key] = df

    return pd.concat([v for k, v in df_dict.items()], axis=1)
