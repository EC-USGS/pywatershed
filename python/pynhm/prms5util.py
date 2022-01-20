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
    input_df = functools.reduce(lambda left, right: pd.merge(left, right, on="date"), templist)
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
    output_df = functools.reduce(lambda left, right: pd.merge(left, right, on="date"), templist)
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

