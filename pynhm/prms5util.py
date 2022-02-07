import datetime
import functools
import os

import numpy as np
import pandas as pd


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


def load_prms_parameters(pfn):
    line_num = 0
    vals = {}
    dims = {}
    param_dims = {}
    param_type = {}

    with open(pfn) as f:
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
                print("**** read parameters exception line = ", line)
                print(
                    "**** read parameters exception line_num = ", str(line_num)
                )
                print("**** Unexpected error:", sys.exc_info()[0])

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
                        print("parameter ", param_name, " is already in ", pfn)
                    else:
                        vals[param_name] = vs

            except:
                raise ValueError(
                    f"read parameters exception line_num = {line_num}"
                )

    return dims, vals, param_dims, param_type


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
