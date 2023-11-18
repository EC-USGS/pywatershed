import sys

import holoviews as hv
import hvplot.pandas
import hvplot.xarray  # noqa
import numpy as np
import pandas as pd
import panel as pn
import pywatershed
import xarray as xr

sys.path.append("../common")
import metrics

pd.options.plotting.backend = "holoviews"

fileish = pywatershed.constants.fileish


def get_diff_tol_df(
    var_df: pd.DataFrame, var_name: str, tol: np.float64
) -> pd.DataFrame:
    """A dataframe of diffs w.r.t tol
    Define diffs as both abs diffs or relative abs diffs exceede tolerance.
    Data frame includes mask for input dataframe.
    """
    var_df["diff"] = var_df[f"{var_name}_pynhm"] - var_df[f"{var_name}_prms"]
    var_df["rel_diff"] = var_df["diff"] / var_df[f"{var_name}_prms"]
    diff_mask = abs(var_df["diff"]) > tol
    rel_diff_mask = abs(var_df["rel_diff"]) > tol
    diff_df = var_df.copy()
    diff_df["diff"] = np.where(
        diff_mask & rel_diff_mask, var_df["diff"], np.nan
    )
    diff_df["rel_diff"] = np.where(
        diff_mask & rel_diff_mask, var_df["rel_diff"], np.nan
    )
    diff_df["mask"] = np.where(diff_mask & rel_diff_mask, True, False)
    diff_df = diff_df.drop(columns=[f"{var_name}_pynhm", f"{var_name}_prms"])
    return diff_df


def get_proc_var_dict(model: pywatershed.Model, ds: xr.Dataset) -> dict:
    """Dictionary of process: var for a model given a dataset containing actual
    model,obs variables from that model to filter available variables
    """
    avail_vars = list(ds.variables)
    drop_vars = ["doy", "nhm_id", "time", "nhm_seg"]
    avail_vars = list(
        set(
            [
                "_".join(vv.split("_")[0:-1])
                for vv in avail_vars
                if vv not in drop_vars
            ]
        )
    )
    proc_var_dict = {}
    for proc_key, proc_val in model.processes.items():
        proc_var_dict[proc_key] = [
            var for var in proc_val.get_variables() if var in avail_vars
        ]
        if not len(proc_var_dict[proc_key]):
            del proc_var_dict[proc_key]
    return proc_var_dict


def get_comp_ds(model: pywatershed.Model, var_name: str) -> xr.Dataset:
    """Get a dataset for comparing variables from pywatershed to PRMS"""
    prms_file = model.input_dir / f"{var_name}.nc"
    pynhm_file = model._netcdf_dir / f"{var_name}.nc"
    if not prms_file.exists():
        print(f"PRMS file does not exist: {prms_file}")
    if not pynhm_file.exists():
        print(f"pynhm file does not exist: {pynhm_file}")
    if (not prms_file.exists()) or (not pynhm_file.exists()):
        return None

    comp_ds = xr.merge(
        [
            xr.open_dataset(prms_file, decode_timedelta=False).rename(
                {var_name: f"{var_name}_prms"}
            ),
            xr.open_dataset(pynhm_file, decode_timedelta=False).rename(
                {var_name: f"{var_name}_pynhm"}
            ),
        ]
    )
    comp_ds.attrs["Description"] = f"Variable comparison for PRMS and pynhm"
    comp_ds.attrs[var_name] = pywatershed.meta.get_vars(var_name)[var_name]

    return comp_ds


def get_stat_location_id(stat_df, stat_name, tap_value, space_coord):
    """Get the spatial identifier nearest the tapped stat/error value"""
    stat_df["diff"] = stat_df[stat_name] - tap_value
    return stat_df[abs(stat_df["diff"]) == abs(stat_df["diff"]).min()][
        space_coord
    ].values[0]


def err_panel(
    model: pywatershed.Model,
    fig_width: int = 600,
    tol: np.float64 = np.finfo(np.single).resolution,
    gis_dir: fileish = None,
    err_data_hru: np.ndarray = None,
    err_data_seg: np.ndarray = None,
) -> pn.Column:
    """A panel dashboard for exploring error distributions and drilling down
    to individual locations for details on their timeseries and errors
    """
    var_list = [
        var
        for proc_key, proc_val in model.processes.items()
        for var in proc_val.get_variables()
    ]
    # var_dict = {proc_key: list(proc_val.get_variables()) for proc_key, proc_val in model.processes.items()}

    data_list = [get_comp_ds(model, var) for var in var_list]
    data_ds = xr.merge(
        [data for data in data_list if data], combine_attrs="drop_conflicts"
    ).load()
    # del data_list

    var_name = pn.widgets.Select(
        groups=get_proc_var_dict(model, data_ds), width=200, name="Variable"
    )
    stat_name = pn.widgets.Select(
        options=list(metrics.stat_dict.keys()),
        width=200,
        name="Error Statistic",
    )

    space_coord = "nhm_id" if "nhm_id" in data_ds.coords else "nhm_seg"
    space_id_sel = pn.widgets.Select(
        options=sorted(data_ds[space_coord].values.tolist()),
        width=200,
        name="Spatial Coord ID",
    )

    def update_space_id_value(event):
        # with open("hello.txt", mode="w") as file:
        #    file.write(f"{event.new}")
        space_coord, time_coord = get_var_info(var_name.value)
        stat_df = get_stat_df(stat_name.value, var_name.value)
        space_id = get_stat_location_id(
            stat_df, stat_name.value, event.new, space_coord
        )
        space_id_sel.value = space_id
        return

    pn.extension()
    tap = hv.streams.Tap(y=0)

    def get_var_info(var_name):
        """Get space and time coordinate info for a variable
        Depends on: data_ds being in scope
        """
        var_ds = data_ds[[f"{var_name}_pynhm", f"{var_name}_prms"]]
        time_coord = "time" if "time" in var_ds.coords else "doy"
        space_coord = "nhm_id" if "nhm_id" in var_ds.coords else "nhm_seg"
        return space_coord, time_coord

    def get_stat_df(stat_name, var_name):
        """Get a data frame of the error statistic
        Depends on: get_var_info, data_ds in scope
        """
        space_coord, time_coord = get_var_info(var_name)
        # var_ds = data_ds[[f"{var_name}_pynhm", f"{var_name}_prms"]]

        stat_df = (
            metrics.stat_dict[stat_name](
                mod=data_ds[f"{var_name}_pynhm"],
                obs=data_ds[f"{var_name}_prms"],
                dim=time_coord,
            )
            .to_dataframe(name=stat_name)
            .reset_index()
        )
        if space_coord == "nhm_id":
            err_data_hru[:] = stat_df[stat_name]
        else:
            err_data_seg[:] = stat_df[stat_name]

        stat_df["x"] = ""
        stat_df = stat_df.set_index("x")

        return stat_df

    def get_ts_df(var_name, space_coord, space_id):
        """Get a timeseries dataframe for a tapped coordinate location"""
        return (
            data_ds[[f"{var_name}_prms", f"{var_name}_pynhm"]]
            .sel(indexers={space_coord: space_id})
            .drop(space_coord)
            .to_dataframe()
        )

    @pn.depends(stat_name, var_name)
    def err_plot(stat_name, var_name):
        space_coord, time_coord = get_var_info(var_name)
        stat_df = get_stat_df(stat_name, var_name)

        var_meta = pywatershed.meta.get_vars(var_name)[var_name]
        title = (
            f"{var_meta['desc']}\n"
            f"domain 1979-1980: {stat_name.upper()} by {space_coord}\n"
            f"(N = n_stat_hru) ({var_meta['units']})"
        )
        violin = stat_df.hvplot.violin(
            y=stat_name,
            ylabel=stat_name.upper(),
            xlabel="HRU",
            width=fig_width,
            title=title,
        )
        scatter = stat_df.hvplot.scatter(
            y=stat_name,
            x="x",
            hover_cols=[stat_name, space_coord],
            use_index=False,
            c="red",
            size=15.2,
        )
        tap.source = scatter
        return violin * scatter

    @pn.depends(space_id_sel, var_name)
    def timeseries(space_id, var_name):
        space_coord, time_coord = get_var_info(var_name)
        ts_df = get_ts_df(var_name, space_coord, space_id)
        title = f"{space_coord} = {space_id}"
        return ts_df.hvplot(title=title, shared_axes=False)

    @pn.depends(space_id_sel, var_name)
    def scatter(space_id, var_name):
        space_coord, time_coord = get_var_info(var_name)
        ts_df = get_ts_df(var_name, space_coord, space_id)
        title = f"{space_coord} = {space_id}"
        return ts_df.hvplot.scatter(
            x=f"{var_name}_prms",
            y=f"{var_name}_pynhm",
            hover_cols=["time"],
            title=title,
            width=int(fig_width / 2),
            shared_axes=False,
        )

    @pn.depends(space_id_sel, var_name)
    def timeseries_diffs(space_id, var_name):
        space_coord, time_coord = get_var_info(var_name)
        diff_df = get_diff_tol_df(
            get_ts_df(var_name, space_coord, space_id), var_name, tol
        )
        diff_df = diff_df.drop(columns=["mask"])
        title = f"{space_coord} = {space_id}"
        return diff_df.hvplot(title=title, shared_axes=False)

    @pn.depends(space_id_sel, var_name)
    def scatter_diffs(space_id, var_name):
        space_coord, time_coord = get_var_info(var_name)
        ts_df = get_ts_df(var_name, space_coord, space_id)
        diff_df = get_diff_tol_df(ts_df, var_name, tol)
        ts_df["diff"] = diff_df["diff"]
        ts_df["rel_diff"] = diff_df["rel_diff"]
        plot_df = ts_df[diff_df["mask"]]
        title = f"{space_coord} = {space_id}"
        return plot_df.hvplot.scatter(
            x=f"{var_name}_prms",
            y=f"{var_name}_pynhm",
            hover_cols=["time", "diff", "rel_diff"],
            title=title,
            width=int(fig_width / 2),
            shared_axes=False,
        )

    tap.param.watch(update_space_id_value, "y")
    tap.param.trigger("y")

    plt = pn.Column(
        pn.Row(stat_name, var_name, space_id_sel),
        pn.Row(err_plot),
        pn.Row(timeseries, scatter),
        pn.Row(timeseries_diffs, scatter_diffs),
    )
    return plt
