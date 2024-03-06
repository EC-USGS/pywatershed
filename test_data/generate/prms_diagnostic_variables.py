import pathlib as pl

import pywatershed as pws
from pywatershed.constants import dnearzero, nearzero, nan, zero

import numpy as np
import xarray as xr

"""This module is for generating PRMS diagnostic variables"""

# For generating timeseries of previous states
previous_vars = {
    "gwres_stor": pws.PRMSGroundwater,
    "dprst_stor_hru": pws.PRMSRunoff,
    "hru_impervstor": pws.PRMSRunoff,
    "hru_intcpstor": pws.PRMSCanopy,
    "pref_flow_stor": pws.PRMSSoilzone,
    "slow_stor": pws.PRMSSoilzone,
    "soil_lower": pws.PRMSSoilzone,
    "soil_rechr": pws.PRMSSoilzone,
    "ssres_stor": pws.PRMSSoilzone,
    "freeh2o": pws.PRMSSnow,
    "pk_ice": pws.PRMSSnow,
}

prev_rename = {
    "gwres_stor": "gwres_stor_old",
    "hru_impervstor": "hru_impervstor_old",
    "hru_intcpstor": "hru_intcpstor_old",
    "dprst_stor_hru": "dprst_stor_hru_old",
}

change_rename = {}


def diagnose_simple_vars_to_nc(
    var_name: str,
    data_dir: pl.Path,
    control_file: pl.Path,
    domain_dir: pl.Path = None,
    output_dir: pl.Path = None,
):
    """Diagnose variables from PRMS output with only dependencies on var_name

    These variables are "simple" because they do not have any dependencies
    besides the "var_name" variable which is a PRMS output variable.

    This is not to say that PRMS dosent output diagnostic variables, it
    just dosent output these natively.

    Args:
        var_name: str name of the variable to create
        data_dir: where the netcdf file will be written, also where to look
            for necessary inputs to create the diagnostic variable.
        control_file: the pathlib representation of the control file for
            the run that produced the data.
        domain_dir: defaults to the parent dir of output_dir, this is where
            domain files for the domain can be found
        output_dir: defaults to data_dir unless otherwise specified

    """

    if domain_dir is None:
        domain_dir = data_dir.parent

    if output_dir is None:
        output_dir = data_dir

    nc_path = data_dir / f"{var_name}.nc"
    control = pws.Control.load_prms(control_file, warn_unused_options=False)
    param_file = control_file.parent / control.options["parameter_file"]

    if var_name in previous_vars.keys():
        # the _prev is the desired suffix but PRMS legacy is inconsistent
        current = xr.open_dataarray(nc_path)
        prev = current * nan
        prev[:] = np.roll(current.values, 1, axis=0)

        # the final value (-1) was wrapped to the zeroth position
        # get the initial conditions for the first time by initializing the
        # model This works based on the control file, so could handle restart.
        control.options = control.options | {
            "input_dir": domain_dir / "output",
        }
        params = pws.parameters.PrmsParameters.load(param_file)
        proc_class = previous_vars[var_name]
        inputs = {}
        for input_name in proc_class.get_inputs():
            # taken from process._initialize_var
            meta = pws.meta.find_variables([input_name])
            dims = [params.dims[vv] for vv in meta[input_name]["dims"]]
            init_type = meta[input_name]["type"]

            if len(dims) == 1:
                inputs[input_name] = np.full(dims, nan, dtype=init_type)
            else:
                inputs[input_name] = pws.base.timeseries.TimeseriesArray(
                    var_name=input_name,
                    control=control,
                    array=np.full(dims, nan, dtype=init_type),
                    time=np.arange(
                        np.datetime64("2017-01-01"),
                        np.datetime64("2017-01-03"),
                    ),
                )

        proc = proc_class(
            control=control,
            discretization=None,
            parameters=params,
            **inputs,
        )
        # these are the init_values
        prev[0, :] = proc[var_name]
        del proc

        change = current - prev

        # write the previous file
        if var_name in prev_rename.keys():
            nc_name = f"{prev_rename[var_name]}.nc"
        else:
            nc_name = f"{var_name}_prev.nc"
        out_nc_path = output_dir / nc_name
        prev.rename(out_nc_path.stem).to_netcdf(out_nc_path)
        assert nc_path.exists()

        # write the change file
        if var_name in change_rename.keys():
            nc_name = f"{change_rename[var_name]}.nc"
        else:
            nc_name = f"{var_name}_change.nc"
        out_nc_path = output_dir / nc_name
        change.rename(out_nc_path.stem).to_netcdf(out_nc_path)
        assert nc_path.exists()
        return nc_path

    if var_name in ["sroff", "ssres_flow", "gwres_flow"]:
        params = pws.parameters.PrmsParameters.load(param_file)
        ds = xr.open_dataset(nc_path)
        ds = ds.rename({var_name: f"{var_name}_vol"})
        ds = ds * params.data_vars["hru_in_to_cf"]
        ds.to_netcdf(output_dir / f"{var_name}_vol.nc")
        ds.close()

    if var_name == "infil":
        params = pws.parameters.PrmsParameters.load(param_file).parameters
        imperv_frac = params["hru_percent_imperv"]
        if (
            "dprst_flag" in control.options.keys()
            and control.options["dprst_flag"]
        ):
            dprst_frac = params["dprst_frac"]
        else:
            dprst_frac = zero
        perv_frac = 1.0 - imperv_frac - dprst_frac
        ds = xr.open_dataset(nc_path.with_suffix(".nc"))
        ds = ds.rename(infil="infil_hru")
        ds["infil_hru"] = ds["infil_hru"] * perv_frac
        ds.to_netcdf(output_dir / "infil_hru.nc")
        ds.close()


def diagnose_final_vars_to_nc(
    var_name: str,
    data_dir: pl.Path,
    control_file: pl.Path,
    domain_dir: pl.Path = None,
    output_dir: pl.Path = None,
):
    """Diagnose variables from multiple PRMS outputs

    In this case "var_name" is not a PRMS output variable, the variable is a
    new variable which requires much/all of the PRMS output to be already
    present in netcdf output format. This is why this is a final diagnostic

    Args:
        var_name: str name of the variable to create, not a PRMS variable.
        output_dir: where the netcdf file will be written, also where to look
            for necessary inputs to create the diagnostic variable.
        domain_dir: defaults to the parent dir of output_dir, this is where
            domain files (parameter & control) for the domain can be found

    """
    if domain_dir is None:
        domain_dir = data_dir.parent

    if output_dir is None:
        output_dir = data_dir

    control = pws.Control.load_prms(control_file, warn_unused_options=False)
    param_file = control_file.parent / control.options["parameter_file"]

    if var_name == "through_rain":
        out_file = output_dir / "through_rain.nc"

        data_vars = [
            "net_ppt",
            "pptmix_nopack",
            "snowmelt",
            "pkwater_equiv",
            "pk_ice_change",
            "freeh2o_change",
            "snow_evap",
            "net_snow",
            "net_rain",
        ]

        data = {}
        for vv in data_vars:
            data_file = data_dir / f"{vv}.nc"
            data[vv] = xr.open_dataarray(data_file)

        cond1 = data["net_ppt"] > zero
        cond2 = data["pptmix_nopack"] != 0
        cond3 = data["snowmelt"] < nearzero
        cond4 = data["pkwater_equiv"] < dnearzero
        cond5 = data["snow_evap"] < nearzero
        cond6 = data["net_snow"] < nearzero
        cond7 = data["snow_evap"] > -1 * (
            data["pk_ice_change"] + data["freeh2o_change"]
        )

        through_rain = data["net_rain"] * zero
        # these are in reverse order
        through_rain[:] = np.where(
            cond1 & cond3 & cond4 & cond6, data["net_rain"], zero
        )
        through_rain[:] = np.where(
            cond1 & cond3 & cond4 & cond5, data["net_ppt"], through_rain
        )
        through_rain[:] = np.where(
            cond1 & cond2, data["net_rain"], through_rain
        )

        # This condition does not exist in PRMS as far as I can tell
        # but is necessary for mass balance
        # This is when it rains on snow (no new snow) and then snow_evap
        # consumes the pack during the timestep.
        through_rain[:] = np.where(
            cond1 & cond6 & cond7,
            zero,
            through_rain,
        )

        through_rain.to_dataset(name="through_rain").to_netcdf(out_file)
        through_rain.close()

        for vv in data_vars:
            data[vv].close()

        assert out_file.exists()

    if var_name in [
        "channel_sroff_vol",
        "channel_ssres_flow_vol",
        "channel_gwres_flow_vol",
        "seg_lateral_inflow",
    ]:
        if var_name != "seg_lateral_inflow":
            return

        data_vars = [
            "sroff_vol",
            "ssres_flow_vol",
            "gwres_flow_vol",
            "seg_outflow",
            "seg_inflow",
        ]
        data = {}
        for vv in data_vars:
            data_file = data_dir / f"{vv}.nc"
            data[vv] = xr.open_dataarray(data_file)

        control = pws.Control.load_prms(
            control_file, warn_unused_options=False
        )
        s_per_time = control.time_step_seconds
        params = pws.parameters.PrmsParameters.load(param_file)

        ntimes = control.n_times
        nseg = params.dims["nsegment"]
        nhru = params.dims["nhru"]
        hru_segment = params.parameters["hru_segment"] - 1
        to_segment = (params.parameters["tosegment"] - 1).astype("int64")

        outflow_mask = np.full((len(to_segment)), False)
        for iseg in range(nseg):
            if to_segment[iseg] < 0:
                outflow_mask[iseg] = True
                continue

        seg_stor_change = (
            data["seg_inflow"] - data["seg_outflow"]
        ) * s_per_time

        channel_outflow_vol = (
            np.where(outflow_mask, data["seg_outflow"], zero)
        ) * s_per_time

        seg_lateral_inflow = np.zeros((ntimes, nseg))
        channel_sroff_vol = np.zeros((ntimes, nhru))
        channel_ssres_flow_vol = np.zeros((ntimes, nhru))
        channel_gwres_flow_vol = np.zeros((ntimes, nhru))

        for ihru in range(nhru):
            iseg = hru_segment[ihru]
            if iseg < 0:
                # This is bad, selective handling of fluxes is not cool,
                # mass is being discarded in a way that has to be coordinated
                # with other parts of the code.
                # This code shuold be removed evenutally.
                continue

            else:
                channel_sroff_vol[:, ihru] = data["sroff_vol"].values[:, ihru]
                channel_ssres_flow_vol[:, ihru] = data[
                    "ssres_flow_vol"
                ].values[:, ihru]
                channel_gwres_flow_vol[:, ihru] = data[
                    "gwres_flow_vol"
                ].values[:, ihru]

            # cubicfeet to cfs
            lateral_inflow = (
                channel_sroff_vol[:, ihru]
                + channel_ssres_flow_vol[:, ihru]
                + channel_gwres_flow_vol[:, ihru]
            ) / (s_per_time)

            seg_lateral_inflow[:, iseg] += lateral_inflow

        for hvar in [
            "channel_sroff_vol",
            "channel_ssres_flow_vol",
            "channel_gwres_flow_vol",
            "seg_lateral_inflow",
            "seg_stor_change",
            "channel_outflow_vol",
        ]:
            if hvar in [
                "seg_lateral_inflow",
                "seg_stor_change",
                "channel_outflow_vol",
            ]:
                dum = data["seg_outflow"].rename(hvar) * np.nan
            else:
                dum = data["sroff_vol"].rename(hvar) * np.nan

            dum[:, :] = locals()[hvar]
            out_file = data_dir / f"{hvar}.nc"
            dum.to_netcdf(out_file)
            print(out_file)
            assert out_file.exists()
