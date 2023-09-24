import pathlib as pl

import pywatershed as pws
from pywatershed.constants import epsilon32, nan, zero

import numpy as np
import xarray as xr

"""This module is for generating PRMS diagnostic variables"""

# For generating timeseries of previous states
previous_vars = {
    "gwres_stor": pws.PRMSGroundwater,
    "dprst_stor_hru": pws.PRMSRunoff,
    "hru_impervstor": pws.PRMSRunoff,
    "pref_flow_stor": pws.PRMSSoilzone,
    "slow_stor": pws.PRMSSoilzone,
    "soil_lower": pws.PRMSSoilzone,
    "soil_moist": pws.PRMSSoilzone,
    "soil_rechr": pws.PRMSSoilzone,
    "ssres_stor": pws.PRMSSoilzone,
    "freeh2o": pws.PRMSSnow,
    "pk_ice": pws.PRMSSnow,
}

prev_rename = {
    "gwres_stor": "gwres_stor_old",
}

change_rename = {}


def diagnose_simple_vars_to_nc(
    var_name: str,
    data_dir: pl.Path,
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
        domain_dir: defaults to the parent dir of output_dir, this is where
            domain files (parameter & control) for the domain can be found
        output_dir: defaults to data_dir unless otherwise specified

    """

    if domain_dir is None:
        domain_dir = data_dir.parent

    if output_dir is None:
        output_dir = data_dir

    nc_path = data_dir / f"{var_name}.nc"

    if var_name in previous_vars.keys():
        # the _prev is the desired suffix but PRMS legacy is inconsistent
        current = xr.open_dataarray(nc_path)
        prev = current * nan
        prev[:] = np.roll(current.values, 1, axis=0)

        # the final value (-1) was wrapped to the zeroth position
        # get the initial conditions for the first time by initializing the
        # model This works based on the control file, so could handle restart.
        control = pws.Control.load(domain_dir / "control.test")
        control.options = control.options | {
            "input_dir": domain_dir / "output",
        }
        params = pws.parameters.PrmsParameters.load(
            domain_dir / "myparam.param"
        )
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
        params = pws.parameters.PrmsParameters.load(
            domain_dir / "myparam.param"
        )
        ds = xr.open_dataset(nc_path)
        ds = ds.rename({var_name: f"{var_name}_vol"})
        ds = ds * params.data_vars["hru_in_to_cf"]
        ds.to_netcdf(output_dir / f"{var_name}_vol.nc")
        ds.close()

    if var_name == "infil":
        params = pws.parameters.PrmsParameters.load(
            domain_dir / "myparam.param"
        ).parameters
        imperv_frac = params["hru_percent_imperv"]
        dprst_frac = params["dprst_frac"]
        perv_frac = 1.0 - imperv_frac - dprst_frac
        ds = xr.open_dataset(nc_path.with_suffix(".nc"))
        ds = ds.rename(infil="infil_hru")
        ds["infil_hru"] = ds["infil_hru"] * perv_frac
        ds.to_netcdf(output_dir / "infil_hru.nc")
        ds.close()


def diagnose_final_vars_to_nc(
    var_name: str,
    data_dir: pl.Path,
    domain_dir: pl.Path = None,
    output_dir: pl.Path = None,
):
    """Diagnose variables from multiple PRMS outputs

    In this case "var_name" is not a PRMS output variable, the variable is a
    new variable which requires much/all of the PRMS output to be already
    present in netcdf output format. This is why this is a final diagnostic

    Currently only: "through_rain"

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

    if var_name == "through_rain":
        out_file = output_dir / "through_rain.nc"

        data_vars = [
            "net_ppt",
            "pptmix_nopack",
            "snowmelt",
            "pkwater_equiv",
            "snow_evap",
            "net_snow",
            "net_rain",
        ]

        data = {}
        for vv in data_vars:
            data_file = data_dir / f"{vv}.nc"
            data[vv] = xr.open_dataset(data_file)[vv]

        nearzero = 1.0e-6

        cond1 = data["net_ppt"] > zero
        cond2 = data["pptmix_nopack"] != 0
        cond3 = data["snowmelt"] < nearzero
        cond4 = data["pkwater_equiv"] < epsilon32
        cond5 = data["snow_evap"] < nearzero
        cond6 = data["net_snow"] < nearzero

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

        through_rain.to_dataset(name="through_rain").to_netcdf(out_file)
        through_rain.close()

        for vv in data_vars:
            data[vv].close()

        assert out_file.exists()
