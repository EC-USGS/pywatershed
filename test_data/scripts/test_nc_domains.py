from pathlib import Path
from filelock import FileLock

import numpy as np

import pywatershed as pws
from pywatershed import CsvFile, Soltab
from pywatershed.parameters import PrmsParameters
from pywatershed.constants import epsilon32, nan, zero

import pytest
import xarray as xr

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


@pytest.fixture
def netcdf_file(csv_file) -> Path:
    """Convert CSV files from model output to NetCDF"""

    var_name = csv_file.stem
    nc_path = csv_file.with_suffix(".nc")
    data_dir = csv_file.parent
    domain_dir = data_dir.parent
    CsvFile(csv_file).to_netcdf(nc_path)
    assert nc_path.exists()

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
        params = PrmsParameters.load(domain_dir / "myparam.param")
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
        nc_path = data_dir / nc_name
        prev.rename(nc_path.stem).to_netcdf(nc_path)
        assert nc_path.exists()

        # write the change file
        if var_name in change_rename.keys():
            nc_name = f"{change_rename[var_name]}.nc"
        else:
            nc_name = f"{var_name}_change.nc"
        nc_path = data_dir / nc_name
        change.rename(nc_path.stem).to_netcdf(nc_path)
        assert nc_path.exists()

    if var_name in ["sroff", "ssres_flow", "gwres_flow"]:
        params = PrmsParameters.load(domain_dir / "myparam.param")
        ds = xr.open_dataset(nc_path)
        ds = ds.rename({var_name: f"{var_name}_vol"})
        ds = ds * params.data_vars["hru_in_to_cf"]
        ds.to_netcdf(data_dir / f"{var_name}_vol.nc")
        ds.close()

    if var_name == "infil":
        params = PrmsParameters.load(domain_dir / "myparam.param").parameters
        imperv_frac = params["hru_percent_imperv"]
        dprst_frac = params["dprst_frac"]
        perv_frac = 1.0 - imperv_frac - dprst_frac
        ds = xr.open_dataset(nc_path.with_suffix(".nc"))
        ds = ds.rename(infil="infil_hru")
        ds["infil_hru"] = ds["infil_hru"] * perv_frac
        ds.to_netcdf(data_dir / "infil_hru.nc")
        ds.close()

    return nc_path


def make_netcdf_files(netcdf_file):
    print(f"Created NetCDF from CSV: {netcdf_file}")


@pytest.fixture(scope="session")
def soltab_netcdf_file(tmp_path_factory, soltab_file) -> Path:
    """Convert soltab files to NetCDF, one file for each variable"""

    # the nhm_ids are not available in the solta_debug file currently, so get
    # them from the domain parameters
    domain_dir = soltab_file.parent
    output_dir = domain_dir / "output"
    output_dir.mkdir(exist_ok=True)

    # with FileLock(root_tmpdir / "soltab_nc.lock"):
    #     yield  # postpone the work until session cleanup

    # this is a hack that should probably rely on the yaml if/when this
    # fails
    params = PrmsParameters.load(domain_dir / "myparam.param")
    nhm_ids = params.parameters["nhm_id"]
    soltab = Soltab(soltab_file, nhm_ids=nhm_ids)

    soltab.to_netcdf(output_dir=output_dir)
    print(f"Created NetCDF files from soltab file {soltab_file}:")

    for var in soltab.variables:
        nc_path = output_dir / f"{var}.nc"
        assert nc_path.exists()


def make_soltab_netcdf_files(soltab_netcdf_file):
    print(f"Creating NetCDF files for soltab file {soltab_netcdf_file}")


@pytest.fixture(scope="session")
def final_netcdf_file(tmp_path_factory, simulation) -> Path:
    """Create the final NetCDF file (through_rain.nc) from other NetCDFs"""

    root_tmpdir = tmp_path_factory.getbasetemp().parent
    data_dir = Path(simulation["ws"]) / "output"
    data_dir.mkdir(exist_ok=True)
    nc_path = data_dir / "through_rain.nc"

    with FileLock(root_tmpdir / "final_nc.lock"):
        yield  # do this in session cleanup

        if nc_path.is_file():
            return

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

        through_rain.to_dataset(name="through_rain").to_netcdf(
            data_dir / "through_rain.nc"
        )
        through_rain.close()

        for vv in data_vars:
            data[vv].close()


def make_final_netcdf_files(final_netcdf_file):
    print(f"Creating final NetCDF file {final_netcdf_file}")
