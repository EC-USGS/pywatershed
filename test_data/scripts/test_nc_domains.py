from pathlib import Path
from filelock import FileLock

import numpy as np

from pywatershed import CsvFile, Soltab
from pywatershed.parameters import PrmsParameters
from pywatershed.constants import epsilon64, zero

import pytest
import xarray as xr


# For generating timeseries of previous states
previous_vars = [
    "dprst_stor_hru",
    "freeh2o",
    "hru_impervstor",
    "pk_ice",
    "pref_flow_stor",
    "slow_stor",
    "soil_lower",
    "soil_moist",
    "soil_rechr",
    "ssres_stor",
]


@pytest.fixture
def netcdf_file(csv_file) -> Path:
    """Convert CSV files from model output to NetCDF"""

    nc_path = csv_file.with_suffix(".nc")
    data_dir = csv_file.parent
    domain_dir = data_dir.parent
    CsvFile(csv_file).to_netcdf(nc_path)
    assert nc_path.exists()

    if csv_file.stem in previous_vars:
        nc_name = csv_file.stem + "_prev.nc"
        nc_path = data_dir / nc_name
        csv = CsvFile({nc_path.stem: csv_file})
        csv._get_data()  # why so private?
        orig_dates = csv._data["date"].copy()
        csv._data = np.roll(csv._data, 1, axis=0)
        csv._data["date"] = orig_dates
        # Here we will eventually want to supply the desired initial conditions
        for hh in range(len(csv._data[0])):
            if hh == 0:
                continue
            csv._data[0][hh] = np.zeros([1])[0]
        csv.to_netcdf(nc_path)
        assert nc_path.exists()

    if nc_path.stem in ["sroff", "ssres_flow", "gwres_flow"]:
        params = PrmsParameters.load(domain_dir / "myparam.param")
        var = nc_path.stem
        ds = xr.open_dataset(nc_path)
        ds = ds.rename({var: f"{var}_vol"})
        ds = ds * params.data_vars["hru_in_to_cf"]
        ds.to_netcdf(data_dir / f"{var}_vol.nc")
        ds.close()

    if nc_path.stem == "infil":
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
    root_tmpdir = tmp_path_factory.getbasetemp().parent
    nc_paths = []

    with FileLock(root_tmpdir / "soltab_nc.lock"):
        yield soltab_file  # wait til session cleanup

        # this is a hack that should probably rely on the yaml if/when this
        # fails
        params = PrmsParameters.load(domain_dir / "myparam.param")
        nhm_ids = params.parameters["nhm_id"]
        soltab = Soltab(soltab_file, nhm_ids=nhm_ids)

        already_created = (output_dir / "soltab_potsw.nc").is_file()
        if not already_created:
            soltab.to_netcdf(output_dir=output_dir)
            print(f"Created NetCDF files from soltab file {soltab_file}:")

        for var in soltab.variables:
            nc_path = output_dir / f"{var}.nc"
            nc_paths.append(nc_path)
            assert nc_path.exists()
            print(f"\t{nc_path}")


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
        yield nc_path  # wait til session cleanup

        if nc_path.is_file():
            return

        data_vars = [
            "net_rain",
            "pk_ice_prev",
            "freeh2o_prev",
            "newsnow",
            "pptmix_nopack",
        ]

        data = {}
        for vv in data_vars:
            data_file = data_dir / f"{vv}.nc"
            result = xr.open_dataset(data_file)[vv]
            data[vv] = result

        try:
            wh_through = (
                ((data["pk_ice_prev"] + data["freeh2o_prev"]) <= epsilon64)
                & ~(data["newsnow"] == 1)
            ) | (data["pptmix_nopack"] == 1)

            through_rain = data["net_rain"].copy()
            through_rain[:] = np.where(wh_through, data["net_rain"], zero)
            through_rain.to_dataset(name="through_rain").to_netcdf(nc_path)
            through_rain.close()
            assert nc_path.exists()
        finally:
            for vv in data_vars:
                data[vv].close()


def make_final_netcdf_files(final_netcdf_file):
    print(f"Creating final NetCDF file {final_netcdf_file}")
