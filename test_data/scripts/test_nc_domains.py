from time import sleep

import numpy as np
import pytest
import xarray as xr

from pywatershed import CsvFile, Soltab
from pywatershed.constants import epsilon32, zero
from pywatershed.parameters import PrmsParameters


def test_csv_to_netcdf(csv_files):
    nc_path = csv_files.with_suffix(".nc")
    CsvFile(csv_files).to_netcdf(nc_path)
    assert nc_path.exists()


def test_soltab_to_netcdf(soltab_file):
    # the nhm_ids are not available in the solta_debug file currently, so get
    # them from the domain parameters
    domain_dir = soltab_file.parent
    # this is a hack that should probably rely on the yaml if/when this fails
    params = PrmsParameters.load(domain_dir / "myparam.param")
    nhm_ids = params.parameters["nhm_id"]

    output_dir = domain_dir / "output"
    soltab = Soltab(soltab_file, output_dir=output_dir, nhm_ids=nhm_ids)

    for var in soltab.variables:
        assert (output_dir / f"{var}.nc").exists()


@pytest.mark.order(after=["test_csv_to_netcdf"])
def test_csv_to_previous_netcdf(csv_files_prev):
    nc_name = csv_files_prev.with_suffix("").name + "_prev.nc"
    nc_path = csv_files_prev.parent / nc_name
    csv = CsvFile({nc_path.stem: csv_files_prev})
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


@pytest.mark.order(after=["test_csv_to_previous_netcdf"])
def test_misc_netcdf(misc_nc_files_input):
    if misc_nc_files_input.name in ["sroff", "ssres_flow", "gwres_flow"]:
        data_dir = misc_nc_files_input.parent
        domain_dir = misc_nc_files_input.parent.parent
        params = PrmsParameters.load(domain_dir / "myparam.param")
        var = misc_nc_files_input.name
        ds = xr.open_dataset(data_dir / f"{var}.nc")
        ds = ds.rename({var: f"{var}_vol"})
        ds = ds * params.data_vars["hru_in_to_cf"]
        ds.to_netcdf(data_dir / f"{var}_vol.nc")
        ds.close()

    if misc_nc_files_input.name == "infil":
        domain_dir = misc_nc_files_input.parent.parent
        params = PrmsParameters.load(domain_dir / "myparam.param").parameters
        imperv_frac = params["hru_percent_imperv"]
        dprst_frac = params["dprst_frac"]
        perv_frac = 1.0 - imperv_frac - dprst_frac
        ds = xr.open_dataset(misc_nc_files_input.with_suffix(".nc"))
        ds = ds.rename(infil="infil_hru")

        # not necessary
        # perv_frac = np.tile(perv_frac, (ds["infil_hru"].shape[0], 1))
        ds["infil_hru"] = ds["infil_hru"] * perv_frac

        ds.to_netcdf(misc_nc_files_input.parent / "infil_hru.nc")
        ds.close()

    assert True


@pytest.mark.order(after=["test_misc_netcdf", "test_csv_to_previous_netcdf"])
def test_misc_final(misc_nc_final_input):
    if misc_nc_final_input.name == "through_rain":
        data_dir = misc_nc_final_input.parent
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

            if False:
                result = xr.open_dataset(data_file)[vv]
            else:
                # None of this case should be necessary if the test
                # really runs after the other tests as specified
                # But the following will apparently trigger the race conditions
                # rm ../*/output/*.nc; pytest -vv -n=auto test_nc_domains.py
                result = None
                while result is None:
                    try:
                        result = xr.open_dataset(data_file)[vv]
                    except FileNotFoundError:
                        sleep(0.1)
                        pass
                    except OSError:
                        sleep(0.1)
                        pass

            data[vv] = result

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
            misc_nc_final_input.parent / "through_rain.nc"
        )
        through_rain.close()
        for vv in data_vars:
            data[vv].close()

    if misc_nc_final_input.name in [
        "pk_ice",
        "freeh2o",
        "hru_impervstor",
        "dprst_stor_hru",
    ]:
        var_name = misc_nc_final_input.name
        var_prev_name = f"{var_name}_prev"
        parent_dir = misc_nc_final_input.parent
        current_file = parent_dir / f"{var_name}.nc"
        prev_file = parent_dir / f"{var_prev_name}.nc"
        current = xr.open_dataset(current_file)[var_name]
        previous = xr.open_dataset(prev_file)[var_prev_name]
        var_change_name = f"{var_name}_change"
        change = (current - previous).rename(var_change_name).to_dataset()
        change_file = parent_dir / f"{var_change_name}.nc"
        change.to_netcdf(change_file)
        assert change_file.exists()
