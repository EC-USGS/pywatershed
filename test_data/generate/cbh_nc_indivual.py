import pathlib as pl

import xarray as xr

domains = ["drb_2yr", "ucb_2yr", "hru_1"]

cbh_files = [pl.Path(f"../{domain}/cbh.nc") for domain in domains]

vars = ["prcp", "tmin", "tmax"]

for cbh_file in cbh_files:
    cbh = (
        xr.open_dataset(cbh_file)
        .swap_dims({"hru": "nhm_id"})
        .rename(datetime="time")
        .set_coords("time")
    )

    for var in vars:
        cbh_var = cbh[var]
        cbh_var.to_netcdf(cbh_file.parent / f"{var}.nc")


print("")
