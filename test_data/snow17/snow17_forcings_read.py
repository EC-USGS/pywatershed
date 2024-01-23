from datetime import datetime, timedelta
import pathlib as pl

import numpy as np
import pandas as pd
import xarray as xr

mm2in = 1.0 / (2.54 * 10.0)

forcing_input_dir = pl.Path(".")
forcing_files = sorted(
    [
        "forcing.snow17bmi.HHWM8IL.csv",
        "forcing.snow17bmi.HHWM8IU.csv",
    ]
)

forcings = []
count = 0
for ff in forcing_files:
    hru_id = ff.split(".")[-2]
    df = pd.read_csv(forcing_input_dir / ff)
    times = np.arange(
        datetime(df.year[0], df.mo[0], df.dy[0]),
        datetime(df.year.values[-1], df.mo.values[-1], df.dy.values[-1])
        + timedelta(1),
        timedelta(1),
    )
    ntimes = len(times)
    assert len(df) == ntimes

    ds = xr.Dataset(
        data_vars=dict(
            tavgc=(["time", "hru_id"], df.tavg_degc.values.reshape(ntimes, 1)),
            net_ppt=(
                ["time", "hru_id"],
                (df["prec_mm_s-1"].values * 24 * 60 * 60 * mm2in).reshape(
                    ntimes, 1
                ),
            ),
        ),
        coords=dict(time=times, hru_id=np.array([count])),
        attrs=dict(description=f"Snow17 forcing data for HRU {hru_id}"),
    )

    ds.tavgc.attrs["units"] = "deg C"
    ds.net_ppt.attrs["units"] = "in (/day)"

    # ds
    nc_filename = f"snow17_forcing_data_{hru_id}.nc"
    ds.squeeze().to_netcdf(nc_filename)
    forcings += [ds]
    count += 1


combined = xr.concat(forcings, dim="hru_id")
assert (combined.tavgc[:, 0] == forcings[0].tavgc).all()
assert (combined.tavgc[:, 1] == forcings[1].tavgc).all()
assert (combined.net_ppt[:, 0] == forcings[0].net_ppt).all()
assert (combined.net_ppt[:, 1] == forcings[1].net_ppt).all()

f_dir = pl.Path("two_forcing_data")
f_dir.mkdir(exist_ok=True)

combined.net_ppt.to_dataset().to_netcdf(f_dir / "net_ppt.nc")
combined.tavgc.to_dataset().to_netcdf(f_dir / "tavgc.nc")
