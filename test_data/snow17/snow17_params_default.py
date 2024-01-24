# Create default snow17 parameters for the the snow17 test case and for the
# drb_2yr domain.

# The default values are taken from here
# https://github.com/UW-Hydro/tonic/blob/master/tonic/models/snow17/snow17.py
# which mostly cites Shamir and Georgakakos 2007 default values for American
# River basin in Ca.
import pathlib as pl

import numpy as np
import pywatershed as pws

snow17_defaults = {
    "lat": 50,
    "elevation": 0,
    "scf": 1.0,
    "rvs": 1,
    "uadj": 0.04,
    "mbase": 1.0,
    "mfmax": 1.05,
    "mfmin": 0.6,
    "tipm": 0.1,
    "nmf": 0.15,
    "plwhc": 0.04,
    "pxtemp": 1.0,
    "pxtemp1": -1.0,
    "pxtemp2": 3.0,
}

dom_dict = {"drb": 765, "two": 2}

for ddk, ddv in dom_dict.items():
    nhru = ddv

    dims = {
        "nhru": nhru,
    }

    coords = {
        "hru_id": np.array(range(nhru)),
    }

    data_vars = {}
    metadata = {}
    metadata["hru_id"] = {"dims": ["nhru"]}
    for kk, vv in snow17_defaults.items():
        data_vars[kk] = np.ones(nhru) * vv
        metadata[kk] = {"dims": ["nhru"]}

    params = pws.parameters.PrmsParameters(
        dims=dims,
        coords=coords,
        data_vars=data_vars,
        metadata=metadata,
        validate=True,
    )

    netcdf_file = f"{ddk}_snow17_params_default.nc"
    if ddk == "drb":
        netcdf_file = pl.Path("../drb_2yr") / netcdf_file

    params.to_netcdf(netcdf_file, use_xr=True)
