# this creates individual netcdf file from individual cbh files
import pathlib as pl

from pywatershed import PrmsParameters
from pywatershed.utils.cbh_utils import cbh_files_to_netcdf

dom_dir = pl.Path("../../../data/pywatershed/conus_2yr/")
param_file = dom_dir / "myparam.param"
params = PrmsParameters.load(dom_dir / "myparam.param")

cbh_files = {
    "prcp": dom_dir / "prcp.cbh",
    "rhavg": dom_dir / "rhavg.cbh",
    "tmax": dom_dir / "tmax.cbh",
    "tmin": dom_dir / "tmin.cbh",
}

# Do them individually so we dont have to separate later.
for kk, vv in cbh_files.items():
    out_file = dom_dir / f"{kk}.nc"
    if out_file.exists():
        print(f"output file already exists, skipping: {out_file}")
        continue
    else:
        print(f"creating {out_file}")

    cbh_files_to_netcdf({kk: vv}, params, out_file)
