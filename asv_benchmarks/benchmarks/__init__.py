import os
import pathlib as pl
import sys

print("sys.version: ", sys.version)

# TODO remove backwards compatiability with pynhm once
#      reported slows downs are sorted output
# TODO remove backwards compatability with <0.2.0 once
#      it is released.

try:
    # backwards compatability
    import pynhm as pws

    _is_pws = False
except ModuleNotFoundError:
    import pywatershed as pws

    _is_pws = True

if not "constants" in pws.__dict__.keys():
    del pws
    import pywatershed as pws

    _is_pws = False

# The package is installed without data. The test data are relative to the repo
# used to do the install. So use the asv env var to find that.
try:
    asv_conf_dir = pl.Path(os.environ["ASV_CONF_DIR"]).resolve()
except KeyError:
    asv_conf_dir = pl.Path(".")


assert asv_conf_dir.exists()
pws_root = asv_conf_dir / "../pywatershed"
test_data_dir = asv_conf_dir / "../test_data"


def parameterized(names, params):
    def decorator(func):
        func.param_names = names
        func.params = params
        return func

    return decorator
