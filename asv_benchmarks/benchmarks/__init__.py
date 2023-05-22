import os
import pathlib as pl

import pynhm as pws

# For backwards compatabilty with pynhm
if "constants" in pws.__dict__.keys():
    _is_pws = False
else:
    import pywatershed as pws

    _is_pws = True


# The package is installed without data. The test data are relative to the repo
# used to do the install. So use the asv env var to find that.
asv_conf_dir = pl.Path(os.environ["ASV_CONF_DIR"]).resolve()
assert asv_conf_dir.exists()
pws_root = asv_conf_dir / "../pywatershed"
test_data_dir = asv_conf_dir / "../test_data"


def parameterized(names, params):
    def decorator(func):
        func.param_names = names
        func.params = params
        return func

    return decorator
