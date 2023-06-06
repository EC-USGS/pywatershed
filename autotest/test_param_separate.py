import pathlib as pl

import numpy as np
import pytest
import xarray as xr

import pywatershed
from pywatershed.base.control import Control
from pywatershed.base.data_model import xr_ds_to_dd
from pywatershed.parameters import PrmsParameters
from pywatershed.utils import separate_domain_params_to_ncdf


nhm_processes = [
    pywatershed.PRMSSolarGeometry,
    pywatershed.PRMSAtmosphere,
    pywatershed.PRMSCanopy,
    pywatershed.PRMSSnow,
    pywatershed.PRMSRunoff,
    pywatershed.PRMSSoilzone,
    pywatershed.PRMSGroundwater,
    pywatershed.PRMSChannel,
]


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


@pytest.mark.parametrize("use_xr", [True])  # TODO: add False
def test_param_sep(domain, control, use_xr, tmp_path):
    tmp_path = pl.Path(tmp_path)

    domain_name = domain["domain_name"]
    prms_param_file = domain["param_file"]

    print(tmp_path)
    proc_nc_files = separate_domain_params_to_ncdf(
        prms_param_file,
        domain_name,
        tmp_path,
        process_list=nhm_processes,
        use_xr=use_xr,
    )
    assert len(proc_nc_files) == len(nhm_processes)

    # check roundtrip to file and back
    for proc_class, proc_file in proc_nc_files.items():
        params_file = xr_ds_to_dd(xr.open_dataset(proc_file))
        assert set(params_file["data_vars"].keys()) == set(
            proc_class.get_parameters()
        )

        # checking dim equality is not necessary, it is covered by the
        # metadata check
        for param in proc_class.get_parameters():
            np.testing.assert_equal(
                params_file["data_vars"][param],
                control.params.parameters[param],
            )
            np.testing.assert_equal(
                params_file["metadata"][param]["dims"],
                control.params.metadata[param]["dims"],
            )

    return
