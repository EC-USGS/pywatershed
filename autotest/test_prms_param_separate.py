import pathlib as pl

import numpy as np
import pytest

import pywatershed
from pywatershed.base.data_model import open_datasetdict
from pywatershed.parameters import PrmsParameters
from pywatershed.utils import separate_domain_params_dis_to_ncdf

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


@pytest.mark.parametrize("use_xr", [True])  # TODO: add False
def test_param_sep(domain, params, use_xr, tmp_path):
    tmp_path = pl.Path(tmp_path)

    domain_name = domain["domain_name"]
    prms_param_file = domain["param_file"]

    print(tmp_path)
    proc_nc_files = separate_domain_params_dis_to_ncdf(
        prms_param_file,
        domain_name,
        tmp_path,
        process_list=nhm_processes,
        use_xr=use_xr,
    )
    assert len(proc_nc_files) == len(nhm_processes) + 3  # hru, seg, both

    # check roundtrip: to file and back
    for proc_name, proc_file in proc_nc_files.items():
        print(proc_file)
        params_file = open_datasetdict(proc_file)

        if isinstance(proc_name, type):
            # some coords in the file may not be in the class' parameters. In
            # that case, assert these coords are in file coords
            coords_not_in_params = set(params_file.variables.keys()).difference(
                set(proc_name.get_parameters())
            )
            assert not len(
                coords_not_in_params.difference(set(params_file["coords"].keys()))
            )

        # checking dim equality is not necessary, it is covered by the
        # metadata check
        # loop on what's in the file and not the process parameters because
        # some process parameters come from dis now. The sufficiency of
        # the parameters on file will come from running processes with these
        # files
        for param in params_file.variables.keys():
            np.testing.assert_equal(
                params_file.variables[param],
                params.parameters[param],
            )
            np.testing.assert_equal(
                params_file["metadata"][param]["dims"],
                params.metadata[param]["dims"],
            )

    return
