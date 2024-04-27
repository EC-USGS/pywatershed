import pickle
import random

import numpy as np
import pytest
import xarray as xr

from pywatershed import Starfit
from pywatershed.base.adapter import Adapter, AdapterNetcdf
from pywatershed.base.control import Control

# from pywatershed.base.flow_graph import FlowGraph
from pywatershed.base.parameters import Parameters
from pywatershed.constants import __pywatershed_root__ as pws_root
from pywatershed.constants import nan, zero
from pywatershed.hydrology.starfit import StarfitFlowNodeMaker
from pywatershed.parameters import StarfitParameters
from utils_compare import compare_in_memory, compare_netcdfs

# NB: THere is no real comparison of output files because the answer files
#     have different units. Could create a class to manage this but
#     I've checked that the memory values match the file values.
#     do_compare_output_files write files but does not check them
do_compare_output_files = True
do_compare_in_memory = True
rtol = atol = 1.0e-7


starfit_inds_test = random.sample(range(167), 7)
# starfit_inds_test = list(range(167))  # this ran previously, takes ~4min

data_dir = pws_root / "hydrology/starfit_minimal"


@pytest.fixture(scope="function")
def parameters():
    # TODO: make this work with the original source files
    param_files = {
        "grand_ids": data_dir / "domain_grand_ids.nc",
        "resops_domain": data_dir / "ResOpsUS.nc",
        "istarf_conus": data_dir / "starfit/ISTARF-CONUS.nc",
        "grand_dams": data_dir / "grand_dams.nc",
    }

    # passing the parameter names is optional
    starfit_parameters = Starfit.get_parameters()
    parameters_ds = StarfitParameters.from_netcdf(
        **param_files,
        param_names=starfit_parameters,
    ).to_xr_ds()

    # the first parameter/reservoir
    merge_list = []
    for ii in starfit_inds_test:
        merge_list += [parameters_ds.isel(nreservoirs=slice(ii, ii + 1))]
    parameters_ds = xr.concat(merge_list, dim="nreservoirs")
    parameters = StarfitParameters.from_ds(parameters_ds)

    return parameters


@pytest.fixture(scope="function")
def control(parameters):
    control = Control(
        parameters.variables["start_time"].min(),
        np.datetime64("2019-09-30 00:00:00"),
        np.timedelta64(24, "h"),
    )
    control.options["budget_type"] = "warn"
    return control


# @pytest.mark.parametrize("calc_method", calc_methods)
def test_starfit_flow_node_compare_starfit(control, parameters, tmp_path):
    with open(data_dir / "starfit_outputs.pickle", "rb") as handle:
        ans_dict = pickle.load(handle)
    answers = {}
    answers["lake_storage"] = ans_dict["Ssim_list"]
    answers["lake_release"] = ans_dict["Rsim_list"]
    answers["lake_spill"] = ans_dict["SPILLsim_list"]
    ans_grand_ids = ans_dict["grand_ids"]

    starfit_inflow_name = Starfit.get_inputs()[0]
    nc_path = data_dir / f"{starfit_inflow_name}.nc"
    input_variables = AdapterNetcdf(nc_path, starfit_inflow_name, control)

    class NodeInflowAdapter(Adapter):
        def __init__(
            self,
            starfit_inflows: Adapter,
            variable: str = "inflows",
        ):
            self._variable = variable
            self._starfit_inflows = starfit_inflows

            self._nnodes = len(starfit_inds_test)
            self._current_value = np.zeros(self._nnodes) * nan
            return

        def advance(self) -> None:
            self._starfit_inflows.advance()
            self._current_value[:] = self._starfit_inflows.current[
                tuple(starfit_inds_test),
            ]
            return

    inflows_node = NodeInflowAdapter(input_variables)

    node_maker = StarfitFlowNodeMaker(
        discretization=None,
        parameters=parameters,
        # calc_method=calc_method
    )

    nodes = [
        node_maker.get_node(control, ii)
        for ii, zz in enumerate(starfit_inds_test)
    ]

    for istep in range(control.n_times):
        control.advance()
        inflows_node.advance()
        for inode, node in enumerate(nodes):
            node.advance()
            node.prepare_timestep()
            node.calculate_subtimestep(0, inflows_node.current[inode], zero)
            # subtimesteps dont matter to starfit node.
            node.calculate_subtimestep(1, inflows_node.current[inode], zero)

        # check as we go
        ymd = control.current_datetime.strftime("%Y-%m-%d")
        for ii, si in enumerate(starfit_inds_test):
            assert parameters.parameters["grand_id"][ii] == ans_grand_ids[si]
            for var in ["lake_storage", "lake_release", "lake_spill"]:
                actual = nodes[ii][f"_{var}"]
                ans_series = answers[var][si]
                if ymd in ans_series.index:
                    np.testing.assert_allclose(
                        actual, ans_series[ymd], rtol=rtol, atol=atol
                    )
                else:
                    assert np.isnan(actual)
