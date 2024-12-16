import pathlib as pl

import numpy as np
import pytest
import xarray as xr

from pywatershed.base.adapter import Adapter, AdapterNetcdf
from pywatershed.base.control import Control

# from pywatershed.base.flow_graph import FlowGraph
from pywatershed.constants import cm_to_cf, cms_to_cfs, nan, zero
from pywatershed.hydrology.starfit import StarfitFlowNodeMaker
from pywatershed.parameters import Parameters, StarfitParameters

# NB:
#   Here we are comparing a daily offline starfit against an hourly
#   StarfitNode. The reference output is the mean value from offline runs run
#   from 1995-2001 in the file
#   ../test_data/starfit/starfit_mean_output_1995-2001.nc
#   We only advance the hourly StarfitNode one substepper day. It's
#   resulting flow rates are identical but the change in storage is 1/24
#   of the daily value, so we check this. We have to track previous storage
#   to do this and get the delta storages.
#   TODO: There is no comparison of output files at the moment.
do_compare_output_files = True
do_compare_in_memory = True
rtol = atol = 1.0e-7

end_time = np.datetime64("2001-12-31 00:00:00")
# This is 117 reservoirs where
#    (parameters_ds.start_time == np.datetime64("1995-01-01 00:00:00"))
#    & (parameters_ds.end_time >= np.datetime64("2001-12-31 00:00:00"))
# fmt: off
starfit_inds_test = [
    0,   1,   2,   3,   4,   5,   6,   8,   9,   10,  11,  12,  13,
    15,  16,  18,  20,  21,  22,  23,  24,  25,  26,  28,  29,  30,
    31,  32,  33,  36,  37,  38,  40,  43,  44,  47,  48,  49,  51,
    52,  53,  55,  56,  59,  62,  63,  64,  65,  67,  68,  69,  70,
    71,  72,  74,  75,  76,  77,  86,  87,  89,  90,  91,  92,  93,
    94,  95,  96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106,
    107, 108, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120,
    122, 123, 130, 134, 137, 139, 140, 141, 145, 148, 149, 152, 154,
    155, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166
] # noqa
# fmt: on


@pytest.fixture(scope="function")
def parameters():
    parameter_file = pl.Path(
        "../test_data/starfit/starfit_original_parameters.nc"
    )
    parameters_ds = Parameters.from_netcdf(parameter_file).to_xr_ds()
    # wh_model_reservoirs = np.where(
    #     (parameters_ds.start_time == np.datetime64("1995-01-01 00:00:00"))
    #     & (parameters_ds.end_time >= np.datetime64("2001-12-31 00:00:00"))
    # )
    # asdf
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
        end_time,
        np.timedelta64(24, "h"),
    )
    control.options["budget_type"] = "error"
    return control


@pytest.fixture(scope="function")
def answers():
    ans_file = pl.Path("../test_data/starfit/starfit_mean_output_1995-2001.nc")
    ans = xr.open_dataset(ans_file)
    return ans.isel(grand_id=starfit_inds_test)


@pytest.mark.parametrize(
    "io_in_cfs", [True, False], ids=("io_in_cfs", "io_in_cms")
)
def test_starfit_flow_node_compare_starfit(
    control, parameters, answers, io_in_cfs, tmp_path
):
    if io_in_cfs:
        param_ds = parameters.to_xr_ds()
        param_ds["initial_storage"] *= cm_to_cf
        parameters = StarfitParameters.from_ds(param_ds)

    inflow_file = "../test_data/starfit/lake_inflow.nc"
    input_variables = AdapterNetcdf(inflow_file, "lake_inflow", control)
    nreservoirs = len(starfit_inds_test)

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
        io_in_cfs=io_in_cfs,
        nhrs_substep=24,
    )

    nodes = [
        node_maker.get_node(control, ii)
        for ii, zz in enumerate(starfit_inds_test)
    ]

    # we'll fill up timeseries arrays
    results = {
        "lake_storage": np.zeros([control.n_times, nreservoirs]) * np.nan,
        "lake_spill": np.zeros([control.n_times, nreservoirs]) * np.nan,
        "lake_release": np.zeros([control.n_times, nreservoirs]) * np.nan,
    }

    for istep in range(control.n_times):
        control.advance()
        inflows_node.advance()

        if io_in_cfs:
            inflows_node.current[:] *= cms_to_cfs

        for inode, node in enumerate(nodes):
            node.advance()
            node.prepare_timestep()
            for ss in range(1):
                node.calculate_subtimestep(
                    ss, inflows_node.current[inode], zero
                )
            node.finalize_timestep()

        # fill up timeseries arrays
        for ii, si in enumerate(starfit_inds_test):
            for var in results.keys():
                results[var][istep, ii] = nodes[ii][f"_{var}"][0]

    for var in results.keys():
        actual = results[var].mean(0)
        ans = answers[f"{var}_mean"].values
        if io_in_cfs:
            ans *= cms_to_cfs  # same for storage

        # <
        np.testing.assert_allclose(actual, ans, rtol=rtol, atol=atol)
