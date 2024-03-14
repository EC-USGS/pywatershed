import numpy as np
import pytest

from pywatershed import PRMSChannel
from pywatershed.base.adapter import AdapterNetcdf, adapter_factory
from pywatershed.base.control import Control
from pywatershed.base.flow_graph import FlowGraph
from pywatershed.base.parameters import Parameters
from pywatershed.hydrology.prms_channel_flow_graph import (
    HruSegmentInflowAdapter,
    PRMSChannelFlowNodeMaker,
)
from pywatershed.parameters import PrmsParameters
from utils_compare import compare_in_memory, compare_netcdfs

# NB: THere is no real comparison of output files because the answer files
#     have different units. Could create a class to manage this but
#     I've checked that the memory values match the file values.
#     do_compare_output_files write files but does not check them
do_compare_output_files = True
do_compare_in_memory = True
rtol = atol = 1.0e-7

calc_methods = ("numpy", "numba")

rename_vars = {
    "channel_outflow_vol": "outflows",
    "seg_upstream_inflow": "node_upstream_inflows",
    "seg_outflow": "node_outflows",
    "seg_stor_change": "node_storage_changes",
}

# the above values in rename_vars need converted to volume in these cases
convert_to_vol = ["outflows", "node_storage_changes"]


@pytest.fixture(scope="function")
def control(simulation):
    if "drb_2yr:nhm" not in simulation["name"]:
        pytest.skip("Only testing prms channel flow graph for drb_2yr:nhm")
    control = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    del control.options["netcdf_output_dir"]
    del control.options["netcdf_output_var_names"]
    return control


@pytest.fixture(scope="function")
def discretization(simulation):
    dis_hru_file = simulation["dir"] / "parameters_dis_hru.nc"
    dis_seg_file = simulation["dir"] / "parameters_dis_seg.nc"
    return Parameters.merge(
        Parameters.from_netcdf(dis_hru_file, encoding=False),
        Parameters.from_netcdf(dis_seg_file, encoding=False),
    )


@pytest.fixture(scope="function")
def parameters(simulation, control):
    param_file = simulation["dir"] / "parameters_PRMSChannel.nc"
    return PrmsParameters.from_netcdf(param_file)


@pytest.mark.parametrize("calc_method", calc_methods)
def test_prms_channel_flow_graph_compare_prms(
    simulation, control, discretization, parameters, tmp_path, calc_method
):
    node_maker_dict = {
        "prms_channel": PRMSChannelFlowNodeMaker(
            discretization, parameters, calc_method=calc_method
        ),
    }

    # combine lateral inflows to a single non-volumetric inflow
    output_dir = simulation["output_dir"]
    input_variables = {}
    for key in PRMSChannel.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = AdapterNetcdf(nc_path, key, control)

    inflow_prms = HruSegmentInflowAdapter(parameters, **input_variables)

    # FlowGraph
    nsegment = parameters.dims["nsegment"]
    node_maker_name = ["prms_channel"] * nsegment
    node_maker_index = np.arange(nsegment)
    to_graph_index = discretization.parameters["tosegment"] - 1

    flow_graph = FlowGraph(
        control,
        inflow_prms,
        node_maker_dict,
        node_maker_name,
        node_maker_index,
        to_graph_index,
        budget_type="error",
    )

    if do_compare_output_files:
        nc_parent = tmp_path / simulation["name"].replace(":", "_")
        flow_graph.initialize_netcdf(nc_parent)

    if do_compare_in_memory:
        answers = {}
        for var in PRMSChannel.get_variables():
            if var not in rename_vars.keys():
                continue
            new_var = rename_vars[var]
            var_pth = output_dir / f"{var}.nc"
            answers[new_var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

        # answers for the exchange
        var = "seg_lateral_inflow"
        var_pth = output_dir / f"{var}.nc"
        lateral_inflow_answers = adapter_factory(
            var_pth, variable_name=var, control=control
        )

    for istep in range(control.n_times):
        control.advance()
        flow_graph.advance()
        flow_graph.calculate()
        flow_graph.output()

        # check exchange
        lateral_inflow_answers.advance()
        np.testing.assert_allclose(
            inflow_prms.current,
            lateral_inflow_answers.current,
            rtol=rtol,
            atol=atol,
        )

        for var in answers.values():
            var.advance()

        if do_compare_in_memory:
            answers_conv_vol = {}
            for key, val in answers.items():
                if key in convert_to_vol:
                    answers_conv_vol[key] = val.current / (24 * 60 * 60)
                else:
                    answers_conv_vol[key] = val.current

            compare_in_memory(
                flow_graph,
                answers_conv_vol,
                atol=atol,
                rtol=rtol,
                skip_missing_ans=False,
                fail_after_all_vars=False,
            )

    flow_graph.finalize()
