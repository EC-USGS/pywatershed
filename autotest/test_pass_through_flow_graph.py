import pathlib as pl

import numpy as np
import pytest

from pywatershed import PRMSChannel
from pywatershed.base.adapter import Adapter, AdapterNetcdf, adapter_factory
from pywatershed.base.control import Control
from pywatershed.base.flow_graph import FlowGraph
from pywatershed.base.parameters import Parameters
from pywatershed.constants import nan, zero
from pywatershed.hydrology.pass_through_node import PassThroughNodeMaker
from pywatershed.hydrology.prms_channel_flow_graph import (
    HruSegmentInflowAdapter,
    PRMSChannelFlowNodeMaker,
)
from pywatershed.parameters import PrmsParameters
from utils_compare import compare_in_memory, compare_netcdfs

do_compare_output_files = False
do_compare_in_memory = True
rtol = atol = 1.0e-7

# test the exchange too?
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
        pytest.skip("Only testing passthrough flow graph for drb_2yr:nhm")
    return Control.load(simulation["control_file"], warn_unused_options=False)


@pytest.fixture(scope="function")
def discretization_prms(simulation):
    dis_hru_file = simulation["dir"] / "parameters_dis_hru.nc"
    dis_seg_file = simulation["dir"] / "parameters_dis_seg.nc"
    return Parameters.merge(
        Parameters.from_netcdf(dis_hru_file, encoding=False),
        Parameters.from_netcdf(dis_seg_file, encoding=False),
    )


@pytest.fixture(scope="function")
def parameters_prms(simulation, control):
    param_file = simulation["dir"] / "parameters_PRMSChannel.nc"
    return PrmsParameters.from_netcdf(param_file)


def test_prms_channel_pass_through_compare_prms(
    simulation,
    control,
    discretization_prms,
    parameters_prms,
    tmp_path,
):
    # combine PRMS lateral inflows to a single non-volumetric inflow
    output_dir = simulation["output_dir"]
    input_variables = {}
    for key in PRMSChannel.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = AdapterNetcdf(nc_path, key, control)

    inflows_prms = HruSegmentInflowAdapter(parameters_prms, **input_variables)

    class GraphInflowAdapter(Adapter):
        def __init__(
            self,
            variable: str,
            prms_inflows: Adapter,
        ):
            self._variable = variable
            self._prms_inflows = prms_inflows

            self._nnodes = len(self._prms_inflows.current) + 1
            self._current_value = np.zeros(self._nnodes) * nan
            return

        def advance(self) -> None:
            self._prms_inflows.advance()
            self._current_value[0:-1] = self._prms_inflows.current
            self._current_value[-1] = zero  # no inflow at the pass through
            return

    inflows_graph = GraphInflowAdapter("inflow_vol", inflows_prms)

    # FlowGraph

    node_maker_dict = {
        "prms_channel": PRMSChannelFlowNodeMaker(
            discretization_prms, parameters_prms
        ),
        "pass_throughs": PassThroughNodeMaker(),
    }
    nnodes = parameters_prms.dims["nsegment"] + 1
    node_maker_name = ["prms_channel"] * nnodes
    node_maker_name[-1] = "pass_throughs"
    node_maker_index = np.arange(nnodes)
    node_maker_index[-1] = 0
    to_graph_index = np.zeros(nnodes, dtype=np.int64)
    dis_params = discretization_prms.parameters
    to_graph_index[0:-1] = dis_params["tosegment"] - 1
    nhm_seg_intervene_above = 1829  # 1829
    wh_intervene_above_nhm = np.where(
        dis_params["nhm_seg"] == nhm_seg_intervene_above
    )
    wh_intervene_below_nhm = np.where(
        (dis_params["tosegment"] - 1) == wh_intervene_above_nhm[0][0]
    )
    # have to map to the graph from an index found in prms_channel
    wh_intervene_above_graph = np.where(
        (np.array(node_maker_name) == "prms_channel")
        & (node_maker_index == wh_intervene_above_nhm[0][0])
    )
    wh_intervene_below_graph = np.where(
        (np.array(node_maker_name) == "prms_channel")
        & np.isin(node_maker_index, wh_intervene_below_nhm)
    )

    to_graph_index[-1] = wh_intervene_above_graph[0][0]
    to_graph_index[wh_intervene_below_graph] = nnodes - 1

    flow_graph = FlowGraph(
        control,
        inflows_graph,
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
        # flow_graph.output()

        # check exchange
        lateral_inflow_answers.advance()
        assert (
            inflows_graph.current[0:-1] == lateral_inflow_answers.current
        ).all()

        answers_from_graph = {}
        for key, val in answers.items():
            val.advance()
            answers_from_graph[key] = np.zeros(len(val.current) + 1)
            answers_from_graph[key][0:-1] = val.current
            # just use this as correct and verify the rest of the graph is
            # unchanged
            answers_from_graph[key][-1] = flow_graph[key][-1]

        if do_compare_in_memory:
            answers_conv_vol = {}
            for key, val in answers_from_graph.items():
                if key in convert_to_vol:
                    answers_conv_vol[key] = val / (24 * 60 * 60)
                else:
                    answers_conv_vol[key] = val

            compare_in_memory(
                flow_graph,
                answers_conv_vol,
                atol=atol,
                rtol=rtol,
                fail_after_all_vars=False,
            )

    flow_graph.finalize()
