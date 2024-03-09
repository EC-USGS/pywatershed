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
    AdapterExchangeHruSegment,
    PRMSChannelFlowNodeMaker,
)
from pywatershed.parameters import PrmsParameters
from utils_compare import compare_in_memory, compare_netcdfs

do_compare_output_files = False
do_compare_in_memory = True
rtol = atol = 1.0e-7

# test the exchange too?
rename_vars = {
    "channel_outflow_vol": "outflow_vol",
    # "seg_upstream_inflow": "node_upstream_inflow",
    "seg_outflow": "node_outflow",
    "seg_stor_change": "node_storage_change",
}


# we only need to test this on one domain.


@pytest.fixture(scope="function")
def control(simulation):
    if "drb:nhm" not in simulation["name"]:
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
    # PassThroughNode

    node_maker_dict = {
        "prms_channel": PRMSChannelFlowNodeMaker(
            discretization_prms, parameters_prms
        ),
        "pass_throughs": PassThroughNodeMaker(),
    }

    # flow exchange: combine PRMS lateral inflows to a single volumetric inflow
    output_dir = simulation["output_dir"]
    input_variables = {}
    for key in PRMSChannel.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = AdapterNetcdf(nc_path, key, control)

    inflow_exchange_prms = AdapterExchangeHruSegment(
        "inflow_vol", parameters_prms, **input_variables
    )

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

    inflow_exchange = GraphInflowAdapter("inflow_vol", inflow_exchange_prms)

    # FlowGraph
    nnodes = parameters_prms.dims["nsegment"] + 1
    node_maker_name = ["prms_channel"] * nnodes
    node_maker_name[-1] = "pass_throughs"
    node_maker_index = np.arange(nnodes)
    node_maker_index[-1] = 0
    to_graph_index = np.zeros(nnodes, dtype=np.int64)
    to_graph_index[0:-1] = discretization_prms.parameters["tosegment"] - 1
    intervene_above_node = 34  # random choice
    wh_intervene = np.where(to_graph_index == intervene_above_node)
    to_graph_index[-1] = intervene_above_node
    to_graph_index[wh_intervene] = nnodes - 1

    flow_graph = FlowGraph(
        control,
        inflow_exchange,
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
        print(istep)
        control.advance()
        flow_graph.advance()
        flow_graph.calculate()
        # flow_graph.output()

        # check exchange
        lateral_inflow_answers.advance()
        assert (
            inflow_exchange.current[0:-1] == lateral_inflow_answers.current
        ).all()

        answers_on_graph = {}
        for key, val in answers.items():
            val.advance()
            answers_on_graph[key] = np.zeros(len(val.current) + 1)
            answers_on_graph[key][0:-1] = val.current
            # just use this as correct and verify the rest of the graph is
            # unchanged
            answers_on_graph[key][-1] = flow_graph[key][-1]
            print(key, flow_graph[key][-1])

        if do_compare_in_memory:
            compare_in_memory(
                flow_graph,
                answers_on_graph,
                atol=atol,
                rtol=rtol,
                fail_after_all_vars=False,
            )

    flow_graph.finalize()
