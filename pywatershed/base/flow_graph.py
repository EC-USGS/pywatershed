from typing import Literal

import networkx as nx
import numpy as np

from pywatershed.base.accessor import Accessor
from pywatershed.base.adapter import adaptable
from pywatershed.base.conservative_process import ConservativeProcess
from pywatershed.base.control import Control
from pywatershed.constants import nan, zero
from pywatershed.parameters import Parameters


class FlowNode(Accessor):
    """
    A FlowNode base class
    Its calculate or calculate_subtimestep methods take an inflow and compute,
    at a minimum, outflow and storage change which are public quantities.
    Could also calculate head, storage, etc...
    """

    def __init__(self, control):
        raise Exception("This must be overridden")

    def prepare_timestep(self):
        raise Exception("This must be overridden")

    def calculate_subtimestep(self, isubstep, inflow_upstream, inflow_lateral):
        raise Exception("This must be overridden")

    def advance(self):
        raise Exception("This must be overridden")

    def finalize_timestep(self):
        raise Exception("This must be overridden")

    @property
    def outflow(self):
        raise Exception("This must be overridden")

    @property
    def storage_change(self):
        raise Exception("This must be overridden")

    @property
    def sink_source(self):
        raise Exception("This must be overridden")


class FlowNodeMaker(Accessor):
    """
    FlowNodeMaker creates FlowNodes of a certain type.
    * take inputs...


    Notes:
      * initialize the the data needed to return the nodes
      * method to give a node at index
      * the node maker is only used on the init of the FlowGraph
      * does not track mass balance

    """

    def __init__(
        self,
        discretization: Parameters = None,
        parameters: Parameters = None,
    ):
        self.name = "FlowNodeMaker"
        return

    def get_node(control, index):
        raise Exception("This must be overridden")


class FlowGraph(ConservativeProcess):
    """
    Assembles and computes FlowNodes over a list/connectivity of FlowNodes
    and tracks their mass balance.
    The inputs inflows, node_maker_name, node_maker_index, and to_graph_index
    are collated. The order of execution of the graph is not the same as the
    supplied order, the execution order is solved from to_graph_index.
    Note that initial conditions are set by the node makers via their
    parameters.

    Parameters:
        node_maker_name: list or np.array of the maker name for each node
        node_maker_index: list or np.array of the index of the makerto
            use for each node
        to_graph_index: the index of the downstream index in the FlowGraph
            with -1 indicating an outflow node. This must specify a DAG.

    Args:
        control: a Control object
        inflows: An adaptable of inflows to the graph, often referred to as
            "lateral" flows.
        node_maker_dict: {name: flow_node_maker_instance, ...}
        budget_type: If we use a budget, what action to take on budget
            imbalance.
    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,  # could use this, but not necsesary
        parameters: Parameters,  # unnecessary? use for to_graph_index
        inflows: adaptable,
        # node_maker_dict is not really a parameter, it's a composition
        # mechanism
        node_maker_dict: dict,
        budget_type: Literal[None, "warn", "error"] = None,
        verbose: bool = None,
    ):
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "FlowGraph"

        self._set_inputs(locals())
        self._set_options(locals())

        for fnm in self._node_maker_dict.values():
            assert isinstance(fnm, FlowNodeMaker)

        self._init_graph()

        # If/when FlowGraph handles nodes which dont tautologically balance
        # could allow the basis to be unit.
        self._set_budget(basis="global")

        return

    @staticmethod
    def get_dimensions() -> tuple:
        return ("nnodes",)

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "node_maker_name",
            "node_maker_index",
            "to_graph_index",
        )

    @staticmethod
    def get_inputs() -> tuple:
        return ("inflows",)

    @staticmethod
    def get_init_values() -> dict:
        return {
            "outflows": nan,
            "node_upstream_inflows": nan,
            "node_outflows": nan,
            "node_storage_changes": nan,
            "node_storages": nan,
            "node_sink_source": nan,
        }

    @classmethod
    def get_variables(cls) -> tuple:
        """Get a tuple of (public) variable names for this Process."""
        return list(cls.get_init_values().keys())

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "inflows",
            ],
            "outputs": ["outflows"],
            "storage_changes": ["node_storage_changes", "_neg_sink_source"],
        }

    def get_outflow_mask(self):
        return self._outflow_mask

    @property
    def outflow_mask(self):
        return self._outflow_mask

    def _set_initial_conditions(self) -> None:
        self._node_upstream_inflow_sub = np.zeros(self.nnodes) * nan
        self._node_upstream_inflow_acc = np.zeros(self.nnodes) * nan
        self._node_outflow_substep = np.zeros(self.nnodes) * nan
        self._neg_sink_source = np.zeros(self.nnodes) * nan
        return

    def _init_graph(self) -> None:
        params = self._params.parameters
        # where do flows exit the graph?
        self._outflow_mask = np.where(
            params["to_graph_index"] == -1, True, False
        )

        # which nodes do not have upstream nodes?
        self._headwater_nodes = set(range(self.nnodes)).difference(
            set(params["to_graph_index"])
        )

        connectivity = []
        for inode in range(self.nnodes):
            tonode = params["to_graph_index"][inode]
            if tonode < 0:
                continue
            connectivity.append(
                (
                    inode,
                    tonode,
                )
            )

        # use networkx to calculate the Directed Acyclic Graph
        self._graph = nx.DiGraph()
        self._graph.add_edges_from(connectivity)

        # Make sure the user isnt suppling disconnected nodes
        assert len(self._graph) == self.nnodes, "Disconnected nodes present."
        # should we allow disconnected nodes? seems not
        # these are handled in PRMSChannel if we change our minds

        if self.nnodes > 1:
            node_order = list(nx.topological_sort(self._graph))
        else:
            node_order = [0]

        self._node_order = np.array(node_order, dtype="int64")

        # any performance for doing a hash table up front?
        # a hash {to_seg: [from_seg_0, ..., from_seg_n]}

        # instatiate the nodes
        self._nodes = []
        for ii, (maker_name, maker_index) in enumerate(
            zip(params["node_maker_name"], params["node_maker_index"])
        ):
            self._nodes += [
                self._node_maker_dict[maker_name].get_node(
                    self.control, maker_index
                )
            ]

    def _advance_variables(self) -> None:
        for node in self._nodes:
            node.advance()

        # no prognostic variables on the graph
        return

    def calculate(self, time_length: float, n_substeps: int = 24) -> None:
        params = self._params.parameters

        for node in self._nodes:
            node.prepare_timestep()

        self._node_upstream_inflow_acc[:] = zero

        for istep in range(n_substeps):
            # This works because the first nodes calculated do
            # not have upstream reaches
            self._node_upstream_inflow_sub[:] = zero

            for inode in self._node_order:
                # The first nodes calculated dont have upstream inflows
                # Eventually pass timestep length and n_substems to nodes
                # Calculate
                self._nodes[inode].calculate_subtimestep(
                    istep,
                    self._node_upstream_inflow_sub[inode],
                    self.inflows[inode],
                )
                # Get the outflows back
                self._node_outflow_substep[inode] = self._nodes[
                    inode
                ].outflow_substep
                # Add this node's outflow its downstream node's inflow
                if params["to_graph_index"][inode] >= 0:
                    self._node_upstream_inflow_sub[
                        params["to_graph_index"][inode]
                    ] += self._node_outflow_substep[inode]

            # <
            # not sure how PRMS-specific this is
            self._node_upstream_inflow_acc += self._node_upstream_inflow_sub

        for node in self._nodes:
            node.finalize_timestep()

        self.node_upstream_inflows[:] = (
            self._node_upstream_inflow_acc / n_substeps
        )

        for ii in range(self.nnodes):
            self.node_outflows[ii] = self._nodes[ii].outflow
            self.node_storage_changes[ii] = self._nodes[ii].storage_change
            self.node_storages[ii] = self._nodes[ii].storage
            self.node_sink_source[ii] = self._nodes[ii].sink_source

        self._neg_sink_source[:] = -1 * self.node_sink_source

        # global mass balance term
        self.outflows[:] = np.where(
            self._outflow_mask, self.node_outflows, zero
        )

        if self.budget is not None:
            self.budget.advance()
            self.budget.calculate()

        return


def inflow_exchange_factory(
    dimension_names: tuple,
    parameter_names: tuple,
    input_names: tuple,
    init_values: dict,
    mass_budget_terms: dict,
    calculation: callable,
):
    class InflowExchange(ConservativeProcess):
        def __init__(
            self,
            control: Control,
            discretization: Parameters,
            parameters: Parameters,
            budget_type: Literal[None, "warn", "error"] = None,
            verbose: bool = None,
            budget_basis="global",
            **kwargs,
        ) -> None:
            super().__init__(
                control=control,
                discretization=discretization,
                parameters=parameters,
            )
            self.name = "HruSegmentFlowExchange"

            self._set_inputs(locals() | kwargs)
            self._set_options(locals())

            self._set_budget(basis=budget_basis)

            return

        @staticmethod
        def get_dimensions() -> tuple:
            return dimension_names

        @staticmethod
        def get_parameters() -> tuple:
            return parameter_names

        @staticmethod
        def get_inputs() -> tuple:
            return input_names

        @staticmethod
        def get_init_values() -> dict:
            return init_values

        @staticmethod
        def get_mass_budget_terms():
            return mass_budget_terms

        def _set_initial_conditions(self) -> None:
            pass

        def _advance_variables(self) -> None:
            pass

        def _calculate(self, simulation_time: float) -> None:
            calculation(self)
            return

    return InflowExchange
