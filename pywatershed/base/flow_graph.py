from typing import Literal, Union

import networkx as nx
import numpy as np

from pywatershed.base.accessor import Accessor
from pywatershed.base.adapter import adaptable
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

    @property
    def outflow(self):
        raise Exception("This must be overridden")

    @property
    def storage_change(self):
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


class FlowGraph(Accessor):
    """
    Assembles and computes FlowNodes over a list/connectivity of FlowNodes
    and tracks their mass balance.
    The inputs inflows, node_maker_name, node_maker_index, and to_graph_index
    are collated. The order of execution of the graph is not the same as the
    supplied order, the execution order is solved from to_graph_index.
    Note that initial conditions are set by the node makers via their
    parameters.

    Args:
        control: a Control object
        inflows: An adaptable of the lateral inflows to the graph
        node_maker_dict: name:FlowNodeMaker
        node_maker_name: list or np.array of the maker name for each node
        node_maker_index: list or np.array of the index of the makerto
            use for each node
        to_graph_index: the index of the downstream index in the FlowGraph
            with -1 indicating an outflow node. This must specify a DAG.

    """

    def __init__(
        self,
        control: Control,
        # discretization: Parameters,  # could use this, but not necsesary
        # parameters: Parameters,  # unnecessary? use for to_graph_index
        inflows: adaptable,
        node_maker_dict: dict,
        node_maker_name: Union[list, np.ndarray],
        node_maker_index: Union[list, np.ndarray],
        to_graph_index: Union[list, np.ndarray],  # put in parameters?
        budget_type: Literal[None, "warn", "error"] = None,
    ):
        # add a budget like process adds a budget?
        # self._set_budget(basis="global")

        self.control = control
        self._inflows = inflows
        self._node_maker_dict = node_maker_dict
        self._node_maker_name = node_maker_name
        self._node_maker_index = node_maker_index
        self._to_graph_index = to_graph_index

        # basic checks
        self.nnodes = len(self._inflows.current)
        assert len(self._node_maker_name) == self.nnodes
        assert len(self._node_maker_index) == self.nnodes
        assert len(self._to_graph_index) == self.nnodes

        for fnm in self._node_maker_dict.values():
            assert isinstance(fnm, FlowNodeMaker)

        self._init_graph()
        self._init_variables()

        # private variables
        self._seg_upstream_inflow = np.zeros(self.nnodes) * zero
        self._seg_upstream_inflow_acc = np.zeros(self.nnodes) * zero
        self._seg_outflow_substep = np.zeros(self.nnodes) * zero

        return

    @staticmethod
    def get_dimensions() -> tuple:
        return ("nnodes",)

    @staticmethod
    def get_parameters() -> tuple:
        return ("to_graph_index",)

    @staticmethod
    def get_inputs() -> tuple:
        return ("inflow_vol",)

    @staticmethod
    def get_init_values() -> dict:
        return {
            "outflow_vol": nan,
            # "node_upstream_inflow": zero,
            "node_outflow": zero,
            "node_storage_change": zero,
        }

    def get_variables(cls) -> tuple:
        """Get a tuple of (public) variable names for this Process."""
        return list(cls.get_init_values().keys())

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": ["inflow_vol"],
            "outputs": ["outflow_vol"],
            "storage_changes": ["node_storage_change"],
        }

    def get_outflow_mask(self):
        return self._outflow_mask

    @property
    def outflow_mask(self):
        return self._outflow_mask

    def _init_graph(self) -> None:
        # where do flows exit the graph?
        self._outflow_mask = np.where(self._to_graph_index == -1, True, False)

        # which nodes do not have upstream nodes?
        self._headwater_nodes = set(range(self.nnodes)).difference(
            set(self._to_graph_index)
        )

        connectivity = []
        for inode in range(self.nnodes):
            tonode = self._to_graph_index[inode]
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
            zip(self._node_maker_name, self._node_maker_index)
        ):
            self._nodes += [
                self._node_maker_dict[maker_name].get_node(
                    self.control, maker_index
                )
            ]

    def _init_variables(self) -> None:
        """Initialize node variables tracked by the FlowGraph"""
        for key, val in self.get_init_values().items():
            self[key] = np.zeros(self.nnodes) + val

    def advance(self) -> None:
        # need to advance the inputs?
        self._inflows.advance()

        # advance all the nodes
        for node in self._nodes:
            node.advance()

        # no prognostic variables on the graph

        return

    def calculate(self, n_substeps=24) -> None:
        for node in self._nodes:
            node.prepare_timestep()

        self._seg_upstream_inflow_acc[:] = zero

        for istep in range(n_substeps):
            # This works because the first nodes calculated do
            # not have upstream reaches
            self._seg_upstream_inflow[:] = zero

            for jseg in self._node_order:
                # The first nodes calculated dont have upstream inflows
                # Eventually pass timestep length and n_substems to nodes
                # Calculate
                self._nodes[jseg].calculate_subtimestep(
                    istep,
                    self._seg_upstream_inflow[jseg],
                    self._inflows.current[jseg],
                )
                # Get the outflows back
                self._seg_outflow_substep[jseg] = self._nodes[
                    jseg
                ].outflow_substep
                # Add this node's outflow its downstream node's inflow
                self._sum_outflows_to_inflow(jseg)

            # <
            self._seg_upstream_inflow_acc += self._seg_upstream_inflow

        for node in self._nodes:
            node.finalize_timestep()

        # collect public variables from nodes
        # self.node_upstream_inflow[:] = (
        #     self._seg_upstream_inflow_acc / n_substeps
        # )

        for ii in range(self.nnodes):
            self.node_outflow[ii] = self._nodes[ii].outflow
            self.node_storage_change[ii] = self._nodes[ii].storage_change

        # global mass balance term
        s_per_time = self.control.time_step_seconds
        self.outflow_vol[:] = (
            np.where(self._outflow_mask, self.node_outflow, zero)
        ) * s_per_time

        return

    def _sum_outflows_to_inflow(self, jseg):
        if self._to_graph_index[jseg] >= 0:
            self._seg_upstream_inflow[
                self._to_graph_index[jseg]
            ] += self._seg_outflow_substep[jseg]

        return

    # dont have a netcdf output method

    def finalize(self) -> None:
        # eventually close netcdf
        return
