from typing import Literal

import networkx as nx

from ..parameters import Parameters
from .accessor import Accessor
from .conservative_process import ConservativeProcess
from .control import Control


class FlowNode(Accessor):
    """
    A FlowNode base class
    Its calculate or calculate_subtimestep methods take an inflow and compute,
    at a minimum, outflow and storage change which are public quantities.
    Could also calculate head, storage, etc...
    """

    def calculate_subtimestep(self, inflow):
        raise Exception("This must be overridden")

    def advance(self):
        raise Exception("This must be overridden")

    @property
    def outflow(self):
        raise Exception("This must be overridden")

    @property
    def storage_change(self):
        raise Exception("This must be overridden")


class FlowGraph(ConservativeProcess):
    """
    A FlowGraph base class
    Has a list/connectivity of FlowNodes and a method to compute these in
    their sequence.
    Has methods to insert or replace nodes. The inserted nodes must
    calculate in a similar manner, taking the same inputs and providing the
    same outputs.

    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        metadata_patches: dict[dict] = None,
        metadata_patch_conflicts: Literal["ignore", "warn", "error"] = "error",
    ):
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
            metadata_patches=metadata_patches,
            metadata_patch_conflicts=metadata_patch_conflicts,
        )
        self.name = "FlowGraph"
        self._init_data()
        self._init_graph()
        self._graph_finalized = False

        return

    def _init_data(self):
        raise Exception("This must be overridden")
        return

    def _init_graph(self):
        raise Exception("This must be overridden")
        return

    def finalize_graph(self):
        self._graph_finalized = True
        raise Exception("This must be overridden")
        raise

    # def insert_nodes(self, connectivity: list):
    def insert_node(self, new_node_id, up_id, down_id):
        self._check_graph_editable()
        assert new_node_id not in list(self._graph.nodes)
        assert (up_id, down_id) in list(self._graph.edges)
        self.graph.remove_edge(up_id, down_id)
        nx.add_path(self._graph, [up_id, new_node_id, down_id])
        return

    def replace_nodes(self, replacements: list):
        self._check_graph_editable()
        pass

    def _check_graph_editable(self) -> None:
        if self._graph_finalized:
            raise RuntimeWarning("A finalized graph can not be edited.")
        else:
            return

    # def save_graph()
    # def load_graph()

    # ?
    # def get_outflow_mask(self):
    #     """This allows a global mass balance to be calculated"""
    #     return self._outflow_mask

    # ?
    # @property
    # def outflow_mask(self):
    #     return self._outflow_mask

    def _calculate_lateral_inflows(self):
        raise Exception("This must be overridden")

    def _sum_outflows_to_inflow(self, inode):
        raise Exception("This must be overridden")

    def _calculate(self, istep):
        if not self._graph_finalized:
            self.finalize_graph()
        raise Exception("This must be overridden")
