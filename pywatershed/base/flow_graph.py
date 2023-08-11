from typing import Literal

from ..parameters import Parameters
from .accessor import Accessor
from .conservative_process import ConservativeProcess
from .control import Control


class FlowNode(Accessor):
    """FlowNode base class"""

    def __init__(self, control, calc_function):
        pass

    def calculate(self, inflow):
        # self._outflow =
        # self._storage_change =
        # self._storage =
        pass

    def advance(self):
        pass

    @property
    def outflow(self):
        return self._outflow

    @property
    def storage(self):
        return self._storage

    @property
    def storage_change(self):
        return self._storage_change


class FlowGraph(ConservativeProcess):
    """FlowDAG base class"""

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

        # self._construct_graph()

        return

    def _construct_graph(self):
        # self._to_segment =  # int64
        # self._outflow_mask
        # self._segment_order
        # self.flow_nodes objects
        return

    def insert_nodes(self, connectivity: list):
        pass

    def replace_nodes(self, replacements: list):
        pass

    def calculate_graph(self, istep):
        # rename eventually

        # freeze graph upon calculate

        for step in self.substeps:
            self._zero_upstream_inflows()

            for seg in self.segments:
                outflow = self.flow_nodes[seg].calculate(self.inflows[seg])
                self._sum_outflows_to_inflow(seg, outflow)

        self._finalize_time_step()

    def _sum_outflows_to_inflow(self, seg, outflow):
        if self._tosegment[seg] >= 0:
            self.seg_upstream_inflow[self._tosegment[seg]] += self._outflow_ts[
                seg
            ]

    def _zero_upstream_inflows():
        self.inflows[:] = zero
        return

    # def save()
    # def load()


"""
A FlowNode
takes inflow, has a function to compute storage change and outflow


A FlowGraph
Has a connectivity of FlowNodes
channel is a FlowGraph

has method to compute graph

has methods to insert or replace nodes.
river_lake_network = prms_channel.insert_nodes(starfit, connectivity)
river_lake_network = prms_channel.replace_nodes(starfit, replacements)
"""
