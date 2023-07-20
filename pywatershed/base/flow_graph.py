from typing import Literal

from ..parameters import Parameters
from .accessor import Accessor
from .conservative_process import ConservativeProcess
from .control import Control


class FlowNode(Accessor):
    """FlowNode base class"""

    def __init__(self, control, calc_function):
        pass

    def calculate(inflow):
        # self._outflow =
        # self._storage_change =
        # self._storage =
        pass

    def advance():
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
        return

    def insert_notes(self, connectivity: list):
        pass

    def replace_nodes(self, replacements: list):
        pass

    def compute(self):
        pass


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
