from typing import Literal

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
        self._initialize_data()
        self._construct_graph()

        return

    def _initialize_data(self):
        raise Exception("This must be overridden")
        return

    def _construct_graph(self):
        raise Exception("This must be overridden")
        return

    def insert_nodes(self, connectivity: list):
        pass

    def replace_nodes(self, replacements: list):
        pass

    # def freeze_graph(self):
    #     pass

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
        # self.freeze_graph()
        raise Exception("This must be overridden")
