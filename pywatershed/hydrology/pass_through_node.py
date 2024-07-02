from ..base.control import Control
from ..base.flow_graph import FlowNode, FlowNodeMaker
from ..constants import nan, zero


class PassThroughNode(FlowNode):
    """A FlowNode instance that gives what it takes and dosent store.

    See :class:`FlowGraph` for a worked example using PassThroughNode.
    """

    def __init__(self, control: Control):
        """Initialize a PassThroughNode.

        Args:
            control: A control object.
        """
        self.control = control
        return

    def prepare_timestep(self):
        # self._simulation_time = simulation_time  # add argument?
        self._accum_inflow = zero
        return

    def calculate_subtimestep(
        self,
        isubstep: int,
        inflow_upstream: float,
        inflow_lateral: float,
    ):
        self._inflow_subtimestep = inflow_upstream + inflow_lateral
        self._accum_inflow += self._inflow_subtimestep
        self._seg_outflow = self._accum_inflow / (isubstep + 1)
        return

    def finalize_timestep(self):
        return

    def advance(self):
        return

    @property
    def outflow(self):
        return self._seg_outflow

    @property
    def outflow_substep(self):
        return self._inflow_subtimestep

    @property
    def storage_change(self):
        return zero

    @property
    def storage(self):
        return nan

    @property
    def sink_source(self):
        return zero


class PassThroughNodeMaker(FlowNodeMaker):
    """A FlowNodeMaker of PassThroughNodes.

    See :class:`FlowGraph` for a worked example using PassThroughNode.
    """

    def __init__(
        self,
    ) -> None:
        """Initialize a PassThroughNodeMaker."""
        self.name = "PassThroughNodeMaker"

    def get_node(self, control: Control, index: int):
        return PassThroughNode(control)
