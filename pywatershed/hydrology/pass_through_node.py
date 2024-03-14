from ..base.flow_graph import FlowNode, FlowNodeMaker
from ..constants import zero


class PassThroughNode(FlowNode):
    def __init__(self, control):
        self.control = control
        return

    def prepare_timestep(self):
        # self._simulation_time = simulation_time  # add argument?
        return

    def calculate_subtimestep(self, isubstep, inflow_upstream, inflow_lateral):
        self._seg_outflow = inflow_upstream
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
        return self._seg_outflow

    @property
    def storage_change(self):
        return zero


class PassThroughNodeMaker(FlowNodeMaker):
    def __init__(
        self,
    ) -> None:
        self.name = "PassThroughNodeMaker"

    def get_node(self, control, index):
        return PassThroughNode(control)
