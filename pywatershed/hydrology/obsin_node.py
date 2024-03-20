from ..base.flow_graph import FlowNode, FlowNodeMaker
from ..constants import zero


class ObsInNode(FlowNode):
    def __init__(
        self,
        control,
        node_obs_data,
    ):
        self.control = control
        self._node_obs_data = node_obs_data
        return

    def prepare_timestep(self):
        ymd = self.control.current_datetime.strftime("%Y-%m-%d")
        self._seg_outflow = self._node_obs_data[ymd]
        self._sink_source_sum = zero
        return

    def calculate_subtimestep(self, isubstep, inflow_upstream, inflow_lateral):
        inflow = inflow_upstream + inflow_lateral
        self._sink_source_sum += self._seg_outflow - inflow
        self._sink_source = self._sink_source_sum / (isubstep + 1)
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

    @property
    def sink_source(self):
        """Average sink and source through the last subtimestep
        Sink is negative, indicating that incoming flow is being discarded (if
        it were being stored, the storage change would be the opposite sign).
        Source is positive, indicating that incoming flow is being augmented
        (if it were being stored, the storage change would be the opposite
        sign).
        """
        return self._sink_source


class ObsInNodeMaker(FlowNodeMaker):
    def __init__(
        self,
        parameters,
        obs_data,
    ) -> None:
        self.name = "PassThroughNodeMaker"
        self._parameters = parameters
        self._obs_data = obs_data

    def get_node(self, control, index):
        node_poi_id = self._parameters.parameters["poi_gage_id"][index]
        return ObsInNode(control, self._obs_data[node_poi_id])
