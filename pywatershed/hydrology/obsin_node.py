import pandas as pd

from pywatershed.base.control import Control
from pywatershed.base.flow_graph import FlowNode, FlowNodeMaker
from pywatershed.base.parameters import Parameters
from pywatershed.constants import nan, zero


class ObsInNode(FlowNode):
    """A FlowNode that takes inflows but returns observed/specified flows.

    Args:
      control: a Control object.
      node_obs_data: A pandas Series object of observations at this location
        given by pyPRMS.Streamflow.
    """

    def __init__(
        self,
        control: Control,
        node_obs_data: pd.Series,
    ):
        self.control = control
        self._node_obs_data = node_obs_data
        return

    def prepare_timestep(self):
        ymd = self.control.current_datetime.strftime("%Y-%m-%d")
        self._seg_outflow = self._node_obs_data[ymd]
        self._sink_source_sum = zero
        return

    def calculate_subtimestep(
        self,
        isubstep: int,
        inflow_upstream: float,
        inflow_lateral: float,
    ):
        inflow = inflow_upstream + inflow_lateral
        if self._seg_outflow >= zero:
            self._sink_source_sum += self._seg_outflow - inflow
        else:
            self._seg_outflow = inflow
            self._sink_source_sum += zero

        # <
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
    def storage(self):
        return nan

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
    """A FlowNodeMaker for ObsInNode

    Args:
      parameters: A pywatershed Parameters object.
      obs_data: A pandas DataFrame of observations given by pyPRMS.Streamflow.
    """

    def __init__(
        self,
        parameters: Parameters,
        obs_data: pd.DataFrame,
    ) -> None:
        self.name = "PassThroughNodeMaker"
        self._parameters = parameters
        self._obs_data = obs_data

    def get_node(self, control: Control, index: int):
        node_poi_id = self._parameters.parameters["poi_gage_id"][index]
        return ObsInNode(control, self._obs_data[node_poi_id])
