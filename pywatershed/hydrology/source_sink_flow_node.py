import numpy as np
import pandas as pd

from pywatershed.base.control import Control
from pywatershed.base.flow_graph import FlowNode, FlowNodeMaker
from pywatershed.base.parameters import Parameters
from pywatershed.constants import nan, zero


class SourceSinkFlowNode(FlowNode):
    """A FlowNode that adds or removes flow above some minimum flow parameter.

    See :class:`FlowGraph` for examples and discussion.
    """

    def __init__(
        self,
        control: Control,
        flow_min: np.float64,
        source_sink_data: pd.Series,
        # JLM TODO: budget
    ):
        """Initialize an SourceSinkFlowNode.

        Args:
          control: a Control object.
          flow_min: A floating point value for the minium flow.
          source_sink_data: A pandas Series object of sources/sinks at this
            location. See SourceSinkFlowNodeMaker for a description of the
            pd.DataFrame passed to supply this data.
        """
        self.control = control
        self._flow_min = flow_min
        self._source_sink_data = source_sink_data
        return

    # @staticmethod
    # def get_mass_budget_terms():
    #     """Get a dictionary of variable names for mass budget terms."""
    #     return {
    #         "inputs": [
    #             "_lake_inflow",
    #         ],
    #         "outputs": [
    #             "_lake_release",
    #             "_lake_spill",
    #         ],
    #         "storage_changes": [
    #             "_lake_storage_change_flow_units",
    #         ],
    #     }

    def prepare_timestep(self):
        ymd = self.control.current_datetime.strftime("%Y-%m-%d")
        self._source_sink_requested = self._source_sink_data[ymd]
        self._sink_source_sum = zero
        return

    def calculate_subtimestep(
        self,
        isubstep: int,
        inflow_upstream: float,
        inflow_lateral: float,
    ):
        inflow = inflow_upstream + inflow_lateral
        source_sink = self._source_sink_requested
        min_flow = self._flow_min

        if source_sink >= zero:
            # a source is always applied
            outflow = inflow + source_sink

        elif (source_sink < zero) and (inflow < min_flow):
            # sink is not applied when inflow < min_flow
            outflow = inflow
            source_sink = zero
            source_sink = zero

        elif (source_sink < zero) and (inflow >= min_flow):
            if inflow + source_sink < min_flow:
                source_sink = min_flow - inflow
                outflow = min_flow
            else:
                outflow = inflow + source_sink

        # <
        self._seg_outflow = outflow
        self._sink_source_sum += source_sink
        self._sink_source = self._sink_source_sum / (isubstep + 1)

        return

    def finalize_timestep(self):
        return

    def advance(self):
        # if self.budget is not None:
        #     self.budget.advance()
        #     self.budget.calculate()

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


class SourceSinkFlowNodeMaker(FlowNodeMaker):
    """A FlowNodeMaker for SourceSinkFlowNode.

    See :class:`FlowGraph` for related examples and discussion.
    """

    def __init__(
        self,
        parameters: Parameters,
        source_sink_df: pd.DataFrame,
    ) -> None:
        """Initialize a ObsInFlowNodeMaker.

        Args:
          parameters: A pywatershed Parameters object.
          source_sink_df: A pandas DataFrame of observations with a time
            index which can be selected by '%Y-%m-%d' strftime of a datetime64.
            The column names are not used and may be anything, for example the
            nhm_seg of the upstream segment. However, the columns order MUST
            be collated with the input data vectors. The sign convention is:
            sources are positive and sinks are negative. That is, the sign is
            from the perspective of the node.

        """
        self.name = "SourceSinkNodeMaker"
        self._parameters = parameters
        self._source_sink_df = source_sink_df

    def get_node(self, control: Control, index: int):
        flow_min = self._parameters.parameters["flow_min"][index]
        col_name = self._source_sink_df.columns[index]
        data_ts = self._source_sink_df[col_name]
        return SourceSinkFlowNode(control, flow_min, data_ts)
