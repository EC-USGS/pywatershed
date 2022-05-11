from typing import Union

import networkx as nx
import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control
from ..constants import zero

adaptable = Union[str, np.ndarray, Adapter]

FT2_PER_ACRE = 43560.0
INCHES_PER_FOOT = 12.0
SECS_PER_HOUR = 3600.0


class PRMSChannel(StorageUnit):
    """PRMS channel flow

    Args:
        params: parameter object

    """

    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        sroff: adaptable,
        ssres_flow: adaptable,
        gwres_flow: adaptable,
        verbose: bool = False,
    ) -> "PRMSChannel":

        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
        )

        self.name = "PRMSChannel"

        self._input_variables_dict = {}
        self._input_variables_dict["sroff"] = adapter_factory(
            sroff,
            "sroff",
        )
        self._input_variables_dict["ssres_flow"] = adapter_factory(
            ssres_flow,
            "ssres_flow",
        )
        self._input_variables_dict["gwres_flow"] = adapter_factory(
            gwres_flow, "gwres_flow"
        )

        # process channel data
        self._process_channel_data()

        return

    def _process_channel_data(self):
        # convert prms data to zero-based
        self.hru_segment -= 1

        # calculate connectivity
        self.tosegment -= 1

        connectivity = []
        for iseg in range(self.nsegment):
            tosegment = self.tosegment[iseg]
            if tosegment < 0:
                continue
            connectivity.append(
                (
                    iseg,
                    tosegment,
                )
            )

        graph = nx.DiGraph()
        graph.add_edges_from(connectivity)
        segment_order = list(nx.topological_sort(graph))
        if len(segment_order) < 1:
            segment_order = [i for i in range(self.nsegment)]
        self._segment_order = segment_order

        # calculate the Muskingum parameters
        velocity = (
            (1.0 / self.mann_n)
            * np.sqrt(self.seg_slope)
            * self.seg_depth ** (2.0 / 3.0)
        )
        Kcoef = np.full(self.nsegment, 0.01, dtype=float)
        Kcoef = np.where(
            velocity > 0., self.seg_length / (velocity * 60.0 * 60.0), Kcoef
        )
        self._Kcoef = np.where(Kcoef > 24.0, 24.0, Kcoef)

        self._ts = np.ones(self.nsegment, dtype=float)
        self._tsi = np.ones(self.nsegment, dtype=int)

        # todo: vectorize this
        for iseg in range(self.nsegment):
            k = self._Kcoef[iseg]
            if k < 1.0:
                self._tsi[iseg] = -1
            elif k < 2.0:
                self._ts[iseg] = 1.0
                self._tsi[iseg] = 1
            elif k < 3.0:
                self._ts[iseg] = 2.0
                self._tsi[iseg] = 2
            elif k < 4.0:
                self._ts[iseg] = 3.0
                self._tsi[iseg] = 3
            elif k < 6.0:
                self._ts[iseg] = 4.0
                self._tsi[iseg] = 4
            elif k < 8.0:
                self._ts[iseg] = 6.0
                self._tsi[iseg] = 6
            elif k < 12.0:
                self._ts[iseg] = 8.0
                self._tsi[iseg] = 8
            elif k < 24.0:
                self._ts[iseg] = 12.0
                self._tsi[iseg] = 12
            else:
                self._ts[iseg] = 24.0
                self._tsi[iseg] = 24

        d = self._Kcoef - (self._Kcoef * self.x_coef) + (0.5 * self._ts)
        self._c0 = (-(self._Kcoef * self.x_coef) + (0.5 * self._ts)) / d
        self._c1 = ((self._Kcoef * self.x_coef) + (0.5 * self._ts)) / d
        self._c2 = (
            self._Kcoef - (self._Kcoef * self.x_coef) - (0.5 * self._ts)
        ) / d

        # Short travel time
        idx = self._c2 < 0.0
        self._c1[idx] = self._c1[idx] + self._c2[idx]
        self._c2[idx] = 0.0

        # Long travel time
        idx = self._c0 < 0.0
        self._c1[idx] = self._c1[idx] + self._c0[idx]
        self._c0[idx] = 0.0

        # local flow variables
        self._seg_inflow = np.zeros(self.nsegment, dtype=float)
        self._seg_inflow0 = np.zeros(self.nsegment, dtype=float)
        self._seg_outflow0 = np.zeros(self.nsegment, dtype=float)
        self._inflow_ts = np.zeros(self.nsegment, dtype=float)
        self._outflow_ts = np.zeros(self.nsegment, dtype=float)
        self._seg_current_sum = np.zeros(self.nsegment, dtype=float)

        # initialize variables
        self.seg_outflow = self.segment_flow_init
        for iseg in range(self.nsegment):
            jseg = self.tosegment[iseg]
            if jseg < 0:
                continue
            self._seg_inflow[jseg] = self.seg_outflow[iseg]

        return

    def set_initial_conditions(self):
        # initialize channel segment storage
        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get channel segment parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nhru",
            "nsegment",
            "hru_area",
            "hru_segment",
            "mann_n",
            "seg_depth",
            "seg_length",
            "seg_slope",
            "segment_type",
            "tosegment",
            "tosegment_nhm",
            "x_coef",
            "segment_flow_init",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get channel segment input variables

        Returns:
            variables: input variables

        """
        return (
            "sroff",
            "ssres_flow",
            "gwres_flow",
        )

    @staticmethod
    def get_variables() -> tuple:
        """Get channel segment output variables

        Returns:
            variables: output variables

        """
        return (
            "seg_lateral_inflow",
            "seg_upstream_inflow",
            "seg_outflow",
        )

    @staticmethod
    def get_init_values() -> dict:
        """Get channel segment initial values

        Returns:
            dict: initial values for named variables
        """
        # No GW res values need initialized prior to calculation.
        return {
            "seg_lateral_inflow": zero,
            "seg_upstream_inflow": zero,
            "seg_outflow": zero,
        }

    def _advance_variables(self) -> None:
        """Advance the channel segment variables
        Returns:
            None
        """
        return

    def calculate(self, simulation_time):
        """Calculate channel segment terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

        hru_area = self.hru_area
        # FT2_PER_ACRE / INCHES_PER_FOOT / Timestep_seconds
        time_step_seconds = self.control.time_step / np.timedelta64(1, "s")
        in_to_cfs = hru_area * FT2_PER_ACRE
        in_to_cfs /= INCHES_PER_FOOT
        in_to_cfs /= time_step_seconds

        # calculate lateral flow term
        self.seg_lateral_inflow[:] = 0.0
        for ihru in range(self.nhru):
            iseg = self.hru_segment[ihru]
            if iseg < 0:
                continue
            lateral_inflow = (
                self.sroff[ihru]
                + self.ssres_flow[ihru]
                + self.gwres_flow[ihru]
            ) * in_to_cfs[ihru]
            self.seg_lateral_inflow[iseg] += lateral_inflow

        # initialize variables for the day
        self._seg_inflow0 = self._seg_inflow.copy()
        self._seg_outflow0 = self.seg_outflow.copy()
        self._seg_inflow[:] = 0.0
        self.seg_outflow[:] = 0.0
        self._inflow_ts[:] = 0.0
        self._seg_current_sum[:] = 0.0

        for ihr in range(24):
            self.seg_upstream_inflow[:] = 0.0

            for iseg in range(self.nsegment):
                jseg = self._segment_order[iseg]

                # current inflow to the segment is the time-weighted average
                # of the outflow of the upstream segments and the lateral HRU
                # inflow plus any gains
                seg_current_inflow = self.seg_lateral_inflow[jseg]

                # todo: evaluate if obsin_segment needs to be implemented -
                #  would be needed needed if headwater basins are not included
                #  in a simulation

                seg_current_inflow += self.seg_upstream_inflow[jseg]
                self._seg_inflow[jseg] += seg_current_inflow
                self._inflow_ts[jseg] += seg_current_inflow
                self._seg_current_sum[jseg] += self.seg_upstream_inflow[jseg]

                remainder = (ihr + 1) % self._tsi[jseg]
                if remainder == 0:
                    # segment routed on current hour
                    self._inflow_ts[jseg] /= self._ts[jseg]

                    if self._tsi[jseg] > 0:

                        # todo: evaluated if denormal results should be dealt with

                        # Muskingum routing equation
                        self._outflow_ts[jseg] = (
                            self._inflow_ts[jseg] * self._c0[jseg]
                            + self._seg_inflow0[jseg] * self._c1[jseg]
                            + self._outflow_ts[jseg] * self._c2[jseg]
                        )
                    else:
                        # travel time is 1 hour or less so outflow is set
                        # equal to the inflow - outflow_ts is the value for
                        # the previous hour
                        self._outflow_ts[jseg] = self._inflow_ts[jseg]

                    # previous inflow is equal to inflow_ts from the previous
                    # routed time step
                    self._seg_inflow0[jseg] = self._inflow_ts[jseg]

                    # upstream inflow is used, reset it to zero so a new
                    # average can be calculated next routing time step
                    self._inflow_ts[jseg] = 0.0

                # todo: evaluate if obsout_segment needs to be implemented -
                #  would be needed needed fixing ourflow to observed data is
                #  required in a simulation

                # todo: water use

                # segment outflow (the mean daily flow rate for each segment)
                # will be the average of hourly values
                self.seg_outflow[jseg] += self._outflow_ts[jseg]

                # previous segment outflow is equal to the inflow_ts on the
                # previous routed timestep
                self._seg_outflow0[jseg] = self._outflow_ts[jseg]

                # add current time step flow rate to the upstream flow rate
                # for the segemnt this segment is connected to
                toseg = self.tosegment[jseg]
                if toseg >= 0:
                    self.seg_upstream_inflow[toseg] += self._outflow_ts[jseg]

        self.seg_outflow /= 24.0
        self._seg_inflow /= 24.0
        self.seg_upstream_inflow = self._seg_current_sum.copy() / 24.0

        return
