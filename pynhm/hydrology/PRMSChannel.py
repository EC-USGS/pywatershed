from typing import Tuple, Union

import networkx as nx
import numpy as np
from numba import jit

from pynhm.base.storageUnit import StorageUnit

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import SegmentType, zero

FT2_PER_ACRE = 43560.0
INCHES_PER_FOOT = 12.0
SECS_PER_HOUR = 3600.0


class PRMSChannel(StorageUnit):
    """PRMS channel flow (muskingum_mann)

    The muskingum module was originally developed for the Precipitation Runoff
    Modeling System (PRMS) by Mastin and Vaccaro (2002) and developed further
    by Markstrom and others (2008). This module has been modified from past
    versions to make it more stable for stream network routing in watersheds
    with stream segments with varying travel times. Although this module runs
    on the same daily time step as the rest of the modules in PRMS, it has an
    internal structure which allows for a different computational time step for
    each segment in the stream network, ensuring that the simulation produces
    stable values. Flow values computed at these finer time steps are
    aggregated by the Muskingum module to provide consistent daily time step
    values, regardless of the segment length.

    Delta t, which is the travel time (in hours), is rounded down to an even
    divisor of 24 hours (for example 24, 12, 6, 4, 3, 2, and 1). PRMS is
    restricted to daily time steps, so Delta t segment can never be more than
    one day in length. This means that the travel time of any segment in the
    stream network (K_coef) must be less than one day. An implication of this
    is that the routed streamflow in each segment is computed using different
    solution time steps. Consequently, streamflow must be aggregated when
    flowing from segments with shorter Delta t segment to segments with longer
    Delta t. Likewise, streamflow must be disaggregated when flowing from
    segments with longer Delta to shorter Delta t. In either case, flow
    from upstream segments is averaged and summed to the appropriate value of
    Delta t.

    The muskingum_mann method is a modified version of the original muskingum
    function in PRMS that was introduced in PRMS version 5.2.1 (1/20/2021).
    The muskingum_mann method provides a method to calculate K_coef values
    using mann_n, seg_length, seg_depth (bank full), and seg_slope. The
    velocity at bank full segment depth is calculated using Manning's equation

        velocity = ((1/n) sqrt(seg_slope) seg_depth**(2/3)

    K_coef ,in hours, is then calculated using

        K_coef = seg_length / (velocity * 60 * 60)

    K_coef values computed greater than 24.0 are set to 24.0, values computed
    less than 0.01 are set to 0.01, and the value for lake HRUs is set to 24.0.

    Args:
        control: control object
        sroff: surface runoff adapter object
        ssres_flow: subsurface (gravity) reservoir lateral flow adapter object
        gwres_flow: groundwater reservoir baseflow adapter object
        verbose: verbose output boolean (default is False)


    """

    def __init__(
        self,
        control: Control,
        sroff: adaptable,
        ssres_flow: adaptable,
        gwres_flow: adaptable,
        verbose: bool = False,
        budget_type: str = None,
    ) -> "PRMSChannel":

        super().__init__(
            control=control,
            verbose=verbose,
        )
        self.name = "PRMSChannel"

        self.set_inputs(locals())
        self.set_budget(budget_type)

        # process channel data
        self._initialize_channel_data()

        return

    def _initialize_channel_data(self) -> None:
        """Initialize internal variables from raw channel data"""
        # convert prms data to zero-based
        self.hru_segment -= 1
        self.tosegment -= 1

        # calculate connectivity
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

        # use networkx to calculate the Directed Acyclic Graph
        if self.nsegment > 1:
            graph = nx.DiGraph()
            graph.add_edges_from(connectivity)
            segment_order = list(nx.topological_sort(graph))
        else:
            segment_order = [0]
        self._segment_order = np.array(segment_order, dtype=int)

        # calculate the Muskingum parameters
        velocity = (
            (
                (1.0 / self.mann_n)
                * np.sqrt(self.seg_slope)
                * self.seg_depth ** (2.0 / 3.0)
            )
            * 60.0
            * 60.0
        )
        self.seg_slope = np.where(
            self.seg_slope < 1e-7, 0.0001, self.seg_slope
        )  # not in prms6

        # initialize Kcoef to 24.0 for segments with zero velocities
        # this is different from PRMS, which relied on divide by zero resulting
        # in a value of infinity that when evaluated relative to a maximum
        # desired Kcoef value of 24 would be reset to 24. This approach is
        # equivalent and avoids the occurence of a divide by zero.
        Kcoef = np.full(self.nsegment, 24.0, dtype=float)

        # only calculate Kcoef for cells with velocities greater than zero
        idx = velocity > 0.0
        Kcoef[idx] = self.seg_length[idx] / velocity[idx]
        Kcoef = np.where(
            self.segment_type == SegmentType.LAKE.value, 24.0, Kcoef
        )
        Kcoef = np.where(Kcoef < 0.01, 0.01, Kcoef)
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
        d = np.where(np.abs(d) < 1e-6, 0.0001, d)
        self._c0 = (-(self._Kcoef * self.x_coef) + (0.5 * self._ts)) / d
        self._c1 = ((self._Kcoef * self.x_coef) + (0.5 * self._ts)) / d
        self._c2 = (
            self._Kcoef - (self._Kcoef * self.x_coef) - (0.5 * self._ts)
        ) / d

        # Short travel time
        idx = self._c2 < 0.0
        self._c1[idx] += self._c2[idx]
        self._c2[idx] = 0.0

        # Long travel time
        idx = self._c0 < 0.0
        self._c1[idx] += self._c0[idx]
        self._c0[idx] = 0.0

        # local flow variables
        self._seg_inflow = np.zeros(self.nsegment, dtype=float)
        self._seg_inflow0 = np.zeros(self.nsegment, dtype=float)
        self._seg_outflow0 = np.zeros(self.nsegment, dtype=float)
        self._inflow_ts = np.zeros(self.nsegment, dtype=float)
        self._outflow_ts = np.zeros(self.nsegment, dtype=float)
        self._seg_current_sum = np.zeros(self.nsegment, dtype=float)

        # initialize internal sef_inflow variable
        for iseg in range(self.nsegment):
            jseg = self.tosegment[iseg]
            if jseg < 0:
                continue
            self._seg_inflow[jseg] = self.seg_outflow[iseg]

        return

    def set_initial_conditions(self) -> None:
        # initialize channel segment storage
        self.seg_outflow = self.segment_flow_init
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
            "obsin_segment",
            "obsout_segment",
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
        self._seg_inflow0[:] = self._seg_inflow.copy()
        self._seg_outflow0[:] = self.seg_outflow.copy()
        return

    def calculate(self, simulation_time: float) -> None:
        """Calculate channel segment terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

        hru_area = self.hru_area
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

        self.seg_upstream_inflow, self.seg_outflow = muskingum_routing(
            self._segment_order,
            self.tosegment,
            self.seg_lateral_inflow,
            self.seg_upstream_inflow,
            self._seg_inflow,
            self._seg_inflow0,
            self.seg_outflow,
            self._seg_outflow0,
            self._inflow_ts,
            self._outflow_ts,
            self._seg_current_sum,
            self._tsi,
            self._ts,
            self._c0,
            self._c1,
            self._c2,
        )

        return


@jit(nopython=True)
def muskingum_routing(
    segment_order: np.ndarray,
    to_segment: np.ndarray,
    seg_lateral_inflow: np.ndarray,
    seg_upstream_inflow: np.ndarray,
    seg_inflow: np.ndarray,
    seg_inflow0: np.ndarray,
    seg_outflow: np.ndarray,
    seg_outflow0: np.ndarray,
    inflow_ts: np.ndarray,
    outflow_ts: np.ndarray,
    seg_current_sum: np.ndarray,
    tsi: np.ndarray,
    ts: np.ndarray,
    c0: np.ndarray,
    c1: np.ndarray,
    c2: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Muskingum routing function that calculates the upstream inflow and
    outflow for each segment

    Args:
        segment_order: segment routing order
        to_segment: downstream segment for each segment
        seg_lateral_inflow: segment lateral inflow
        seg_upstream_inflow: initial segment upstream inflow (set to 0)
        seg_inflow: segment inflow variable (internal calculations)
        seg_inflow0: previous segment inflow variable (internal calculations)
        seg_outflow: initial segment outflow variable (set to 0)
        seg_outflow0: previous outflow variable (internal calculations)
        inflow_ts: inflow timeseries variable (internal calculations)
        outflow_ts: outflow timeseries variable (internal calculations)
        seg_current_sum: summation variable (internal calculations)
        tsi: integer flood wave travel time
        ts: float version of integer flood wave travel time
        c0: Muskingum c0 variable
        c1: Muskingum c1 variable
        c2: Muskingum c2 variable

    Returns:
        seg_upstream_inflow: inflow for each segment for the current day
        seg_outflow: outflow for each segment for the current day

    """
    # initialize variables for the day
    # seg_inflow0[:] = seg_inflow.copy()
    # seg_outflow0[:] = seg_outflow.copy()
    seg_inflow[:] = 0.0
    seg_outflow[:] = 0.0
    inflow_ts[:] = 0.0
    seg_current_sum[:] = 0.0

    for ihr in range(24):
        seg_upstream_inflow[:] = 0.0

        for jseg in segment_order:
            # current inflow to the segment is the time-weighted average
            # of the outflow of the upstream segments and the lateral HRU
            # inflow plus any gains
            seg_current_inflow = seg_lateral_inflow[jseg]

            # todo: evaluate if obsin_segment needs to be implemented -
            #  would be needed needed if headwater basins are not included
            #  in a simulation

            seg_current_inflow += seg_upstream_inflow[jseg]
            seg_inflow[jseg] += seg_current_inflow
            inflow_ts[jseg] += seg_current_inflow
            seg_current_sum[jseg] += seg_upstream_inflow[jseg]

            remainder = (ihr + 1) % tsi[jseg]
            if remainder == 0:
                # segment routed on current hour
                inflow_ts[jseg] /= ts[jseg]

                if tsi[jseg] > 0:

                    # todo: evaluated if denormal results should be dealt with

                    # Muskingum routing equation
                    outflow_ts[jseg] = (
                        inflow_ts[jseg] * c0[jseg]
                        + seg_inflow0[jseg] * c1[jseg]
                        + outflow_ts[jseg] * c2[jseg]
                    )
                else:
                    # travel time is 1 hour or less so outflow is set
                    # equal to the inflow - outflow_ts is the value for
                    # the previous hour
                    outflow_ts[jseg] = inflow_ts[jseg]

                # previous inflow is equal to inflow_ts from the previous
                # routed time step
                seg_inflow0[jseg] = inflow_ts[jseg]

                # upstream inflow is used, reset it to zero so a new
                # average can be calculated next routing time step
                inflow_ts[jseg] = 0.0

            # todo: evaluate if obsout_segment needs to be implemented -
            #  would be needed needed fixing ourflow to observed data is
            #  required in a simulation

            # todo: water use

            # segment outflow (the mean daily flow rate for each segment)
            # will be the average of hourly values
            seg_outflow[jseg] += outflow_ts[jseg]

            # previous segment outflow is equal to the inflow_ts on the
            # previous routed timestep
            seg_outflow0[jseg] = outflow_ts[jseg]

            # add current time step flow rate to the upstream flow rate
            # for the segment this segment is connected to
            to_seg = to_segment[jseg]
            if to_seg >= 0:
                seg_upstream_inflow[to_seg] += outflow_ts[jseg]

    seg_outflow /= 24.0
    seg_inflow /= 24.0
    seg_upstream_inflow = seg_current_sum.copy() / 24.0

    return seg_upstream_inflow, seg_outflow
