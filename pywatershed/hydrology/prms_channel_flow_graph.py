from typing import Literal

import numpy as np

from pywatershed.base.adapter import Adapter
from pywatershed.base.control import Control
from pywatershed.base.flow_graph import FlowNode, FlowNodeMaker
from pywatershed.constants import SegmentType, nan, zero
from pywatershed.parameters import Parameters


class PRMSChannelFlowNode(FlowNode):
    def __init__(
        self,
        control: Control,
        tsi: np.int64,
        ts: np.float64,
        c0: np.float64,
        c1: np.float64,
        c2: np.float64,
    ):
        self.control = control

        self._tsi = tsi
        self._ts = ts
        self._c0 = c0
        self._c1 = c1
        self._c2 = c2

        self._outflow_ts = zero
        self._seg_inflow0 = zero
        self._seg_inflow = zero
        self._inflow_ts = zero
        self._seg_current_sum = zero

        return

    def prepare_timestep(self):
        # self._simulation_time = simulation_time  # add argument?

        self._s_per_time = self.control.time_step_seconds

        self._seg_inflow = zero
        self._seg_outflow = zero
        # self._seg_outflow0 = zero
        self._inflow_ts = zero
        self._seg_current_sum = zero

        return

    def calculate_subtimestep(self, ihr, inflow_upstream, inflow_lateral):
        # some way to enforce that ihr is actually an hour?
        seg_current_inflow = inflow_lateral + inflow_upstream
        self._seg_inflow += seg_current_inflow
        self._inflow_ts += seg_current_inflow

        remainder = (ihr + 1) % self._tsi
        if remainder == 0:
            # segment routed on current hour
            self._inflow_ts /= self._ts

            if self._tsi > 0:
                # Muskingum routing equation
                self._outflow_ts = (
                    self._inflow_ts * self._c0
                    + self._seg_inflow0 * self._c1
                    + self._outflow_ts * self._c2
                )
            else:
                self._outflow_ts = self._inflow_ts

            self._seg_inflow0 = self._inflow_ts
            self._inflow_ts = 0.0

        self._seg_outflow += self._outflow_ts

        return

    def finalize_timestep(self):
        # get rid of the magic 24 with argument?
        self._seg_outflow = self._seg_outflow / 24.0
        self._seg_inflow = self._seg_inflow / 24.0
        # self.seg_upstream_inflow = self._seg_current_sum / 24.0

        self.seg_stor_change = (
            self._seg_inflow - self._seg_outflow
        ) * self._s_per_time

        return

    def advance(self):
        self._seg_inflow0 = self._seg_inflow

        return

    @property
    def outflow(self):
        return self._seg_outflow

    @property
    def outflow_substep(self):
        return self._outflow_ts

    @property
    def storage_change(self):
        return self.seg_stor_change


class PRMSChannelFlowNodeMaker(FlowNodeMaker):
    """PRMS channel flow (muskingum_mann) using a FlowGraph.

    See PRMSChannel for details on the method.
    This class bases on FlowGraphNodeMaker which bases on
    ConservativeProcess.

    Args:
        discretization: a discretization of class Parameters
        parameters: a parameter object of class Parameters
        calc_method: one of ["fortran", "numba", "numpy"]. None defaults to
            "numba".
        verbose: Print extra information or not?
    """

    def __init__(
        self,
        NodeClass: FlowNode,
        discretization: Parameters,
        parameters: Parameters,
        calc_method: Literal["fortran", "numba", "numpy"] = None,
        verbose: bool = None,
    ) -> None:
        self.name = "PRMSChannelFlowNodeMaker"
        self._NodeClass = NodeClass
        self._set_data(discretization, parameters)

        self._init_data()

        return

    def get_node(self, control, index):
        # could pass initial conditons here but they arent currently used in
        # PRMS
        return self._NodeClass(
            control=control,
            tsi=self._tsi[index],
            ts=self._ts[index],
            c0=self._c0[index],
            c1=self._c1[index],
            c2=self._c2[index],
        )

    def _set_data(self, discretization, parameters):
        self._parameters = parameters
        self._discretization = discretization
        for param in self.get_parameters():
            if param in self._parameters.parameters.keys():
                self[param] = self._parameters.parameters[param]
            else:
                self[param] = self._discretization.parameters[param]

        self.nsegment = len(self.seg_length)
        return

    @staticmethod
    def get_dimensions() -> tuple:
        return "nsegment"

    @staticmethod
    def get_parameters() -> tuple:
        return (
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

    def _init_data(self) -> None:
        """Initialize internal variables from raw channel data"""

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
        # JLM: This is a bad idea and should throw an error rather than edit
        # inputs in place during run
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
        self._tsi = np.ones(self.nsegment, dtype="int64")

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
        self._seg_inflow0 = np.zeros(self.nsegment, dtype=float) * nan
        self._inflow_ts = np.zeros(self.nsegment, dtype=float)
        self._outflow_ts = np.zeros(self.nsegment, dtype=float)
        self._seg_current_sum = np.zeros(self.nsegment, dtype=float)

        return


class AdapterExchangeHruSegment(Adapter):
    def __init__(
        self,
        variable: str,
        parameters: Parameters,
        sroff_vol: Adapter,
        ssres_flow_vol: Adapter,
        gwres_flow_vol: Adapter,
    ):
        self._variable = variable
        self._parameters = parameters

        self._inputs = {}
        for key in [
            "sroff_vol",
            "ssres_flow_vol",
            "gwres_flow_vol",
        ]:
            self._inputs[key] = locals()[key]

        self._nhru = self._parameters.dims["nhru"]
        self._nsegment = self._parameters.dims["nsegment"]
        # flows from hrus
        self._inflows = np.zeros(self._nhru) * nan
        # flows to segments
        self._current_value = np.zeros(self._nsegment) * nan

        self._hru_segment = self._parameters.parameters["hru_segment"] - 1

        return

    def advance(self):
        for vv in self._inputs.values():
            vv.advance()

        self._inflows = sum([vv.current for vv in self._inputs.values()])

        self._calculate_segment_inflows()

        return

    def _calculate_segment_inflows(self):
        """Map HRU flows on to segments"""
        # This could vary with timestep so leave here
        self._s_per_time = self._inputs["sroff_vol"].control.time_step_seconds

        self._current_value[:] = zero
        for ihru in range(self._nhru):
            iseg = self._hru_segment[ihru]
            if iseg < 0:
                # This is bad, selective handling of fluxes is not cool,
                # mass is being discarded in a way that has to be
                # coordinated with other parts of the code.
                # This code shuold be removed evenutally.
                continue

            # cubicfeet to cfs
            segment_inflow = self._inflows[ihru] / (self._s_per_time)
            self._current_value[iseg] += segment_inflow

        return

    # def _calculate_lateral_inflow_inds(self):
    #     """A hash from segments to HRU"""
    #     self._seg_hrus = {}
    #     for ihru in range(self.nhru):
    #         iseg = self._hru_segment[ihru]  # to be removed via exchange
    #         if iseg in self._seg_hrus.keys():
    #             self._seg_hrus[iseg] += [ihru]
    #         else:
    #             self._seg_hrus[iseg] = [ihru]
