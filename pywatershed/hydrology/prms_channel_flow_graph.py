from typing import Literal

import numba as nb
import numpy as np

from pywatershed.base.adapter import Adapter, adaptable
from pywatershed.base.conservative_process import ConservativeProcess
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
        calc_method: Literal["numba", "numpy"] = None,
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

        if calc_method is None or calc_method == "numba":
            self._calculate_subtimestep = _calculate_subtimestep_numba
        elif calc_method == "numpy":
            self._calculate_subtimestep = _calculate_subtimestep_numpy
        else:
            raise ValueError(f"Invalid choice of calc_method: {calc_method}")
        return

    def prepare_timestep(self):
        # self._simulation_time = simulation_time  # add argument?
        self._seg_inflow = zero
        self._seg_outflow = zero
        # self._seg_outflow0 = zero
        self._inflow_ts = zero
        self._seg_current_sum = zero

        return

    def calculate_subtimestep(self, ihr, inflow_upstream, inflow_lateral):
        (
            self._seg_inflow,
            self._inflow_ts,
            self._outflow_ts,
            self._seg_inflow0,
            self._seg_outflow,
        ) = self._calculate_subtimestep(
            ihr,
            inflow_upstream,
            inflow_lateral,
            self._seg_inflow0,
            self._seg_inflow,  # implied on RHS by +=
            self._inflow_ts,
            self._seg_outflow,
            self._outflow_ts,
            self._tsi,
            self._ts,
            self._c0,
            self._c1,
            self._c2,
        )

    def finalize_timestep(self):
        # get rid of the magic 24 with argument?
        self._seg_outflow = self._seg_outflow / 24.0
        self._seg_inflow = self._seg_inflow / 24.0
        self.seg_stor_change = self._seg_inflow - self._seg_outflow

        return

    def advance(self):
        self._seg_inflow0 = self._seg_inflow

        return

    @property
    def outflow(self):
        """The average outflow over the timestep in cubic feet per second."""
        return self._seg_outflow

    @property
    def outflow_substep(self):
        """The outflow over the sub-timestep in cubic feet per second."""
        return self._outflow_ts

    @property
    def storage_change(self):
        """The volumetric storage change in cubic feet."""
        return self.seg_stor_change

    @property
    def sink_source(self):
        return zero


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
        discretization: Parameters,
        parameters: Parameters,
        calc_method: Literal["numba", "numpy"] = None,
        verbose: bool = None,
    ) -> None:
        self.name = "PRMSChannelFlowNodeMaker"
        self._calc_method = calc_method
        self._set_data(discretization, parameters)
        self._init_data()

        return

    def get_node(self, control, index):
        # could pass initial conditons here but they arent currently used in
        # PRMS
        return PRMSChannelFlowNode(
            control=control,
            tsi=self._tsi[index],
            ts=self._ts[index],
            c0=self._c0[index],
            c1=self._c1[index],
            c2=self._c2[index],
            calc_method=self._calc_method,
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


# <
# module scope functions
def _calculate_subtimestep_numpy(
    ihr,
    inflow_upstream,
    inflow_lateral,
    _seg_inflow0,
    _seg_inflow,
    _inflow_ts,
    _seg_outflow,
    _outflow_ts,
    _tsi,
    _ts,
    _c0,
    _c1,
    _c2,
):
    # some way to enforce that ihr is actually an hour?
    seg_current_inflow = inflow_lateral + inflow_upstream
    _seg_inflow += seg_current_inflow
    _inflow_ts += seg_current_inflow

    remainder = (ihr + 1) % _tsi
    if remainder == 0:
        # segment routed on current hour
        _inflow_ts /= _ts

        if _tsi > 0:
            # Muskingum routing equation
            _outflow_ts = (
                _inflow_ts * _c0 + _seg_inflow0 * _c1 + _outflow_ts * _c2
            )
        else:
            _outflow_ts = _inflow_ts

        _seg_inflow0 = _inflow_ts
        _inflow_ts = 0.0

    _seg_outflow += _outflow_ts

    return (
        _seg_inflow,
        _inflow_ts,
        _outflow_ts,
        _seg_inflow0,
        _seg_outflow,
    )


numba_msg = "prms_channel_flow_graph jit compiling with numba"
print(numba_msg, flush=True)

_calculate_subtimestep_numba = nb.njit(
    nb.types.UniTuple(nb.float64, 5)(
        nb.int64,  # ihr
        nb.float64,  # inflow_upstream
        nb.float64,  # inflow_lateral
        nb.float64,  # _seg_inflow0
        nb.float64,  # _seg_inflow
        nb.float64,  # _inflow_ts
        nb.float64,  # _seg_outflow
        nb.float64,  # _outflow_ts
        nb.int64,  # _tsi
        nb.float64,  # _ts
        nb.float64,  # _c0
        nb.float64,  # _c1
        nb.float64,  # _c2
    ),
    fastmath=True,
    parallel=False,
)(_calculate_subtimestep_numpy)


class HruSegmentFlowAdapter(Adapter):
    """
    Adapt PRMS HRU volumetric outflows to PRMS Segment lateral inflows.
    The calculated lateral flows (in cubic feet per second) are availble from
    the current_value property.

    Args
    ----
    parameters: A PRMSChannel parameter object.
    sroff_vol: An Adapter of surface runoff volume
    ssres_flow_vol: An Adapter of subsurface reservoir outflow volume
    gwres_flow_vol: An Adapter of groundwater reservoir outflow volume
    """

    def __init__(
        self,
        parameters: Parameters,
        sroff_vol: Adapter,
        ssres_flow_vol: Adapter,
        gwres_flow_vol: Adapter,
    ):
        self._variable = "inflows"
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

        # sum of flows from hrus, these are flows and not volumes
        self._inflows = np.zeros(self._nhru) * nan

        # self.current_value (provided in the super) is the lateral flows
        self._current_value = np.zeros(self._nsegment) * nan

        self._hru_segment = self._parameters.parameters["hru_segment"] - 1

        return

    def advance(self):
        """Advance the Segement inflows in time.
        Here we move to flows from flow volumes.
        """
        for vv in self._inputs.values():
            vv.advance()

        # Move to a flow basis here instead of a volume basis.
        # s_per_time could vary with timestep so leave here.
        s_per_time = self._inputs["sroff_vol"].control.time_step_seconds
        self._inflows = (
            sum([vv.current for vv in self._inputs.values()]) / s_per_time
        )

        self._calculate_segment_lateral_inflows()

        return

    def _calculate_segment_lateral_inflows(self):
        """Map HRU inflows to lateral inflows on segments/nodes"""
        self._current_value[:] = zero
        for ihru in range(self._nhru):
            iseg = self._hru_segment[ihru]
            if iseg < 0:
                # This is bad, selective handling of fluxes is not cool,
                # mass is being discarded in a way that has to be
                # coordinated with other parts of the code.
                # This code should be removed evenutally.
                continue

            self._current_value[iseg] += self._inflows[ihru]

        return


class HruSegmentFlowExchange(ConservativeProcess):
    """PRMS class to map HRU outflows to segment inflows.

    This code is contained in the channel part of the PRMS code.

    Args:

    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        sroff_vol: adaptable,
        ssres_flow_vol: adaptable,
        gwres_flow_vol: adaptable,
        budget_type: Literal[None, "warn", "error"] = None,
        verbose: bool = None,
    ) -> None:
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "HruSegmentFlowExchange"

        self._set_inputs(locals())
        self._set_options(locals())

        self._hru_segment = self.hru_segment - 1

        self._set_budget(basis="global")

        return

    @staticmethod
    def get_dimensions() -> tuple:
        return ("nhru", "nnodes")

    @staticmethod
    def get_parameters() -> tuple:
        return ("hru_segment", "node_coord")

    @staticmethod
    def get_inputs() -> tuple:
        return (
            "sroff_vol",
            "ssres_flow_vol",
            "gwres_flow_vol",
        )

    @staticmethod
    def get_init_values() -> dict:
        return {
            "inflows": nan,
            "inflows_vol": nan,
            "sinks": nan,
            "sinks_vol": nan,
        }

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "sroff_vol",
                "ssres_flow_vol",
                "gwres_flow_vol",
            ],
            "outputs": ["inflows_vol", "sinks_vol"],
            "storage_changes": [],
        }

    def _set_initial_conditions(self) -> None:
        pass

    def _advance_variables(self) -> None:
        # could vary with timesstep some day
        return

    def _calculate(self, simulation_time: float) -> None:
        self._simulation_time = simulation_time

        s_per_time = self.control.time_step_seconds
        self._inputs_sum = (
            sum([vv.current for vv in self._input_variables_dict.values()])
            / s_per_time
        )

        self.inflows[:] = zero
        self.sinks[:] = zero  # this is an HRU variable
        for ihru in range(self.nhru):
            iseg = self._hru_segment[ihru]
            if iseg < 0:
                self.sinks[ihru] += self._inputs_sum[ihru]
            else:
                self.inflows[iseg] += self._inputs_sum[ihru]

        self.inflows_vol[:] = self.inflows * s_per_time
        self.sinks_vol[:] = self.sinks * s_per_time

        return
