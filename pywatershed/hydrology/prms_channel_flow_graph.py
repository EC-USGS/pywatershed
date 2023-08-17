from typing import Literal, Tuple
from warnings import warn

import networkx as nx
import numpy as np

from ..base.adapter import adaptable
from ..base.control import Control

from ..base.flow_graph import FlowGraph, FlowNode
from ..constants import SegmentType, nan, zero
from ..parameters import Parameters


class PRMSChannelFlowNode(FlowNode):
    def __init__(self, control, tsi, ts, c0, c1, c2, graph_outflow):
        self.control = control

        self._tsi = tsi
        self._ts = ts
        self._c0 = c0
        self._c1 = c1
        self._c2 = c2
        self._graph_outflow = graph_outflow

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
        self._seg_outflow0 = zero
        self._inflow_ts = zero
        self._seg_current_sum = zero

        return

    def calculate_subtimestep(
        self, ihr, inflow_upstream, inflow_lateral, seg_ind=None
    ):
        # self.seg_upstream_inflow = zero  ## correct ? or in _calculate?
        # self._seg_current_sum += inflow_upstream
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
        self._seg_outflow = self._seg_outflow / 24.0
        self._seg_inflow = self._seg_inflow / 24.0
        # self.seg_upstream_inflow = self._seg_current_sum / 24.0

        self.seg_stor_change = (
            self._seg_inflow - self._seg_outflow
        ) * self._s_per_time

        self.channel_outflow_vol = zero
        if self._graph_outflow:
            self.channel_outflow_vol = self._seg_outflow * self._s_per_time

        return

    def advance(self):
        self._seg_inflow0 = self._seg_inflow.copy()

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


class PRMSChannelFlowGraph(FlowGraph):
    """PRMS channel flow (muskingum_mann) using a FlowGraph.

    See PRMSChannel for details on the method.
    This class bases on FlowGraph while that one bases on ConservativeProcess.

    Args:
        control: a Control object
        discretization: a discretization of class Parameters
        parameters: a parameter object of class Parameters
        sroff_vol: Surface runoff to the stream network for each HRU
        ssres_flow_vol: Interflow volume from gravity and preferential-flow
            reservoirs to the stream network for each HRU
        gwres_flow_vol: Groundwater discharge volume from each GWR to the
            stream network
        budget_type: one of [None, "warn", "error"]
        calc_method: one of ["fortran", "numba", "numpy"]. None defaults to
            "numba".
        verbose: Print extra information or not?
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
        calc_method: Literal["fortran", "numba", "numpy"] = None,
        verbose: bool = None,
    ) -> None:
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "PRMSChannel"

        self._set_inputs(locals())
        self._set_options(locals())

        self._set_budget(basis="global")

        return

    @staticmethod
    def get_dimensions() -> tuple:
        return ("nhru", "nsegment")  # nhru to be removed via exchange

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "hru_area",
            "hru_segment",  # to be removed via exchange
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
        return (
            "sroff_vol",
            "ssres_flow_vol",
            "gwres_flow_vol",
        )

    @staticmethod
    def get_init_values() -> dict:
        return {
            "channel_sroff_vol": nan,
            "channel_ssres_flow_vol": nan,
            "channel_gwres_flow_vol": nan,
            "channel_outflow_vol": nan,
            "seg_lateral_inflow": zero,
            "seg_upstream_inflow": zero,
            "seg_outflow": zero,
            "seg_outflow_substep": zero,
            "seg_stor_change": zero,
        }

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "channel_sroff_vol",
                "channel_ssres_flow_vol",
                "channel_gwres_flow_vol",
            ],
            "outputs": ["channel_outflow_vol"],
            "storage_changes": ["seg_stor_change"],
        }

    def get_outflow_mask(self):
        return self._outflow_mask

    @property
    def outflow_mask(self):
        return self._outflow_mask

    def _set_initial_conditions(self) -> None:
        # initialize channel segment storage
        self.seg_outflow[:] = self.segment_flow_init
        return

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

    # def _init_calc_method(self):
    #     if self._calc_method is None:
    #         self._calc_method = "numba"

    #     avail_methods = ["numpy", "numba", "fortran"]
    #     fortran_msg = ""
    #     if self._calc_method == "fortran" and not has_prmschannel_f:
    #         _ = avail_methods.remove("fortran")
    #         fortran_msg = "\n(Fortran not available as installed)\n"

    #     if self._calc_method.lower() not in avail_methods:
    #         msg = (
    #             f"Invalid calc_method={self._calc_method} for {self.name}. "
    #             f"{fortran_msg}"
    #             f"Setting calc_method to 'numba' for {self.name}"
    #         )
    #         warn(msg)
    #         self._calc_method = "numba"

    #     if self._calc_method.lower() == "numba":
    #         import numba as nb

    #         numba_msg = f"{self.name} jit compiling with numba "
    #         # this method can not be parallelized (? true?)
    #         print(numba_msg, flush=True)

    #         self._muskingum_mann = nb.njit(
    #             nb.types.UniTuple(nb.float64, 4)(
    #                 nb.int64,  # ihr
    #                 nb.float64,  # _seg_inflow0
    #                 nb.float64,  # seg_outflow
    #                 nb.float64,  # _inflow_ts
    #                 nb.float64,  # _outflow_ts
    #                 nb.int64,  # _tsi
    #                 nb.float64,  # _ts
    #                 nb.float64,  # _c0
    #                 nb.float64,  # _c1
    #                 nb.float64,  # _c2
    #             ),
    #             fastmath=True,
    #             parallel=False,
    #         )(_muskingum_mann_numpy)

    #     else:
    #         self._muskingum_mann = _muskingum_mann_numpy

    def _advance_variables(self) -> None:
        self._seg_inflow0[:] = self._seg_inflow

        self._seg_inflow[:] = zero
        self.seg_outflow[:] = zero
        self._inflow_ts[:] = zero

        self._seg_current_sum[:] = zero

        return

    def _init_graph(self) -> None:
        """Initialize internal variables from raw channel data"""

        # convert prms data to zero-based
        self._hru_segment = self.hru_segment - 1  # to be removed via exchange
        self._tosegment = self.tosegment - 1
        self._tosegment = self._tosegment.astype("int64")

        self._outflow_mask = np.full((len(self._tosegment)), False)

        connectivity = []
        for iseg in range(self.nsegment):
            tosegment = self._tosegment[iseg]
            if tosegment < 0:
                self._outflow_mask[iseg] = True
                continue
            connectivity.append(
                (
                    iseg,
                    tosegment,
                )
            )

        # use networkx to calculate the Directed Acyclic Graph
        self._graph = nx.DiGraph()
        self._graph.add_edges_from(connectivity)

        self._calculate_node_order()

        self._flow_nodes = []
        for ii in range(self.nsegment):
            self._flow_nodes.append(
                PRMSChannelFlowNode(
                    control=self.control,
                    tsi=self._tsi[ii],
                    ts=self._ts[ii],
                    c0=self._c0[ii],
                    c1=self._c1[ii],
                    c2=self._c2[ii],
                    graph_outflow=self._outflow_mask[ii],
                )
            )

        # initialize internal self_inflow variable
        for iseg in range(self.nsegment):
            jseg = self._tosegment[iseg]
            if jseg < 0:
                continue
            self._seg_inflow[jseg] = self.seg_outflow[iseg]

        return

    def _calculate_node_order(self):
        if self.nsegment > 1:
            node_order = list(nx.topological_sort(self._graph))
        else:
            node_order = [0]

        # if the domain contains links with no upstream or
        # downstream reaches, we just throw these back at the
        # top of the order since networkx wont handle such nonsense
        wh_mask_set = set(np.where(self._outflow_mask)[0])
        seg_ord_set = set(node_order)
        mask_not_seg_ord = list(wh_mask_set - seg_ord_set)
        if len(mask_not_seg_ord):
            node_order = mask_not_seg_ord + node_order
            # for pp in mask_not_seg_ord:
            #    assert (tosegment[pp] == -1) and (not pp in tosegment)

        self._node_order = np.array(node_order, dtype="int64")
        return

    def finalize_graph(self) -> None:
        self._graph_finalized = True

    def _calculate_lateral_inflows(self):
        # This could vary with timestep so leave here
        self._s_per_time = self.control.time_step_seconds

        self.seg_lateral_inflow[:] = 0.0
        for ihru in range(self.nhru):
            iseg = self._hru_segment[ihru]  # to be removed via exchange
            if iseg < 0:
                # This is bad, selective handling of fluxes is not cool,
                # mass is being discarded in a way that has to be coordinated
                # with other parts of the code.
                # This code shuold be removed evenutally.
                self.channel_sroff_vol[ihru] = zero
                self.channel_ssres_flow_vol[ihru] = zero
                self.channel_gwres_flow_vol[ihru] = zero
                continue

            else:
                self.channel_sroff_vol[ihru] = self.sroff_vol[ihru]
                self.channel_ssres_flow_vol[ihru] = self.ssres_flow_vol[ihru]
                self.channel_gwres_flow_vol[ihru] = self.gwres_flow_vol[ihru]

            # cubicfeet to cfs
            lateral_inflow = (
                self.channel_sroff_vol[ihru]
                + self.channel_ssres_flow_vol[ihru]
                + self.channel_gwres_flow_vol[ihru]
            ) / (self._s_per_time)
            self.seg_lateral_inflow[iseg] += lateral_inflow

        return

    def _calculate(self, simulation_time: float) -> None:
        if not self._graph_finalized:
            self.finalize_graph()

        # what in here should be the responsibility of the nodes?
        self._simulation_time = simulation_time

        self._calculate_lateral_inflows()

        # The upstream inflow and its accumulation dont belog to nodes
        self._seg_upstream_inflow_acc = self.seg_upstream_inflow * zero

        for ii in range(self.nsegment):
            self._flow_nodes[ii].advance()
            self._flow_nodes[ii].prepare_timestep()

        # Should probably pass time lengths to nodes
        for ihr in range(24):
            # This works because the first nodes calculated do
            # not have upstream reaches
            self.seg_upstream_inflow[:] = zero

            for jseg in self._node_order:
                self._flow_nodes[jseg].calculate_subtimestep(
                    ihr,
                    self.seg_upstream_inflow[jseg],
                    self.seg_lateral_inflow[jseg],
                    # jseg,
                )
                self.seg_outflow_substep[jseg] = self._flow_nodes[
                    jseg
                ].outflow_substep
                self._sum_outflows_to_inflow(jseg)

            # <
            self._seg_upstream_inflow_acc += self.seg_upstream_inflow

        for ii in range(self.nsegment):
            self._flow_nodes[ii].finalize_timestep()

        for ii in range(self.nsegment):
            self.seg_outflow[ii] = self._flow_nodes[ii].outflow
            self.seg_stor_change[ii] = self._flow_nodes[ii].storage_change

        # global mass balance term
        s_per_time = self.control.time_step_seconds
        self.channel_outflow_vol[:] = (
            np.where(self._outflow_mask, self.seg_outflow, zero)
        ) * s_per_time

        # <
        # PRMS uses the same variable for two different things
        # and we are continuing the bad practice here.
        # The next timestep apparently needs the average upstream
        # inflow from the previous timestep
        # It appears inconsistent to use substep flows except for the
        # initial substep...
        self.seg_upstream_inflow = self._seg_upstream_inflow_acc / 24.0

        return

    def _sum_outflows_to_inflow(self, jseg):
        if self._tosegment[jseg] >= 0:
            self.seg_upstream_inflow[
                self._tosegment[jseg]
            ] += self.seg_outflow_substep[jseg]

        return
