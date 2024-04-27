from typing import Literal

import numpy as np

from pywatershed.base.adapter import adaptable
from pywatershed.base.conservative_process import ConservativeProcess
from pywatershed.base.control import Control
from pywatershed.base.flow_graph import FlowNode, FlowNodeMaker
from pywatershed.constants import nan, one, zero
from pywatershed.parameters import Parameters

# MCM is million cubic meters


class Starfit(ConservativeProcess):
    """starfit: Storage Targets And Release Function Inference Tool

    Sean W.D. Turner, Jennie Clarice Steyaert, Laura Condon, Nathalie Voisin,
    Water storage and release policies for all large reservoirs of conterminous
    United States, Journal of Hydrology, Volume 603, Part A, 2021, 126843,
    ISSN 0022-1694, https://doi.org/10.1016/j.jhydrol.2021.126843.

    https://github.com/IMMM-SFA/starfit

    Adapted from STARFIT implementation in the MOSART-WM model:
    https://github.com/IMMM-SFA/mosartwmpy/blob/main/mosartwmpy/reservoirs/istarf.py

    Noah Knowles (USGS) and James McCreight (UCAR/USGS)

    Args:
        control: a Control object
        discretization: a discretization of class Parameters
        parameters: a parameter object of class Parameters
        lake_inflow: Daily lake inflow
        budget_type: one of [None, "warn", "error"]
        verbose: Print extra information or not?
        load_n_time_batches: not-implemented

    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        lake_inflow: adaptable,
        budget_type: Literal[None, "warn", "error"] = None,
        verbose: bool = False,
        load_n_time_batches: int = 1,
    ) -> None:
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "Starfit"

        self._set_inputs(locals())
        self._set_options(locals())

        self._set_budget(budget_type)
        return

    @staticmethod
    def get_dimensions() -> dict:
        return ("nreservoirs",)

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "grand_id",
            "initial_storage",
            "start_time",
            "end_time",
            "inflow_mean",
            "NORhi_min",
            "NORhi_max",
            "NORhi_alpha",
            "NORhi_beta",
            "NORhi_mu",
            "NORlo_min",
            "NORlo_max",
            "NORlo_alpha",
            "NORlo_beta",
            "NORlo_mu",
            "Release_min",
            "Release_max",
            "Release_alpha1",
            "Release_alpha2",
            "Release_beta1",
            "Release_beta2",
            "Release_p1",
            "Release_p2",
            "Release_c",
            "GRanD_CAP_MCM",
            "Obs_MEANFLOW_CUMECS",
        )

    @staticmethod
    def get_inputs() -> tuple:
        return ("lake_inflow",)

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "lake_inflow",
            ],
            "outputs": [
                "lake_release",
                "lake_spill",
            ],
            "storage_changes": [
                "lake_storage_change",
            ],
        }

    @staticmethod
    def get_init_values() -> dict:
        return {
            "lake_storage": nan,
            "lake_storage_old": nan,
            "lake_storage_change": nan,
            "lake_inflow": nan,
            "lake_release": nan,
            "lake_spill": nan,
            "lake_availability_status": nan,
        }

    def _set_initial_conditions(self):
        # initialize storage when the appropriate time is reached
        # currently each reservoir can start at a different time.
        self.lake_storage[:] = nan
        self.lake_storage_old[:] = nan

        # HMM, this is sketchy seems like it should be a pre-process
        self.Obs_MEANFLOW_CUMECS = np.where(
            np.isnan(self.Obs_MEANFLOW_CUMECS),
            self.inflow_mean,
            self.Obs_MEANFLOW_CUMECS,
        )

        return

    def _advance_variables(self) -> None:
        """Advance the starfit variables
        Returns:
            None
        """
        self.lake_storage_change[:] = self.lake_storage - self.lake_storage_old
        self.lake_storage_old[:] = self.lake_storage
        return

    def _calculate(self, simulation_time):
        """Calculate starfit for a time step (vectorized)

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

        wh_after_end = np.where(self.control.current_time > self.end_time)
        self.lake_storage[wh_after_end] = nan
        self.lake_release[wh_after_end] = nan
        self.lake_spill[wh_after_end] = nan
        self.lake_availability_status[wh_after_end] = nan

        wh_start = np.where(self.control.current_time == self.start_time)
        self.lake_storage[wh_start] = self.initial_storage[wh_start]
        self.lake_storage_old[wh_start] = self.initial_storage[wh_start]

        (
            self.lake_release[:],
            self.lake_availability_status[:],
        ) = self._calc_istarf_release(
            epiweek=np.minimum(self.control.current_epiweek, 52),
            GRanD_CAP_MCM=self.GRanD_CAP_MCM,
            grand_id=self.grand_id,
            lake_inflow=self.lake_inflow,
            lake_storage=self.lake_storage,
            NORhi_alpha=self.NORhi_alpha,
            NORhi_beta=self.NORhi_beta,
            NORhi_max=self.NORhi_max,
            NORhi_min=self.NORhi_min,
            NORhi_mu=self.NORhi_mu,
            NORlo_alpha=self.NORlo_alpha,
            NORlo_beta=self.NORlo_beta,
            NORlo_max=self.NORlo_max,
            NORlo_min=self.NORlo_min,
            NORlo_mu=self.NORlo_mu,
            Obs_MEANFLOW_CUMECS=self.Obs_MEANFLOW_CUMECS,
            Release_alpha1=self.Release_alpha1,
            Release_alpha2=self.Release_alpha2,
            Release_beta1=self.Release_beta1,
            Release_beta2=self.Release_beta2,
            Release_c=self.Release_c,
            Release_max=self.Release_max,
            Release_min=self.Release_min,
            Release_p1=self.Release_p1,
            Release_p2=self.Release_p2,
        )  # output in m^3/d

        self.lake_release[:] = (
            self.lake_release / 24 / 60 / 60
        )  # convert to m^3/s

        self.lake_storage_change[:] = (
            (self.lake_inflow - self.lake_release) * 24 * 60 * 60 / 1.0e6
        )  # conv to MCM

        # can't release more than storage + inflow. This assumes zero
        # storage = deadpool which may not be accurate, but this situation
        # rarely occurs since STARFIT releases are already designed to keep
        # storage within NOR.
        wh_neg_storage = np.where(
            (self.lake_storage + self.lake_storage_change) < zero
        )
        if len(wh_neg_storage):
            potential_release = self.lake_release[wh_neg_storage] + (
                self.lake_storage[wh_neg_storage]
                + self.lake_storage_change[wh_neg_storage]
            ) * (1.0e6 / 24 / 60 / 60)
            self.lake_release[wh_neg_storage] = np.maximum(
                potential_release,
                potential_release * zero,
            )
            self.lake_storage_change[wh_neg_storage] = (
                (
                    self.lake_inflow[wh_neg_storage]
                    - self.lake_release[wh_neg_storage]
                )
                * 24
                * 60
                * 60
                / 1.0e6
            )

        self.lake_storage[:] = np.maximum(
            self.lake_storage + self.lake_storage_change, zero
        )

        self.lake_spill[:] = nan
        wh_active = np.where(~np.isnan(self.lake_storage))
        self.lake_spill[wh_active] = zero
        wh_spill = np.where(self.lake_storage > self.GRanD_CAP_MCM)
        self.lake_spill[wh_spill] = (
            (self.lake_storage[wh_spill] - self.GRanD_CAP_MCM[wh_spill])
            * 1.0e6
            / 24
            / 60
            / 60
        )
        self.lake_storage[wh_spill] = self.GRanD_CAP_MCM[wh_spill]

        return

    @staticmethod
    def _calc_istarf_release(
        epiweek,
        GRanD_CAP_MCM,
        grand_id,
        lake_inflow,
        lake_storage,
        NORhi_alpha,
        NORhi_beta,
        NORhi_max,
        NORhi_min,
        NORhi_mu,
        NORlo_alpha,
        NORlo_beta,
        NORlo_max,
        NORlo_min,
        NORlo_mu,
        Obs_MEANFLOW_CUMECS,
        Release_alpha1,
        Release_alpha2,
        Release_beta1,
        Release_beta2,
        Release_c,
        Release_max,
        Release_min,
        Release_p1,
        Release_p2,
    ):
        # input is in MCM, this function needs m^3
        storage = lake_storage * 1.0e6
        capacity = GRanD_CAP_MCM * 1.0e6

        # constant
        omega = 1.0 / 52.0

        if not np.isfinite(grand_id).any():
            raise ValueError("Some non-finite grand_ids present")

        max_normal = np.minimum(
            NORhi_max,
            np.maximum(
                NORhi_min,
                NORhi_mu
                + NORhi_alpha * np.sin(2.0 * np.pi * omega * epiweek)
                + NORhi_beta * np.cos(2.0 * np.pi * omega * epiweek),
            ),
        )

        min_normal = np.minimum(
            NORlo_max,
            np.maximum(
                NORlo_min,
                NORlo_mu
                + NORlo_alpha * np.sin(2.0 * np.pi * omega * epiweek)
                + NORlo_beta * np.cos(2.0 * np.pi * omega * epiweek),
            ),
        )

        # TODO could make a better forecast?
        # why not use cumulative volume for the current epiweek, and only
        # extrapolate for the remainder of the week?
        forecasted_weekly_volume = 7.0 * lake_inflow * 24.0 * 60.0 * 60.0
        mean_weekly_volume = 7.0 * Obs_MEANFLOW_CUMECS * 24.0 * 60.0 * 60.0

        standardized_inflow = (
            forecasted_weekly_volume / mean_weekly_volume
        ) - 1.0

        standardized_weekly_release = (
            Release_alpha1 * np.sin(2.0 * np.pi * omega * epiweek)
            + Release_alpha2 * np.sin(4.0 * np.pi * omega * epiweek)
            + Release_beta1 * np.cos(2.0 * np.pi * omega * epiweek)
            + Release_beta2 * np.cos(4.0 * np.pi * omega * epiweek)
        )

        release_min_vol = mean_weekly_volume * (1 + Release_min) / 7.0
        release_max_vol = mean_weekly_volume * (1 + Release_max) / 7.0

        availability_status = (100.0 * storage / capacity - min_normal) / (
            max_normal - min_normal
        )

        release = (
            mean_weekly_volume
            * (
                1
                + (
                    standardized_weekly_release
                    + Release_c
                    + Release_p1 * availability_status
                    + Release_p2 * standardized_inflow
                )
            )
        ) / 7.0

        release_above_normal = (
            storage
            - (capacity * max_normal / 100.0)
            + forecasted_weekly_volume
        ) / 7.0

        release_below_normal = (
            storage
            - (capacity * min_normal / 100.0)
            + forecasted_weekly_volume
        ) / 7.0  # NK: The first part of this sum will be negative.

        release = np.where(
            availability_status > one, release_above_normal, release
        )
        release = np.where(
            availability_status < zero, release_below_normal, release
        )
        release = np.where(release < release_min_vol, release_min_vol, release)
        release = np.where(release > release_max_vol, release_max_vol, release)

        # storage update and boundaries are enforced during the regulation step
        return release, availability_status


class StarfitFlowNode(FlowNode):
    def __init__(
        self,
        control: Control,
        grand_id: np.int64,
        initial_storage: np.float64,
        start_time: np.datetime64,
        end_time: np.datetime64,
        inflow_mean: np.float64,
        NORhi_min: np.float64,
        NORhi_max: np.float64,
        NORhi_alpha: np.float64,
        NORhi_beta: np.float64,
        NORhi_mu: np.float64,
        NORlo_min: np.float64,
        NORlo_max: np.float64,
        NORlo_alpha: np.float64,
        NORlo_beta: np.float64,
        NORlo_mu: np.float64,
        Release_min: np.float64,
        Release_max: np.float64,
        Release_alpha1: np.float64,
        Release_alpha2: np.float64,
        Release_beta1: np.float64,
        Release_beta2: np.float64,
        Release_p1: np.float64,
        Release_p2: np.float64,
        Release_c: np.float64,
        GRanD_CAP_MCM: np.float64,
        Obs_MEANFLOW_CUMECS: np.float64,
        calc_method: Literal["numba", "numpy"] = None,
    ):
        self.control = control

        self._grand_id = grand_id
        self._initial_storage = initial_storage
        self._start_time = start_time
        self._end_time = end_time
        self._inflow_mean = inflow_mean
        self._NORhi_min = NORhi_min
        self._NORhi_max = NORhi_max
        self._NORhi_alpha = NORhi_alpha
        self._NORhi_beta = NORhi_beta
        self._NORhi_mu = NORhi_mu
        self._NORlo_min = NORlo_min
        self._NORlo_max = NORlo_max
        self._NORlo_alpha = NORlo_alpha
        self._NORlo_beta = NORlo_beta
        self._NORlo_mu = NORlo_mu
        self._Release_min = Release_min
        self._Release_max = Release_max
        self._Release_alpha1 = Release_alpha1
        self._Release_alpha2 = Release_alpha2
        self._Release_beta1 = Release_beta1
        self._Release_beta2 = Release_beta2
        self._Release_p1 = Release_p1
        self._Release_p2 = Release_p2
        self._Release_c = Release_c
        self._GRanD_CAP_MCM = GRanD_CAP_MCM
        self._Obs_MEANFLOW_CUMECS = Obs_MEANFLOW_CUMECS

        self._lake_storage = np.zeros(1) * nan
        self._lake_storage_old = np.zeros(1) * nan
        self._lake_storage_change = np.zeros(1) * nan
        # self._lake_inflow = np.zeros(1) * nan
        self._lake_release = np.zeros(1) * nan
        self._lake_spill = np.zeros(1) * nan
        self._lake_availability_status = np.zeros(1) * nan

        if np.isnan(self._Obs_MEANFLOW_CUMECS):
            self._Obs_MEANFLOW_CUMECS = self._inflow_mean

        return

    def prepare_timestep(self):
        return

    def calculate_subtimestep(
        self, ihr, inflow_upstream, inflow_lateral
    ) -> None:
        # this does not have subtimesteps, the outflow is the same on
        # all of them, so dont re-solve
        if ihr > 0:
            return
        self._lake_inflow = inflow_upstream + inflow_lateral
        if (self.control.current_time > self._end_time) or (
            self.control.current_time < self._start_time
        ):
            self._lake_storage[:] = np.zeros(1) * nan
            self._lake_release[:] = np.zeros(1) * nan
            self._lake_spill[:] = np.zeros(1) * nan
            self._lake_availability_status[:] = np.zeros(1) * nan
            return

        if self.control.current_time == self._start_time:
            self._lake_storage[:] = self._initial_storage
            self._lake_storage_old[:] = self._initial_storage

        (
            self._lake_release[:],
            self._lake_availability_status[:],
        ) = Starfit._calc_istarf_release(
            epiweek=np.minimum(self.control.current_epiweek, 52),
            GRanD_CAP_MCM=self._GRanD_CAP_MCM,
            grand_id=self._grand_id,
            lake_inflow=self._lake_inflow,
            lake_storage=self._lake_storage,
            NORhi_alpha=self._NORhi_alpha,
            NORhi_beta=self._NORhi_beta,
            NORhi_max=self._NORhi_max,
            NORhi_min=self._NORhi_min,
            NORhi_mu=self._NORhi_mu,
            NORlo_alpha=self._NORlo_alpha,
            NORlo_beta=self._NORlo_beta,
            NORlo_max=self._NORlo_max,
            NORlo_min=self._NORlo_min,
            NORlo_mu=self._NORlo_mu,
            Obs_MEANFLOW_CUMECS=self._Obs_MEANFLOW_CUMECS,
            Release_alpha1=self._Release_alpha1,
            Release_alpha2=self._Release_alpha2,
            Release_beta1=self._Release_beta1,
            Release_beta2=self._Release_beta2,
            Release_c=self._Release_c,
            Release_max=self._Release_max,
            Release_min=self._Release_min,
            Release_p1=self._Release_p1,
            Release_p2=self._Release_p2,
        )  # output in m^3/d

        self._lake_release[:] = self._lake_release / 24 / 60 / 60  # m^3/s

        self._lake_storage_change[:] = (
            (self._lake_inflow - self._lake_release) * 24 * 60 * 60 / 1.0e6
        )  # MCM: million cubic meters

        # can't release more than storage + inflow. This assumes zero
        # storage = deadpool which may not be accurate, but this situation
        # rarely occurs since STARFIT releases are already designed to keep
        # storage within NOR.
        if (self._lake_storage + self._lake_storage_change) < zero:
            potential_release = self._lake_release + (
                self._lake_storage + self._lake_storage_change
            ) * (1.0e6 / 24 / 60 / 60)  # m^3/s
            self._lake_release[:] = np.maximum(
                potential_release,
                potential_release * zero,
            )  # m^3/s
            self._lake_storage_change[:] = (
                (self._lake_inflow - self._lake_release) * 24 * 60 * 60 / 1.0e6
            )  # MCM: million cubic meters

        self._lake_storage[:] = np.maximum(
            self._lake_storage + self._lake_storage_change, zero
        )  # MCM

        self._lake_spill[:] = nan
        if ~np.isnan(self._lake_storage):
            self._lake_spill[:] = zero

        if self._lake_storage > self._GRanD_CAP_MCM:
            self._lake_spill[:] = (
                (self._lake_storage - self._GRanD_CAP_MCM)
                * 1.0e6
                / 24
                / 60
                / 60
            )  # m^3/s
            self._lake_storage[:] = self._GRanD_CAP_MCM  # MCM

        # m^3/s
        self._outflow = self._lake_release + self._lake_spill

        return

    def finalize_timestep(self):
        return

    def advance(self):
        self._lake_storage_change[:] = (
            self._lake_storage - self._lake_storage_old
        )
        self._lake_storage_old[:] = self._lake_storage
        return

    @property
    def outflow(self):
        return self._outflow

    @property
    def storage_change(self):
        # should this copy?
        return self._lake_storage_change

    @property
    def storage(self):
        # should this copy?
        return self._lake_storage

    @property
    def sink_source(self):
        return zero


class StarfitFlowNodeMaker(FlowNodeMaker):
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

        return

    def get_node(self, control, index):
        return StarfitFlowNode(
            control=control,
            grand_id=self.grand_id[index],
            initial_storage=self.initial_storage[index],
            start_time=self.start_time[index],
            end_time=self.end_time[index],
            inflow_mean=self.inflow_mean[index],
            NORhi_min=self.NORhi_min[index],
            NORhi_max=self.NORhi_max[index],
            NORhi_alpha=self.NORhi_alpha[index],
            NORhi_beta=self.NORhi_beta[index],
            NORhi_mu=self.NORhi_mu[index],
            NORlo_min=self.NORlo_min[index],
            NORlo_max=self.NORlo_max[index],
            NORlo_alpha=self.NORlo_alpha[index],
            NORlo_beta=self.NORlo_beta[index],
            NORlo_mu=self.NORlo_mu[index],
            Release_min=self.Release_min[index],
            Release_max=self.Release_max[index],
            Release_alpha1=self.Release_alpha1[index],
            Release_alpha2=self.Release_alpha2[index],
            Release_beta1=self.Release_beta1[index],
            Release_beta2=self.Release_beta2[index],
            Release_p1=self.Release_p1[index],
            Release_p2=self.Release_p2[index],
            Release_c=self.Release_c[index],
            GRanD_CAP_MCM=self.GRanD_CAP_MCM[index],
            Obs_MEANFLOW_CUMECS=self.Obs_MEANFLOW_CUMECS[index],
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

        self.nreservoirs = len(self.grand_id)
        return

    @staticmethod
    def get_dimensions() -> dict:
        return ("nreservoirs",)

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "grand_id",
            "initial_storage",
            "start_time",
            "end_time",
            "inflow_mean",
            "NORhi_min",
            "NORhi_max",
            "NORhi_alpha",
            "NORhi_beta",
            "NORhi_mu",
            "NORlo_min",
            "NORlo_max",
            "NORlo_alpha",
            "NORlo_beta",
            "NORlo_mu",
            "Release_min",
            "Release_max",
            "Release_alpha1",
            "Release_alpha2",
            "Release_beta1",
            "Release_beta2",
            "Release_p1",
            "Release_p2",
            "Release_c",
            "GRanD_CAP_MCM",
            "Obs_MEANFLOW_CUMECS",
        )
