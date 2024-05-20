from typing import Literal

import numpy as np

from pywatershed.base.adapter import adaptable
from pywatershed.base.budget import Budget
from pywatershed.base.conservative_process import ConservativeProcess
from pywatershed.base.control import Control
from pywatershed.base.flow_graph import FlowNode, FlowNodeMaker
from pywatershed.constants import (
    cf_to_cm,
    cfs_to_cms,
    cm_to_cf,
    cms_to_cfs,
    nan,
    one,
    zero,
)
from pywatershed.parameters import Parameters
from pywatershed.utils.time_utils import datetime_epiweek

# Units note
# a variety of units are used but regression tests against the original code
# fail when composite conversion factors are used, apparently the order matters
# against the original.
# MCM = million cubic meters
# m3/day = m3pd
# m^3/second = m3ps
# m^3/week = m2pw

# daily conversions mass <-> volume
m3ps_to_MCM_day = 24 * 60 * 60 / 1.0e6
MCM_to_m3ps_day = 1.0 / m3ps_to_MCM_day

# hourly conversions mass <-> volume
m3ps_to_MCM_hour = 1 * 60 * 60 / 1.0e6
MCM_to_m3ps_hour = 1.0 / m3ps_to_MCM_hour


# constant
omega = 1.0 / 52.0


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
        budget_type: one of ["defer", None, "warn", "error"] with "defer" being
            the default and defering to control.options["budget_type"] when
            available. When control.options["budget_type"] is not avaiable,
            budget_type is set to "warn".
        verbose: Print extra information or not?
        load_n_time_batches: not-implemented

    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        lake_inflow: adaptable,
        budget_type: Literal["defer", None, "warn", "error"] = "defer",
        verbose: bool = False,
        load_n_time_batches: int = 1,
        io_in_cfs: bool = True,
    ) -> None:
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "Starfit"

        self._set_inputs(locals())
        self._set_options(locals())

        self._set_budget(basis="unit", ignore_nans=True)

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
                "lake_storage_change_flow_units",
            ],
        }

    @staticmethod
    def get_init_values() -> dict:
        return {
            "lake_storage": nan,
            "lake_storage_old": nan,
            "lake_storage_change": nan,
            "lake_storage_change_flow_units": nan,
            "lake_inflow": nan,
            "lake_release": nan,
            "lake_spill": nan,
            "lake_outflow": nan,
            "lake_availability_status": nan,
        }

    def _set_initial_conditions(self):
        # initialize storage when the appropriate time is reached
        # currently each reservoir can start at a different time.
        self.lake_storage[:] = nan
        self.lake_storage_old[:] = nan

        # JLM: this is sketchy seems like it should be a pre-process
        self.Obs_MEANFLOW_CUMECS = np.where(
            np.isnan(self.Obs_MEANFLOW_CUMECS),
            self.inflow_mean,
            self.Obs_MEANFLOW_CUMECS,
        )

        wh_initial_storage_nan = np.isnan(self.initial_storage)

        self.start_time = np.where(
            np.isnat(self.start_time),
            self.control.current_time,  # one day prior to start time
            self.start_time,
        )
        start_epiweeks = np.array(
            [datetime_epiweek(ss) for ss in self.start_time]
        )
        if len(wh_initial_storage_nan):
            min = min_nor(
                self.NORlo_max,
                self.NORlo_min,
                self.NORlo_alpha,
                self.NORlo_beta,
                self.NORlo_mu,
                omega,
                start_epiweeks,
            )
            max = max_nor(
                self.NORhi_max,
                self.NORhi_min,
                self.NORhi_alpha,
                self.NORhi_beta,
                self.NORhi_mu,
                omega,
                start_epiweeks,
            )
            pct_res_cap = (min + max) / 2
            nor_mean_cap = self.GRanD_CAP_MCM * pct_res_cap
            self._initial_storage = np.where(
                wh_initial_storage_nan, nor_mean_cap, self.initial_storage
            )

            # set lake_storage from initial_storage when start is not available
            self.lake_storage[:] = np.where(
                np.isnat(self.start_time), self.initial_storage, nan
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

        # For NaTs, time comparisons (other than !=) are False. If there is no
        # end_time (it is NaT), then this will be false.
        wh_after_end = np.where(self.control.current_time > self.end_time)
        self.lake_storage[wh_after_end] = nan
        self.lake_release[wh_after_end] = nan
        self.lake_spill[wh_after_end] = nan
        self.lake_availability_status[wh_after_end] = nan

        # For NaTs, time comparisons (other than !=) are False. If there is no
        # start_time (it is NaT), then this will be false. For those locations
        # without start_times, lake_storage was already set to initial_storage
        # in _set_initial_conditions.
        wh_start = np.where(self.control.current_time == self.start_time)
        self.lake_storage[wh_start] = self.initial_storage[wh_start]
        self.lake_storage_old[wh_start] = self.initial_storage[wh_start]

        if self._io_in_cfs:
            self.lake_inflow[:] *= cfs_to_cms
            self.lake_storage[:] *= cf_to_cm
            self.lake_storage_old[:] *= cf_to_cm

        # <
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
        )  # convert m3pd to m^3/s. THE NUMERICS ON THIS LINE ARE VERY
        # SENSITIVE TO THE ORDER OF THE DIVISION when reproducing the
        # test data. \_(`@`)_/
        self.lake_storage_change[:] = (
            self.lake_inflow - self.lake_release
        ) * m3ps_to_MCM_day
        # Note: no lake_storage calculation here

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
            ) * (
                MCM_to_m3ps_day
            )  # both terms in m3ps
            self.lake_release[wh_neg_storage] = np.maximum(
                potential_release,
                zero,
            )  # m3ps
            self.lake_storage_change[wh_neg_storage] = (
                self.lake_inflow[wh_neg_storage]
                - self.lake_release[wh_neg_storage]
            ) * m3ps_to_MCM_day

        self.lake_storage[:] = np.maximum(
            self.lake_storage + self.lake_storage_change, zero
        )  # MCM

        self.lake_spill[:] = nan
        wh_active = np.where(~np.isnan(self.lake_storage))
        self.lake_spill[wh_active] = zero
        wh_spill = np.where(self.lake_storage > self.GRanD_CAP_MCM)
        self.lake_spill[wh_spill] = (
            self.lake_storage[wh_spill] - self.GRanD_CAP_MCM[wh_spill]
        ) * MCM_to_m3ps_day
        self.lake_storage[wh_spill] = self.GRanD_CAP_MCM[wh_spill]
        self.lake_storage_change[:] = self.lake_storage - self.lake_storage_old
        self.lake_storage_change_flow_units[:] = (
            self.lake_storage_change * MCM_to_m3ps_day
        )
        self.lake_outflow[:] = self.lake_release + self.lake_spill

        if self._io_in_cfs:
            self.lake_inflow[:] *= cms_to_cfs
            self.lake_release[:] *= cms_to_cfs
            self.lake_spill[:] *= cms_to_cfs
            self.lake_outflow[:] *= cms_to_cfs
            self.lake_storage[:] *= cm_to_cf
            self.lake_storage_old[:] *= cm_to_cf
            self.lake_storage_change_flow_units[:] *= cms_to_cfs

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
        # MCM to m^3
        storage = lake_storage * 1.0e6
        capacity = GRanD_CAP_MCM * 1.0e6

        if not np.isfinite(grand_id).any():
            raise ValueError("Some non-finite grand_ids present")

        max_normal = max_nor(
            NORhi_max,
            NORhi_min,
            NORhi_alpha,
            NORhi_beta,
            NORhi_mu,
            omega,
            epiweek,
        )

        min_normal = min_nor(
            NORlo_max,
            NORlo_min,
            NORlo_alpha,
            NORlo_beta,
            NORlo_mu,
            omega,
            epiweek,
        )

        # TODO could make a better forecast?
        # why not use cumulative volume for the current epiweek, and only
        # extrapolate for the remainder of the week?
        # JLM: these lines are sensitive to the order of calculation when
        # JLM: comparing to the original code so wont use the constants here.
        forecasted_weekly_volume = 7.0 * lake_inflow * 24.0 * 60.0 * 60.0
        mean_weekly_volume = 7.0 * Obs_MEANFLOW_CUMECS * 24.0 * 60.0 * 60.0
        # forecasted_weekly_volume = lake_inflow * m3ps_to_m3pw
        # mean_weekly_volume = Obs_MEANFLOW_CUMECS * m3ps_to_m3pw

        standardized_inflow = (
            forecasted_weekly_volume / mean_weekly_volume
        ) - 1.0

        standardized_weekly_release = (
            Release_alpha1 * np.sin(2.0 * np.pi * omega * epiweek)
            + Release_alpha2 * np.sin(4.0 * np.pi * omega * epiweek)
            + Release_beta1 * np.cos(2.0 * np.pi * omega * epiweek)
            + Release_beta2 * np.cos(4.0 * np.pi * omega * epiweek)
        )

        # m3/week to m3/day
        release_min_vol = mean_weekly_volume * (1 + Release_min) / 7.0
        release_max_vol = mean_weekly_volume * (1 + Release_max) / 7.0

        availability_status = (100.0 * storage / capacity - min_normal) / (
            max_normal - min_normal
        )

        # m3/week to m3/day
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

        # m3/week to m3/day
        release_above_normal = (
            storage
            - (capacity * max_normal / 100.0)
            + forecasted_weekly_volume
        ) / 7.0

        # m3/week to m3/day
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


def max_nor(
    NORhi_max, NORhi_min, NORhi_alpha, NORhi_beta, NORhi_mu, omega, epiweek
):
    return np.minimum(
        NORhi_max,
        np.maximum(
            NORhi_min,
            NORhi_mu
            + NORhi_alpha * np.sin(2.0 * np.pi * omega * epiweek)
            + NORhi_beta * np.cos(2.0 * np.pi * omega * epiweek),
        ),
    )


def min_nor(
    NORlo_max, NORlo_min, NORlo_alpha, NORlo_beta, NORlo_mu, omega, epiweek
):
    return np.minimum(
        NORlo_max,
        np.maximum(
            NORlo_min,
            NORlo_mu
            + NORlo_alpha * np.sin(2.0 * np.pi * omega * epiweek)
            + NORlo_beta * np.cos(2.0 * np.pi * omega * epiweek),
        ),
    )


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
        io_in_cfs: bool = True,
        compute_daily: bool = False,
        budget_type: Literal["defer", None, "warn", "error"] = None,
    ):
        self.name = "StarfitFlowNode"
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
        # calc_method ignored currently

        self._io_in_cfs = io_in_cfs

        def nan1d():
            return np.zeros(1) * nan

        self._lake_inflow = nan1d()
        self._lake_inflow_cms = nan1d()
        self._lake_inflow_accum = nan1d()
        self._lake_inflow_sub = nan1d()
        self._lake_inflow_previous = nan1d()

        self._lake_outflow = nan1d()
        self._lake_outflow_accum = nan1d()
        self._lake_outflow_sub = nan1d()
        self._lake_outflow_sub_next = nan1d()

        self._lake_storage = nan1d()
        self._lake_storage_old = nan1d()
        self._lake_storage_accum = nan1d()

        self._lake_storage_sub = nan1d()
        self._lake_storage_old_sub = nan1d()

        self._lake_storage_change_sub = nan1d()
        self._lake_storage_change = nan1d()
        self._lake_storage_change_flow_units = nan1d()
        self._lake_storage_change_accum = nan1d()

        self._lake_release = nan1d()
        self._lake_release_sub = nan1d()
        self._lake_release_accum = nan1d()

        self._lake_spill = nan1d()
        self._lake_spill_sub = nan1d()
        self._lake_spill_accum = nan1d()

        self._lake_availability_status = nan1d()
        self._lake_availability_status_sub = nan1d()
        self._lake_availability_status_accum = nan1d()

        if np.isnan(self._Obs_MEANFLOW_CUMECS):
            self._Obs_MEANFLOW_CUMECS = self._inflow_mean

        wh_initial_storage_nan = np.isnan(self._initial_storage)
        if self._io_in_cfs:
            self._initial_storage = np.where(
                wh_initial_storage_nan,
                self._initial_storage,
                self._initial_storage * cf_to_cm,
            )

        start_time = np.where(
            np.isnat(self._start_time),
            self.control.current_time,  # one day prior to start time
            self._start_time,
        )

        start_epiweeks = np.array([datetime_epiweek(start_time[()])])

        if wh_initial_storage_nan:
            min = min_nor(
                self._NORlo_max,
                self._NORlo_min,
                self._NORlo_alpha,
                self._NORlo_beta,
                self._NORlo_mu,
                omega,
                start_epiweeks,
            )
            max = max_nor(
                self._NORhi_max,
                self._NORhi_min,
                self._NORhi_alpha,
                self._NORhi_beta,
                self._NORhi_mu,
                omega,
                start_epiweeks,
            )
            pct_res_cap = (min + max) / 2 / 100
            nor_mean_cap = self._GRanD_CAP_MCM * pct_res_cap
            # note the leading dunder (double underscore): "private reserve"
            self.__initial_storage = np.where(
                wh_initial_storage_nan, nor_mean_cap, self._initial_storage
            )
            # set lake_storage from initial_storage when start is not available
            self._lake_storage_sub[:] = np.where(
                np.isnat(self._start_time), self.__initial_storage, nan
            )

        else:
            self.__initial_storage = self._initial_storage

        self._budget_type = budget_type
        if self._budget_type == "defer":
            if "budget_type" in self.control.options.keys():
                self._budget_type = self.control.options["budget_type"]
            else:
                self._budget_type = "warn"
        if self._budget_type is not None:
            self.budget = Budget.from_storage_unit(
                self,
                time_unit="D",
                description=self.name,
                imbalance_fatal=(self._budget_type == "error"),
                basis="unit",
                ignore_nans=False,
            )
        else:
            self.budget = None

        self._compute_daily = compute_daily
        if self._compute_daily:
            self.calculate_subtimestep = self._calculate_subtimestep_daily
        else:
            self.calculate_subtimestep = self._calculate_subtimestep_hourly

        return

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "_lake_inflow",
            ],
            "outputs": [
                "_lake_release",
                "_lake_spill",
            ],
            "storage_changes": [
                "_lake_storage_change_flow_units",
            ],
        }

    def prepare_timestep(self):
        self._lake_inflow_accum[:] = np.array([zero])
        if self._compute_daily and self._io_in_cfs:
            self._lake_storage[:] *= cf_to_cm
            self._lake_storage_old[:] *= cf_to_cm
        else:
            self._lake_inflow_accum[:] = np.array([zero])
            self._lake_outflow_accum[:] = np.array([zero])
            self._lake_storage_accum[:] = np.array([zero])
            self._lake_storage_change_accum[:] = np.array([zero])
            self._lake_release_accum[:] = np.array([zero])
            self._lake_spill_accum[:] = np.array([zero])
            self._lake_availability_status_accum[:] = np.array([zero])

        return

    def finalize_timestep(self):
        if self._io_in_cfs:
            # self._lake_outflow_sub[:] converted above bc/ it is subtimestep
            self._lake_inflow[:] *= cms_to_cfs
            self._lake_release[:] *= cms_to_cfs
            self._lake_spill[:] *= cms_to_cfs
            self._lake_outflow[:] *= cms_to_cfs
            self._lake_storage[:] *= cm_to_cf
            self._lake_storage_old[:] *= cm_to_cf  # necessary
            self._lake_storage_change_flow_units[:] *= cms_to_cfs

        if self.budget is not None:
            self.budget.advance()
            self.budget.calculate()

        return

    def advance(self):
        self._lake_storage_change[:] = (
            self._lake_storage - self._lake_storage_old
        )
        self._lake_storage_old[:] = self._lake_storage
        return

    @property
    def outflow(self):
        return self._lake_outflow

    @property
    def outflow_substep(self):
        return self._lake_outflow_sub

    @property
    def storage_change(self):
        return self._lake_storage_change_flow_units

    @property
    def storage(self):
        return self._lake_storage

    @property
    def sink_source(self):
        return zero

    def _calculate_subtimestep_daily(
        self, isubstep, inflow_upstream, inflow_lateral
    ) -> None:
        # Here we only calculate outflows on the last subtimestep.
        # Ouflows on substeps are all the same flow rates, and
        # storages are only updated at the end of the timestep.

        # how to get the number of substeps in a timestep?  from control. this
        # might not even be a fixed number
        nsubsteps = 24

        # accumulate inflows
        self._lake_inflow_sub[:] = np.array([inflow_upstream + inflow_lateral])
        if self._io_in_cfs:
            self._lake_inflow_sub[:] *= cfs_to_cms
        self._lake_inflow_accum[:] += self._lake_inflow_sub

        # Special case for the very first subtimestep, we use the very first
        # inflow on the first subtimestep as representative of the mean
        # inflow of the previous day so we can calculate an average outflow
        # to use for the first timestep
        if self.control.itime_step == 0 and isubstep == 0:
            # this is the representative for the nonexistent previous day
            self._lake_inflow[:] = self._lake_inflow_accum
            # two previous, we'll assume the same
            self._lake_storage[:] = self._lake_storage_sub
            self._lake_storage_old[:] = self._lake_storage_sub
            self._lake_storage_change[:] = zero
        elif isubstep < (nsubsteps - 1):
            if isubstep == 0:
                # already in cfs
                self._lake_outflow_sub[:] = self._lake_outflow_sub_next
            return
        else:
            if self.control.itime_step == 0:
                # the end of the first timestep doesnt pass through advance()
                self._lake_storage_old[:] = self._lake_storage
            self._lake_inflow[:] = self._lake_inflow_accum / (isubstep + 1)
            if self._io_in_cfs:
                self._lake_outflow_sub[:] *= cfs_to_cms
                self._lake_release_sub[:] *= cfs_to_cms
                self._lake_spill_sub[:] *= cfs_to_cms
            self._lake_outflow[:] = self._lake_outflow_sub
            self._lake_release[:] = self._lake_release_sub
            self._lake_spill[:] = self._lake_spill_sub

            # calculate storage
            self._lake_storage_change_flow_units[:] = (
                self._lake_inflow - self._lake_outflow
            )
            # below, releases are never more than storage for the day, so this
            # should never be negative since inflows are never negative
            self._lake_storage_change[:] = (
                self._lake_storage_change_flow_units * m3ps_to_MCM_day
            )
            self._lake_storage[:] += self._lake_storage_change

        # <
        self._lake_spill_sub[:] = np.array([zero])
        if self._lake_storage > self._GRanD_CAP_MCM:
            self._lake_spill_sub[:] = (
                self._lake_storage - self._GRanD_CAP_MCM
            ) * MCM_to_m3ps_day
            # spill dosent affect the storage until the next timestep

        # now calculate the (avg) outflows for the next timestep
        (
            self._lake_release_sub[:],
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

        self._lake_release_sub *= m3ps_to_MCM_day / 24 / 60 / 60  # m3pd to MCM

        if (self._lake_storage - self._lake_release_sub) < zero:
            self._lake_release_sub[:] = self._lake_storage
        self._lake_release_sub[:] *= MCM_to_m3ps_day

        self._lake_outflow_sub_next[:] = (
            self._lake_release_sub + self._lake_spill_sub
        )

        if self._io_in_cfs:
            self._lake_outflow_sub[:] *= cms_to_cfs
            self._lake_outflow_sub_next[:] *= cms_to_cfs
            self._lake_release_sub[:] *= cms_to_cfs
            self._lake_spill_sub[:] *= cms_to_cfs

        if self.control.itime_step == 0 and isubstep == 0:
            self._lake_outflow_sub[:] = self._lake_outflow_sub_next
        return

    def _calculate_subtimestep_hourly(
        self, isubstep, inflow_upstream, inflow_lateral
    ) -> None:

        self._lake_inflow_sub[:] = np.array([inflow_upstream + inflow_lateral])
        if self._io_in_cfs:
            self._lake_inflow_sub[:] *= cfs_to_cms

        # <
        self._lake_storage_old_sub[:] = self._lake_storage_sub

        (
            self._lake_release_sub[:],
            self._lake_availability_status_sub[:],
        ) = Starfit._calc_istarf_release(
            epiweek=np.minimum(self.control.current_epiweek, 52),
            GRanD_CAP_MCM=self._GRanD_CAP_MCM,
            grand_id=self._grand_id,
            lake_inflow=self._lake_inflow_sub,
            lake_storage=self._lake_storage_sub,
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

        self._lake_release_sub[:] = (
            self._lake_release_sub / 24 / 60 / 60
        )  # m^3/s

        self._lake_storage_change_sub[:] = (
            self._lake_inflow_sub - self._lake_release_sub
        ) * m3ps_to_MCM_hour  # MCM: million cubic meters

        # can't release more than storage + inflow. This assumes zero
        # storage = deadpool which may not be accurate, but this situation
        # rarely occurs since STARFIT releases are already designed to keep
        # storage within NOR.
        if (self._lake_storage_sub + self._lake_storage_change_sub) < zero:
            potential_release = (
                self._lake_release_sub
                + (self._lake_storage_sub + self._lake_storage_change_sub)
                * MCM_to_m3ps_hour
            )
            self._lake_release_sub[:] = np.maximum(
                potential_release,
                potential_release * zero,
            )  # m^3/s
            self._lake_storage_change_sub[:] = (
                self._lake_inflow_sub - self._lake_release_sub
            ) * m3ps_to_MCM_hour  # MCM: million cubic meters

        self._lake_storage_sub[:] = np.maximum(
            self._lake_storage_sub + self._lake_storage_change_sub,
            zero,
        )  # MCM

        self._lake_spill_sub[:] = nan
        if ~np.isnan(self._lake_storage_sub):
            self._lake_spill_sub[:] = zero

        if self._lake_storage_sub > self._GRanD_CAP_MCM:
            self._lake_spill_sub[:] = (
                self._lake_storage_sub - self._GRanD_CAP_MCM
            ) * MCM_to_m3ps_hour
            self._lake_storage_sub[:] = self._GRanD_CAP_MCM

        self._lake_storage_change_sub[:] = (
            self._lake_storage_sub - self._lake_storage_old_sub
        )

        # subtimestep to timestep calculations
        # m^3/s
        self._lake_inflow_accum[:] += self._lake_inflow_sub
        self._lake_inflow[:] = self._lake_inflow_accum / (isubstep + 1)

        self._lake_outflow_sub[:] = (
            self._lake_release_sub + self._lake_spill_sub
        )
        self._lake_outflow_accum[:] += self._lake_outflow_sub
        self._lake_outflow[:] = self._lake_outflow_accum / (isubstep + 1)

        # these variables arent currently retrieved and output by FlowGraph
        # but probably soon
        self._lake_release_accum[:] += self._lake_release_sub
        self._lake_release[:] = self._lake_release_accum / (isubstep + 1)

        self._lake_spill_accum[:] += self._lake_spill_sub
        self._lake_spill[:] = self._lake_spill_accum / (isubstep + 1)

        self._lake_availability_status_accum[
            :
        ] += self._lake_availability_status_sub
        self._lake_availability_status[:] = (
            self._lake_availability_status_accum / (isubstep + 1)
        )
        self._lake_storage_accum[:] += self._lake_storage_sub
        self._lake_storage[:] = self._lake_storage_accum / (isubstep + 1)

        # million volume units
        # self._storage_change_subtimestep
        self._lake_storage_change_accum[:] += self._lake_storage_change_sub
        self._lake_storage_change[:] = self._lake_storage_change_accum / (
            isubstep + 1
        )
        self._lake_storage_change_flow_units[:] = (
            self._lake_storage_change * MCM_to_m3ps_hour
        )

        if self._io_in_cfs:
            self._lake_outflow_sub[:] *= cms_to_cfs

        return


class StarfitFlowNodeMaker(FlowNodeMaker):
    def __init__(
        self,
        discretization: Parameters,
        parameters: Parameters,
        calc_method: Literal["numba", "numpy"] = None,
        io_in_cfs: bool = True,
        verbose: bool = None,
        compute_daily: bool = True,
        budget_type: Literal["defer", None, "warn", "error"] = None,
    ) -> None:
        self.name = "StarfitFlowNodeMaker"
        self._calc_method = calc_method
        self._io_in_cfs = io_in_cfs
        self._compute_daily = compute_daily
        self._budget_type = budget_type

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
            io_in_cfs=self._io_in_cfs,
            compute_daily=self._compute_daily,
            budget_type=self._budget_type,
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
