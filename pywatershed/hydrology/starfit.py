from typing import Literal

import numpy as np

from pywatershed.base.conservative_process import ConservativeProcess

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import nan, one, zero
from ..parameters import Parameters


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
