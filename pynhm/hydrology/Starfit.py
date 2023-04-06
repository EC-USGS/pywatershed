import numpy as np

from pynhm.base.storageUnit import StorageUnit

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import nan, one, zero


class Starfit(StorageUnit):
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
    """

    def __init__(
        self,
        control: Control,
        lake_inflow: adaptable,
        budget_type: str = None,
        calc_method: str = None,
        verbose: bool = False,
        load_n_time_batches: int = 1,
    ) -> "Starfit":
        super().__init__(
            control=control,
            verbose=verbose,
            load_n_time_batches=load_n_time_batches,
        )
        self.name = "Starfit"

        self._calc_method = str(calc_method)

        self._set_inputs(locals())
        self._set_budget(budget_type)
        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get groundwater reservoir parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nreservoirs",
            "nhru",  # to remove
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
        """Get groundwater reservoir input variables

        Returns:
            variables: input variables

        """
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
        """Get groundwater initial values

        Returns:
            dict: initial values for named variables
        """
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
        """Advance the groundwater reservoir variables
        Returns:
            None
        """
        self.lake_storage_change[:] = self.lake_storage - self.lake_storage_old
        self.lake_storage_old[:] = self.lake_storage
        return

    def _calculate(self, simulation_time):
        """Calculate groundwater reservoir terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

        for ires in range(self.nreservoirs):
            if self.control.current_time > self.end_time[ires]:
                self.lake_storage[ires] = nan
                self.lake_release[ires] = nan
                self.lake_spill[ires] = nan
                self.lake_availability_status[ires] = nan
                continue

            if self.control.current_time < self.start_time[ires]:
                continue

            if self.control.current_time == self.start_time[ires]:
                self.lake_storage[ires] = self.initial_storage[ires]
                self.lake_storage_old[ires] = self.initial_storage[ires]

        if self._calc_method.lower() in ["none", "numpy"]:
            (
                self.lake_release,
                self.lake_availability_status,
            ) = self._calculate_numpy(
                # use kws, sort by kw?
                np.minimum(self.control.current_epiweek, 52),
                self.grand_id,
                self.NORhi_min,
                self.NORhi_max,
                self.NORhi_alpha,
                self.NORhi_beta,
                self.NORhi_mu,
                self.NORlo_min,
                self.NORlo_max,
                self.NORlo_alpha,
                self.NORlo_beta,
                self.NORlo_mu,
                self.Release_min,
                self.Release_max,
                self.Release_alpha1,
                self.Release_alpha2,
                self.Release_beta1,
                self.Release_beta2,
                self.Release_p1,
                self.Release_p2,
                self.Release_c,
                self.GRanD_CAP_MCM,
                self.Obs_MEANFLOW_CUMECS,
                self.lake_storage,
                self.lake_inflow,
            )  # output in m^3/d

        else:
            msg = f"Invalid calc_method={self._calc_method} for {self.name}"
            raise ValueError(msg)

        self.lake_release[:] = (
            self.lake_release / 24 / 60 / 60
        )  # convert to m^3/s

        # update storage (dS=I-R)
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
    def _calculate_numpy(
        epiweek,
        reservoir_id,
        upper_min,
        upper_max,
        upper_alpha,
        upper_beta,
        upper_mu,
        lower_min,
        lower_max,
        lower_alpha,
        lower_beta,
        lower_mu,
        release_min_parameter,
        release_max_parameter,
        release_alpha_one,
        release_alpha_two,
        release_beta_one,
        release_beta_two,
        release_p_one,
        release_p_two,
        release_c,
        capacity_MCM,
        inflow_mean,
        storage_MCM,
        inflow,
    ):
        # input is in MCM, this function needs m^3
        storage = storage_MCM * 1.0e6
        capacity = capacity_MCM * 1.0e6

        # constant
        omega = 1.0 / 52.0

        if not np.isfinite(reservoir_id).any():
            raise ValueError("Some non-finite reservoir_ids present")

        max_normal = np.minimum(
            upper_max,
            np.maximum(
                upper_min,
                upper_mu
                + upper_alpha * np.sin(2.0 * np.pi * omega * epiweek)
                + upper_beta * np.cos(2.0 * np.pi * omega * epiweek),
            ),
        )

        min_normal = np.minimum(
            lower_max,
            np.maximum(
                lower_min,
                lower_mu
                + lower_alpha * np.sin(2.0 * np.pi * omega * epiweek)
                + lower_beta * np.cos(2.0 * np.pi * omega * epiweek),
            ),
        )

        # TODO could make a better forecast?
        # why not use cumulative volume for the current epiweek, and only
        # extrapolate for the remainder of the week?
        forecasted_weekly_volume = 7.0 * inflow * 24.0 * 60.0 * 60.0
        mean_weekly_volume = 7.0 * inflow_mean * 24.0 * 60.0 * 60.0

        standardized_inflow = (
            forecasted_weekly_volume / mean_weekly_volume
        ) - 1.0

        standardized_weekly_release = (
            release_alpha_one * np.sin(2.0 * np.pi * omega * epiweek)
            + release_alpha_two * np.sin(4.0 * np.pi * omega * epiweek)
            + release_beta_one * np.cos(2.0 * np.pi * omega * epiweek)
            + release_beta_two * np.cos(4.0 * np.pi * omega * epiweek)
        )

        release_min = mean_weekly_volume * (1 + release_min_parameter) / 7.0
        release_max = mean_weekly_volume * (1 + release_max_parameter) / 7.0

        availability_status = (100.0 * storage / capacity - min_normal) / (
            max_normal - min_normal
        )

        # # above normal
        # if availability_status > 1:
        #     release = (
        #         storage
        #         - (capacity * max_normal / 100.0)
        #         + forecasted_weekly_volume
        #     ) / 7.0

        # # below normal
        # elif availability_status < 0:
        #     release = (
        #         storage
        #         - (capacity * min_normal / 100.0)
        #         + forecasted_weekly_volume
        #     ) / 7.0  # NK: The first part of this sum will be negative.

        # # within normal
        # else:
        #     release = (
        #         mean_weekly_volume
        #         * (
        #             1
        #             + (
        #                 standardized_weekly_release
        #                 + release_c
        #                 + release_p_one * availability_status
        #                 + release_p_two * standardized_inflow
        #             )
        #         )
        #     ) / 7.0

        # # enforce boundaries on release
        # if release < release_min:
        #     # shouldn't release less than min
        #     release = release_min
        # elif release > release_max:
        #     # shouldn't release more than max
        #     release = release_max

        release = (
            mean_weekly_volume
            * (
                1
                + (
                    standardized_weekly_release
                    + release_c
                    + release_p_one * availability_status
                    + release_p_two * standardized_inflow
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
        release = np.where(release < release_min, release_min, release)
        release = np.where(release > release_max, release_max, release)

        # storage update and boundaries are enforced during the regulation step
        return release, availability_status
