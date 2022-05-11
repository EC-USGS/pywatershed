from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control

from ..constants import one, nan, zero, HruType

adaptable = Union[str, np.ndarray, Adapter]


class PRMSSoilzone(StorageUnit):
    """PRMS soil zone

    Args:
        params: parameter object
        atm: atmosphere object

    """

    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        dprst_evap_hru: adaptable,
        dprst_seep_hru: adaptable,
        hru_impervevap: adaptable,
        hru_intcpevap: adaptable,
        infil: adaptable,
        sroff: adaptable,
        potet: adaptable,
        transp_on: adaptable,
        snow_evap: adaptable,
        snowcov_area: adaptable,
        verbose: bool = False,
    ) -> "PRMSSoilzone":

        self.name = "PRMSSoilzone"
        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
            subclass_name=self.name,
        )

        # Adapt every input
        self._input_variables_dict = {}
        for ii in self.inputs:
            self._input_variables_dict[ii] = adapter_factory(locals()[ii], ii)

        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get soil zone parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nhru",
            "cov_type",
            "fastcoef_lin",
            "fastcoef_sq",
            "hru_area",
            "hru_percent_imperv",
            "hru_type",
            "pref_flow_den",
            "sat_threshold",
            "slowcoef_lin",
            "slowcoef_sq",
            "soil_moist_max",
            "soil_moist_init_frac",
            "soil_rechr_init_frac",
            "soil_rechr_max_frac",
            "soil_type",
            "soil2gw_max",
            "ssr2gw_exp",
            "ssr2gw_rate",
            "ssstor_init_frac",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get soil zone input variables

        Returns:
            variables: input variables

        """
        return (
            "snowcov_area",
            "dprst_evap_hru",  ## JLM ?? needs this stuff to calculate evap?
            "dprst_seep_hru",
            "hru_impervevap",  ## JLM ??
            "hru_intcpevap",  # JLM ???
            "infil",
            # soil_moist_chg => model_runoff%soil_moist_chg, &
            # soil_rechr_chg => model_runoff%soil_rechr_chg, &
            "sroff",
            "potet",
            # hru_ppt => model_precip%hru_ppt, & # JLM ??
            "transp_on",
            "snow_evap",
            "snowcov_area",
        )

    @staticmethod
    def get_variables() -> tuple:
        """Get soil zone variables

        Returns:
            variables: variables

        """
        # JLM: Seems like some of these should just be internal/postprocess
        # may not need to set all on self.
        return (
            "cap_infil_tot",
            "cap_waterin",
            "dunnian_flow",
            "hru_actet",
            "lakein_sz",
            "perv_actet",
            "potet_lower",
            "potet_rechr",
            "pref_flow",
            "pref_flow_in",
            "pref_flow_infil",
            "pref_flow_max",
            "pref_flow_stor",
            "pref_flow_thrsh",
            "recharge",
            "slow_flow",
            "slow_stor",
            "snow_free",
            "soil_lower",
            "soil_lower_ratio",
            "soil_moist",
            "soil_moist_tot",
            "soil_rechr",
            "soil_to_gw",
            "soil_to_ssr",
            "soil_zone_max",
            "ssr_to_gw",
            "ssres_flow",
            "ssres_in",
            "ssres_stor",
            "swale_actet",
            "unused_potet",
            "upslope_dunnianflow",
            "upslope_interflow",
        )

    @staticmethod
    def get_restart_variables() -> tuple:
        """Get soil zone restart variables

        Returns:
            variables: restart variables
        """

        return ()

    @staticmethod
    def get_init_values() -> dict:
        """Get soil zone inital values

        Returns:
            dict: inital values for named variables
        """

        return {
            "cap_infil_tot": zero,
            "cap_waterin": zero,
            "dunnian_flow": zero,
            "hru_actet": zero,
            "perv_actet": zero,
            "potet_lower": zero,
            "potet_rechr": zero,
            "pref_flow": zero,
            "pref_flow_in": zero,
            "pref_flow_infil": zero,
            "pref_flow_max": zero,
            "pref_flow_stor": zero,
            "pref_flow_thrsh": zero,
            "recharge": zero,
            "slow_flow": zero,
            "slow_stor": zero,
            "snow_free": nan,  # sm_soilzone
            "soil_lower": nan,  # completely set later
            "soil_lower_ratio": zero,
            "soil_lower_stor_max": nan,  # completely set later
            "soil_moist": nan,  # sm_climateflow
            "soil_moist_tot": nan,  # completely set later
            "soil_rechr": nan,  # sm_climateflow
            "soil_to_gw": zero,
            "soil_to_ssr": zero,
            "soil_zone_max": nan,  # this is completely later
            "ssr_to_gw": zero,
            "ssres_flow": zero,
            "ssres_in": zero,
            "ssres_stor": nan,  # sm_soilzone
            "swale_actet": zero,
            "unused_potet": zero,
        }

    def set_initial_conditions(self):
        """Initialize PRMSSoilzone variables."""

        # Derived parameters
        # JLM: is this awkward here?
        # JLM: it's definitely awkward to edit a parameter. maybe
        #      changes/checks on params should throw errors so
        #      parameter values remain the responsibility of the users
        #      and their deficiencies are transparent?
        self.hru_area_imperv = self.hru_percent_imperv * self.hru_area
        self.hru_area_perv = self.hru_area - self.hru_area_imperv
        self.hru_frac_perv = one - self.hru_percent_imperv
        self.soil_rechr_max = self.soil_rechr_max_frac * self.soil_moist_max
        self.snow_free = one - self.snowcov_area
        # edit a param
        wh_inactive_or_lake = np.where(
            (self.hru_type == HruType.INACTIVE)
            | (self.hru_type == HruType.LAKE)
        )
        self.sat_threshold[wh_inactive_or_lake] = 0.0
        # edit a param
        wh_not_land = np.where(self.hru_type != HruType.LAND)
        self.pref_flow_den[wh_not_land] = zero

        # variables
        if self.control.config["init_vars_from_file"] in [0, 2, 5]:
            # these are set in sm_climateflow
            self.soil_moist = self.soil_moist_init_frac * self.soil_moist_max
            self.soil_rechr = self.soil_rechr_init_frac * self.soil_rechr_max
        else:
            # call ctl_data%read_restart_variable(
            #    'soil_moist', this%soil_moist)
            # call ctl_data%read_restart_variable(
            #    'soil_rechr', this%soil_rechr)
            raise RuntimeError("Soilzone restart capability not implemented")

        # Note that the following is often editing a copy of the parameters,
        # which is a bit odd in that it might be contrary to the users
        # expectations. Move this parameter business to __init__

        # Do these need set on self or added to variables?
        self.grav_dunnian_flow = np.full(self.nhru, zero, "float64")
        self.pfr_dunnian_flow = np.full(self.nhru, zero, "float64")

        # ssres_stor
        if self.control.config["init_vars_from_file"] in [0, 2, 5]:
            self.ssres_stor = self.ssstor_init_frac * self.sat_threshold
            wh_inactive_or_lake = np.where(
                (self.hru_type == HruType.INACTIVE)
                | (self.hru_type == HruType.LAKE)
            )
            self.ssres_stor[wh_inactive_or_lake] = 0.0
            del self.ssstor_init_frac

        else:
            # call ctl_data%read_restart_variable(
            #     'ssres_stor', this%ssres_stor)
            raise RuntimeError("Soilzone restart capability not implemented")

        # <
        # need to set on swale_limit self? move to variables?
        self.swale_limit = np.full(self.nhru, zero, "float64")
        wh_swale = np.where(self.hru_type == HruType.SWALE)
        self.swale_limit[wh_swale] = 3.0 * self.sat_threshold[wh_swale]

        self.pref_flow_thrsh[wh_swale] = self.sat_threshold[wh_swale]
        wh_land = np.where(self.hru_type == HruType.LAND)
        self.pref_flow_thrsh[wh_land] = self.sat_threshold[wh_land] * (
            one - self.pref_flow_den[wh_land]
        )
        self.pref_flow_max[wh_land] = (
            self.sat_threshold[wh_land] - self.pref_flow_thrsh[wh_land]
        )

        # Need to set pref_flow_flag on self? or add to variables?
        wh_land_and_prf_den = np.where(
            (self.hru_type == HruType.LAND) & (self.pref_flow_den > zero)
        )
        self.pref_flow_flag = np.full(self.nhru, False, dtype=int)
        self.pref_flow_flag[wh_land_and_prf_den] = True

        # can this one be combined with the restart read logic above?
        if self.control.config["init_vars_from_file"] in [0, 2, 5]:
            wh_land_or_swale = np.where(
                (self.hru_type == HruType.LAND)
                | (self.hru_type == HruType.SWALE)
            )
            # JLM: to verify this
            self.slow_stor[wh_land_or_swale] = np.minimum(
                self.ssres_stor[wh_land_or_swale],
                self.pref_flow_thrsh[wh_land_or_swale],
            )
            self.pref_flow_stor[wh_land_or_swale] = (
                self.ssres_stor[wh_land_or_swale]
                - self.slow_stor[wh_land_or_swale]
            )

        else:
            # call ctl_data%read_restart_variable(
            #    'pref_flow_stor', self.pref_flow_stor)
            # call ctl_data%read_restart_variable('slow_stor', self.slow_stor)
            raise RuntimeError("Soilzone restart capability not implemented")

        # <
        # need to set soil2gw_flag on self? move to variables?
        self.soil2gw_flag = np.full(self.nhru, False, dtype=int)
        wh_soil2gwmax = np.where(self.soil2gw_max > 0.0)
        self.soil2gw_flag[wh_soil2gwmax] = True

        self.soil_zone_max = (
            self.sat_threshold + self.soil_moist_max * self.hru_frac_perv
        )
        self.soil_moist_tot = (
            self.ssres_stor + self.soil_moist * self.hru_frac_perv
        )

        self.soil_lower = self.soil_moist - self.soil_rechr
        self.soil_lower_stor_max = self.soil_moist_max - self.soil_rechr_max

        wh_soil_lower_stor = np.where(self.soil_lower_stor_max > 0.0)
        self.soil_lower_ratio[wh_soil_lower_stor] = (
            self.soil_lower[wh_soil_lower_stor]
            / self.soil_lower_stor_max[wh_soil_lower_stor]
        )

        # JLM: Things allocated that I dont immediately see
        # allocate(self.gvr2pfr(nhru))
        # self.gvr2pfr = 0.0_dp

        return

    def _advance_variables(self) -> None:
        return

    def calculate(self, simulation_time):
        """Calculate soil zone for a time step"""
        pass
