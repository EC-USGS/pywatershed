from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control

from ..constants import one, zero

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

        # Derived parameters
        self.hru_area_imperv = self.hru_percent_imperv * self.hru_area
        self.hru_area_perv = self.hru_area - self.hru_area_imperv
        self.hru_frac_perv = one - self.hru_percent_imperv
        self.soil_rechr_max = self.soil_rechr_max_frac * self.soil_moist_max

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
            # srunoff_updated_soil => model_runoff%srunoff_updated_soil, &
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
        # Seems like some of these should just be internal
        return (
            "cap_infil_tot",
            "cap_waterin",
            "dunnian_flow",
            "hru_actet",
            "hru_sz_cascadeflow",
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

        return {}

    def set_initial_conditions(self):
        """Initialize PRMSSoilzone variables."""

        if self.control.config["init_vars_from_file"] in [0, 2, 3]:

            pass

        else:

            raise RuntimeError("Soilzone restart capability not implemented")
            # JLM: a list of restart variables dosent shed light on what states
            # actually have memory. could there be a diagnostic part of the
            # advance?

        return

    def _advance_variables(self) -> None:
        # self.pkwater_ante = self.pkwater_equiv.copy()
        return

    def calculate(self, simulation_time):
        """Calculate soil zone for a time step"""
        pass
