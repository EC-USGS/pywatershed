from typing import Union

import numpy as np

from pynhm.atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control

adaptable = Union[str, np.ndarray, Adapter]


class PRMSGroundwater(StorageUnit):
    """PRMS groundwater reservoir

    Args:
        params: parameter object
        atm: atmosphere object

    """

    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        soil_to_gw: adaptable,
        ssr_to_gw: adaptable,
        dprst_seep_hru: adaptable,
        verbose: bool = False,
    ) -> "PRMSGroundwater":

        self.name = "PRMSGroundwater"
        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
            subclass_name=self.name,
        )

        self._input_variables_dict = {}
        self._input_variables_dict["soil_to_gw"] = adapter_factory(
            soil_to_gw, "soil_to_gw"
        )
        self._input_variables_dict["ssr_to_gw"] = adapter_factory(
            ssr_to_gw,
            "ssr_to_gw",
        )
        self._input_variables_dict["dprst_seep_hru"] = adapter_factory(
            dprst_seep_hru,
            "dprst_seep_hru",
        )

        return

    def set_initial_conditions(self):
        # initialize groundwater reservoir storage
        self.gwres_stor = self.gwstor_init.copy()
        self.gwres_stor_old = self.gwstor_init.copy()
        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get groundwater reservoir parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nhru",
            "ngw",
            "hru_area",
            "gwflow_coef",
            "gwsink_coef",
            "gwstor_init",
            "gwstor_min",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get groundwater reservoir input variables

        Returns:
            variables: input variables

        """
        return (
            "soil_to_gw",
            "ssr_to_gw",
            "dprst_seep_hru",
        )

    @staticmethod
    def get_variables() -> tuple:
        """Get groundwater reservoir output variables

        Returns:
            variables: output variables

        """
        return (
            "gwres_flow",
            "gwres_in",
            "gwres_sink",
            "gwres_stor",
        )

    @staticmethod
    def get_init_values() -> dict:
        """Get groundwater reservoir inital values

        Returns:
            dict: inital values for named variables
        """
        # No GW res values need initialized prior to calculation.
        return {}

    def advance(self) -> None:
        """Advance the groundwater reservoir
        Returns:
            None

        """
        self.gwres_stor_old = self.gwres_stor
        self._itime_step += 1

        for key, value in self._input_variables_dict.items():
            value.advance()
            v = getattr(self, key)
            v[:] = value.current

        return

    def calculate(self, simulation_time):
        """Calculate groundwater reservoir terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

        gwarea = self.hru_area

        # calculate volume terms
        # gwstor_min_vol = self.gwstor_min * gwarea
        gwres_stor = self.gwres_stor * gwarea
        soil_to_gw_vol = self.soil_to_gw * gwarea
        ssr_to_gw_vol = self.ssr_to_gw * gwarea
        dprst_seep_hru_vol = self.dprst_seep_hru * gwarea

        # initialize calculation variables
        gwres_in = soil_to_gw_vol + ssr_to_gw_vol + dprst_seep_hru_vol

        # todo: what about route order

        gwres_stor += gwres_in

        gwres_flow = gwres_stor * self.gwflow_coef

        gwres_stor -= gwres_flow

        gwres_sink = gwres_stor * self.gwsink_coef
        idx = np.where(gwres_sink > gwres_stor)
        gwres_sink[idx] = gwres_stor[idx]

        gwres_stor -= gwres_sink

        # output variables
        self.gwres_stor = gwres_stor / gwarea
        self.gwres_in = (
            gwres_in  # for some stupid reason this is left in acre-inches
        )
        self.gwres_flow = gwres_flow / gwarea
        self.gwres_sink = gwres_sink / gwarea

        return
