import numpy as np

from pynhm.base.storageUnit import StorageUnit

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import nan


class PRMSGroundwater(StorageUnit):
    """PRMS groundwater reservoir.

    Args:
    """

    def __init__(
        self,
        control: Control,
        soil_to_gw: adaptable,
        ssr_to_gw: adaptable,
        dprst_seep_hru: adaptable,
        budget_type: str = None,
        verbose: bool = False,
    ) -> "PRMSGroundwater":

        super().__init__(
            control=control,
            verbose=verbose,
        )
        self.name = "PRMSGroundwater"

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
            "gwres_flow",  # todo: privatize this var and keep vol public
            "gwres_flow_vol",
            "gwres_sink",
            "gwres_stor",
            "gwres_stor_old",
            "gwres_stor_change",
        )

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "soil_to_gw",
                "ssr_to_gw",
                "dprst_seep_hru",
            ],
            "outputs": [
                "gwres_flow",
            ],
            "storage_changes": [
                "gwres_stor_change",
            ],
        }

    @staticmethod
    def get_init_values() -> dict:
        """Get groundwater initial values

        Returns:
            dict: initial values for named variables
        """
        return {
            "gwres_flow": nan,
            "gwres_flow_vol": nan,
            "gwres_sink": nan,
            "gwres_stor": nan,
            "gwres_stor_old": nan,
            "gwres_stor_change": nan,
        }

    def _set_initial_conditions(self):
        # initialize groundwater reservoir storage
        self.gwres_stor[:] = self.gwstor_init.copy()
        self.gwres_stor_old[:] = self.gwstor_init.copy()
        return

    def _advance_variables(self) -> None:
        """Advance the groundwater reservoir variables
        Returns:
            None
        """
        self.gwres_stor_old[:] = self.gwres_stor
        return

    def _calculate(self, simulation_time):
        """Calculate groundwater reservoir terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

        gwarea = self.hru_area

        # calculate volume terms
        # denote volume units with leading underscore
        soil_to_gw_vol = self.soil_to_gw * gwarea
        ssr_to_gw_vol = self.ssr_to_gw * gwarea
        dprst_seep_hru_vol = self.dprst_seep_hru * gwarea

        # todo: what about route order

        _gwres_stor = self.gwres_stor * gwarea
        _gwres_stor += soil_to_gw_vol + ssr_to_gw_vol + dprst_seep_hru_vol

        _gwres_flow = _gwres_stor * self.gwflow_coef
        _gwres_stor -= _gwres_flow

        _gwres_sink = _gwres_stor * self.gwsink_coef
        idx = np.where(_gwres_sink > _gwres_stor)
        _gwres_sink[idx] = _gwres_stor[idx]
        _gwres_stor -= _gwres_sink

        # convert most units back to self variables
        # output variables
        self.gwres_stor[:] = _gwres_stor / gwarea
        # for some stupid reason this is left in acre-inches
        self.gwres_flow[:] = _gwres_flow / gwarea
        self.gwres_sink[:] = _gwres_sink / gwarea

        self.gwres_stor_change[:] = self.gwres_stor - self.gwres_stor_old
        self.gwres_flow_vol[:] = (
            self.gwres_flow * self.control.params.hru_in_to_cf
        )

        return
