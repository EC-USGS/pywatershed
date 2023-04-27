import numpy as np

from pywatershed.base.storageUnit import StorageUnit

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import nan


class SWBRootZone(StorageUnit):
    """SWB Root Zone Storage object.

    Args:

    """

    def __init__(
        self,
        control: Control,
        net_rain: adaptable,
        snowmelt: adaptable,
        potet: adaptable,
        budget_type: str = None,
        calc_method: str = None,
        verbose: bool = False,
        load_n_time_batches: int = 1,
    ) -> "SWBRootZone":
        super().__init__(
            control=control,
            verbose=verbose,
            load_n_time_batches=load_n_time_batches,
        )
        self.name = "SWBRootZone"

        self._calc_method = str(calc_method)

        self._set_inputs(locals())
        self._set_budget(budget_type)
        return

    @staticmethod
    def get_dimensions() -> tuple:
        """Get dimensions"""
        return ("nhru",)

    @staticmethod
    def get_parameters() -> tuple:
        """Get root zone reservoir parameters

        Returns:
            parameters: input parameters

        """
        return (
            "base_curve_number",
            "max_net_infiltration_rate",
            "available_water_capacity",
            "rooting_depth",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get root zone reservoir input variables

        Returns:
            variables: input variables

        """
        return ("net_rain", "snowmelt", "potet")

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "net_rain",
                "snowmelt",
            ],
            "outputs": [
                "swb_runoff",
                "net_infiltration",
            ],
            "storage_changes": [
                "storage_change",
            ],
        }

    @staticmethod
    def get_init_values() -> dict:
        """Get SWB rootzone initial values

        Returns:
            dict: initial values for named variables
        """
        return {
            "swb_runoff": nan,
            "storage": nan,
            "storage_old": nan,
            "storage_change": nan,
        }

    def _set_initial_conditions(self):
        # initialize SWB rootzone storage
        self.swb_runoff[:] = 0.0
        self.net_infiltration = 0.0
        # self.storage[:] = ?
        # self.storage_old[:] = ?
        return

    def _advance_variables(self) -> None:
        """Advance the SWBRunoffCurveNumber
        Returns:
            None
        """
        self.storage_old[:] = self.storage
        return

    def _calculate(self, simulation_time):
        """Calculate SWBRunoffCurveNumber for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

        (
            self.runoff[:],
            self.storage[:],
        ) = self._calculate_numpy(
            self.base_curve_number,
            self.net_rain,
            self.snowmelt,
            self.storage_old,
        )
        return

    @staticmethod
    def _calculate_numpy(
        base_curve_number,
        net_rain,
        snowmelt,
        storage_old,
    ):
        # TBD

        return (runoff, storage)
