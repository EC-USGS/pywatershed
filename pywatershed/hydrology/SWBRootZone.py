from pywatershed.base.storageUnit import StorageUnit

from pywatershed.base.adapter import adaptable
from pywatershed.base.control import Control
from pywatershed.constants import nan, zero
import pywatershed.functions.runoff_curve_number as cn


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
            "swb_rootzone_storage_init",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get root zone reservoir input variables
        Returns:
            variables: input variables
        """
        return (
            "net_rain",
            "snowmelt",
            "potet",
        )

    @staticmethod
    def get_init_values() -> dict:
        """Get SWB rootzone initial values
        Returns:
            dict: initial values for named variables
        """
        return {
            "swb_runoff": zero,
            "swb_rootzone_storage": zero,
            "swb_rootzone_storage_old": zero,
            "swb_rootzone_storage_change": zero,
            "swb_net_infiltration": zero,
        }

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "net_rain",
                "snowmelt",
            ],
            "outputs": [
                "swb_runoff",
                "swb_net_infiltration",
            ],
            "storage_changes": [
                "swb_rootzone_storage_change",
            ],
        }

    def _set_initial_conditions(self):
        self.swb_rootzone_storage[:] = self.swb_rootzone_storage_init.copy()
        return

    def _advance_variables(self) -> None:
        """Advance the SWBRunoffCurveNumber
        Returns:
            None
        """
        self.swb_rootzone_storage_old[:] = self.swb_rootzone_storage
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
            self.swb_runoff[:],
            self.swb_rootzone_storage[:],
        ) = self._calculate_numpy(
            base_curve_number=self.base_curve_number,
            net_rain=self.net_rain,
            snowmelt=self.snowmelt,
            storage_old=self.swb_rootzone_storage_old,
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
        storage = storage_old + net_rain
        runoff = snowmelt

        return (runoff, storage)
