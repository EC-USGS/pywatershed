from pywatershed.base.storageUnit import StorageUnit

from pywatershed.base.adapter import adaptable
from pywatershed.base.control import Control
from pywatershed.constants import nan, zero
import pywatershed.functions.runoff_curve_number as cn
import pywatershed.functions.actual_et_thornthwaite_mather as aet
import numpy as np


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
    ):
        super().__init__(
            control=control,
            verbose=verbose,
            load_n_time_batches=load_n_time_batches,
        )
        self.name = "SWBRootZone"

        self._calc_method = str(calc_method)

        self._set_inputs(locals())
        self._set_budget(budget_type)

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
            "swb_soil_storage_init",
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
            "swb_actual_et": zero,
            "swb_soil_storage": zero,
            "swb_soil_storage_old": zero,
            "swb_soil_storage_change": zero,
            "swb_net_infiltration": zero,
            "swb_soil_storage_max": zero,
            "swb_accumulated_potential_water_loss": zero,
            "swb_accumulated_potential_water_loss_old": zero,
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
                "swb_soil_storage_change",
            ],
        }

    def _set_initial_conditions(self):
        self.swb_soil_storage[:] = self.swb_soil_storage_init.copy()
        self.swb_soil_storage_max[:] = self.available_water_capacity * self.rooting_depth
        return

    def _advance_variables(self) -> None:
        """Advance the SWBRunoffCurveNumber
        Returns:
            None
        """
        self.swb_soil_storage_old[:] = self.swb_soil_storage
        return

    def _calculate(self, simulation_time):
        """Calculate SWBRunoffCurveNumber for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

#        self.accumulated_potential_water_loss_old[:] = self.swb_accumulated_potential_water_loss[:].copy()

        (
            self.swb_runoff[:],
            self.swb_soil_storage[:],
            self.swb_accumulated_potential_water_loss[:],
            self.swb_actual_et[:],
            self.swb_net_infiltration[:],
        ) = self._calculate_numpy(
            base_curve_number=self.base_curve_number,
            net_rain=self.net_rain,
            snowmelt=self.snowmelt,
            reference_et=self.potet,
            apwl=self.swb_accumulated_potential_water_loss,
            storage_max=self.swb_soil_storage_max,
            storage_old=self.swb_soil_storage_old,
        )
        return

    @staticmethod
    def _calculate_numpy(
        base_curve_number,
        net_rain,
        snowmelt,
        reference_et,
        apwl,
        storage_max,
        storage_old,
    ):
        apwl_new = np.full_like(apwl, fill_value=0.0)
        actual_et = np.full_like(apwl, fill_value=0.0)
        storage = np.full_like(storage_old, fill_value=0.0)
        net_infiltration = np.full_like(apwl, fill_value=0.0)

        inflow = net_rain + snowmelt
        curve_number_storage_S = cn.calculate_cn_S_inches(base_curve_number)
        runoff = cn.calculate_cn_runoff(inflow=inflow, storage_S=curve_number_storage_S)
        
        P_minus_PE = inflow - runoff - reference_et

        P_minus_PE_ge_0 = np.where(P_minus_PE >= 0.)
        P_minus_PE_lt_0 = np.where(P_minus_PE < 0.)

        # handle cases where P minus PET is >= 0
        if len(P_minus_PE_ge_0[0]) > 0:
            actual_et[P_minus_PE_ge_0] = reference_et[P_minus_PE_ge_0]
            temp_storage = storage_old[P_minus_PE_ge_0] + P_minus_PE[P_minus_PE_ge_0]
            zeros = np.full_like(temp_storage, fill_value=0.0)
            storage[P_minus_PE_ge_0] = np.min((storage_max[P_minus_PE_ge_0], temp_storage))
            net_infiltration[P_minus_PE_ge_0] = np.max((zeros, temp_storage - storage_max[P_minus_PE_ge_0]))

            apwl_new[P_minus_PE_ge_0] = aet.thornthwaite_mather_accumulated_potential_water_loss_inches(max_soil_moisture=storage_max[P_minus_PE_ge_0],
                                                                                                        soil_moisture=storage)        
            actual_et[P_minus_PE_ge_0] = reference_et[P_minus_PE_ge_0]
            
        # now handle cases where P minus PET < 0
        if len(P_minus_PE_lt_0[0]) > 0:
            apwl_new[P_minus_PE_lt_0] = apwl_new[P_minus_PE_lt_0] + np.abs(P_minus_PE[P_minus_PE_lt_0])        
            storage[P_minus_PE_lt_0] = aet.thornthwaite_mather_soil_moisture_inches(max_soil_moisture=storage_max[P_minus_PE_lt_0],
                                                                                    apwl=apwl_new[P_minus_PE_lt_0])
            actual_et[P_minus_PE_lt_0] = storage_old[P_minus_PE_lt_0] - storage[P_minus_PE_lt_0]        
                
        #breakpoint()
                
        return (runoff, storage, apwl_new, actual_et, net_infiltration)
