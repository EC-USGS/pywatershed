from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control
from ..constants import HruType, zero

adaptable = Union[str, np.ndarray, Adapter]

NEARZERO = 1.0e-6
DNEARZERO = np.finfo(float).eps  # EPSILON(0.0D0)

RAIN = 0
SNOW = 1

BARESOIL = 0
GRASSES = 1

OFF = 0
ACTIVE = 1

LAND = HruType.LAND
LAKE = HruType.LAKE


class PRMSRunoff(StorageUnit):
    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        soil_moist: adaptable,
        net_rain: adaptable,
        verbose: bool = False,
    ):

        self.name = "PRMSRunoff"
        verbose = True
        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
            subclass_name=self.name,
        )

        # store dependencies
        self._input_variables_dict = {}
        self._input_variables_dict["soil_moist"] = adapter_factory(
            soil_moist, "soil_moist"
        )
        self._input_variables_dict["net_rain"] = adapter_factory(
            net_rain, "net_rain"
        )
        return

    def set_initial_conditions(self):
        # Where does the initial storage come from? Document here.
        # apparently it's just zero?
        # self.var1_stor[:] = np.zeros([1])[0]
        # self.var1_stor_old = None
        return

    @staticmethod
    def get_parameters() -> tuple:
        """
        Return a list of the parameters required for this process

        """
        return (
            "nhru",
            "carea_max",
            "smidx_coef",
            "smidx_exp",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get canopy reservoir input variables

        Returns:
            variables: input variables

        """
        return ("soil_moist", "net_rain")

    @staticmethod
    def get_variables() -> tuple:
        """Get output variables

        Returns:
            variables: output variables

        """
        return ("contrib_fraction", "infil")

    @staticmethod
    def get_init_values() -> dict:
        """Get initial values

        Returns:
            dict: initial values for named variables
        """

        return {
            "contrib_fraction": zero,
            "infil": zero,
        }

    def advance(self):
        """Advance box
        Returns:
            None

        """
        # self.var1_stor_old = self.var1_stor
        self._itime_step += 1

        for key, value in self._input_variables_dict.items():
            value.advance()
            v = getattr(self, key)
            v[:] = value.current

        return

    def calculate(self, time_length, vectorized=False):
        """Calculate terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """
        if vectorized:
            self.calculate_vectorized(time_length)
        else:
            self.calculate_procedural(time_length)
        return

    def calculate_procedural(self, time_length):

        # Make calculations
        for i in range(self.nhru):

            self.infil[i] = 0.0
            srp = 0.0

            # If rain/snow event with no antecedent snowpack, compute the runoff from the
            # rain first and then proceed with the snowmelt computations.
            self.infil[i] = self.infil[i] + self.net_rain[i]

            precip = self.net_rain[i]
            self.contrib_fraction[i] = self.calculate_contrib_fraction(
                self.carea_max[i],
                self.smidx_coef[i],
                self.smidx_exp[i],
                self.soil_moist[i],
                precip,
            )

            self.infil[i], srp = self.compute_infil_srp(
                self.contrib_fraction[i], self.net_rain[i], self.infil[i], srp
            )

        return

    def calculate_vectorized(self, time_length):

        # Make calculations
        precip = self.net_rain
        self.contrib_fraction = self.calculate_contrib_fraction(
            self.carea_max,
            self.smidx_coef,
            self.smidx_exp,
            self.soil_moist,
            precip,
        )

        return

    @staticmethod
    def calculate_contrib_fraction(
        carea_max, smidx_coef, smidx_exp, soil_moist, precip
    ):
        smidx = soil_moist + (0.5 * precip)
        contrib_frac = smidx_coef * 10.0 ** (smidx_exp * smidx)
        contrib_frac = min(carea_max, contrib_frac)
        return contrib_frac

    @staticmethod
    def compute_infil_srp(contrib_frac, precip, infil, perv_runoff):
        adjusted_precip = contrib_frac * precip
        infil = infil - adjusted_precip
        perv_runoff = perv_runoff + adjusted_precip
        return infil, perv_runoff
