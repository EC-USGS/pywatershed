import numpy as np

from pywatershed.base.budget import Budget
from pywatershed.base.process import Process

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import nan, one, zero
from ..parameters import Parameters

# THis class seems to be a sort of intermediate class that should
# not exist in the future.
# PRMS: solution is calculate hru_actet at the bottom of the chain, in
#      soilzone.
# pywatershed solution: potet and hru_actet are in atmosphere. potet is an
#     "input" and avail_potet is passed around, resulting in hru_actet at the
#      end of each time calculation.


class PRMSEt(Process):
    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        potet: adaptable,
        hru_impervevap: adaptable,
        hru_intcpevap: adaptable,
        snow_evap: adaptable,
        dprst_evap_hru: adaptable,
        perv_actet: adaptable,
        verbose: bool = False,
        budget_type: str = None,
    ) -> "PRMSEt":
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "PRMSEt"

        self._set_inputs(locals())
        self._set_options(locals())

        # Cant set the default budget for ET
        # self.set_budget(budget_type)
        # because the input/output conventions dont match.
        # potet is an "input" but all the other inputs to
        # the class are actually outputs

        # explicitly declare the budget
        if budget_type is None:
            self.budget = None
        else:
            budget_terms = {
                "inputs": {"potet": self.potet},
                "outputs": {
                    "hru_impervevap": self.hru_impervevap,
                    "hru_intcpevap": self.hru_intcpevap,
                    "snow_evap": self.snow_evap,
                    "dprst_evap_hru": self.dprst_evap_hru,
                    "hru_perv_actet": self.hru_perv_actet,
                },
                "storage_changes": {"unused_potet": self.unused_potet},
            }
            self.budget = Budget(
                control=self.control,
                **budget_terms,
                description=self.name,
                imbalance_fatal=(budget_type == "strict"),
            )

        return

    @staticmethod
    def get_dimensions() -> tuple:
        """Get ET dimensions

        Returns:
            dimensions: input dimensions

        """
        return ("nhru",)

    @staticmethod
    def get_parameters() -> tuple:
        """Get ET parameters

        Returns:
            parameters: input parameters

        """
        return (
            "dprst_frac",
            "hru_percent_imperv",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get ET input variables

        Returns:
            variables: input variables

        """
        return (
            "potet",
            "hru_impervevap",
            "hru_intcpevap",
            "snow_evap",
            "dprst_evap_hru",
            "perv_actet",
        )

    @staticmethod
    def get_variables() -> tuple:
        """Get ET variables

        Returns:
            variables: variables

        """
        return ("unused_potet", "hru_perv_actet", "hru_actet")

    @staticmethod
    def get_init_values() -> dict:
        """Get ET inital values

        Returns:
            dict: inital values for named variables
        """

        return {
            "unused_potet": nan,
            "hru_perv_actet": nan,
            "hru_actet": nan,
        }

    def _set_initial_conditions(self):
        """Initialize PRMSEt variables."""

        # Derived parameter
        # "hru_percent_imperv" is actually a fraction
        self.frac_perv = one - self["hru_percent_imperv"] - self["dprst_frac"]

    def _advance_variables(self) -> None:
        # The "change in storage" is the hruactet, which is the amount of
        # water back to the atm. But that amount is not tracked over time,
        # it's not really a storage
        return

    def _calculate(self, simulation_time):
        """Calculate actual ET for a time step"""
        self.hru_actet[:] = self["potet"] - self["available_potet"]
        return

    @property
    def available_potet(self) -> np.ndarray:
        """remaining/unused/available potet property

        Experimental to see if it can be updated on request
        Needs zeros at the beginning of each timestep, see _advance_variables

        """
        self.hru_perv_actet[:] = self["perv_actet"] * self.frac_perv

        loss_vars = [
            "hru_intcpevap",
            "snow_evap",
            "dprst_evap_hru",
            "hru_impervevap",
            "hru_perv_actet",
        ]

        loss = None
        for ll in loss_vars:
            if loss is None:
                loss = self[ll].copy()
            else:
                loss += self[ll]

        self.unused_potet[:] = np.maximum(self["potet"] - loss, zero)

        return self.unused_potet
