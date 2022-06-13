from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.budget import Budget
from ..base.control import Control
from ..constants import nan, one, zero

adaptable = Union[str, np.ndarray, Adapter]


class PRMSEt(StorageUnit):
    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        potet: adaptable,
        hru_impervevap: adaptable,
        hru_intcpevap: adaptable,
        snow_evap: adaptable,
        dprst_evap_hru: adaptable,
        perv_actet: adaptable,
        verbose: bool = False,
        imbalance_fatal: bool = False,
    ) -> "PRMSEt":

        self.name = "PRMSEt"
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

        # budget
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
        self.budget = Budget(**budget_terms, imbalance_fatal=imbalance_fatal)

        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get ET parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nhru",
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

    # @staticmethod
    # def get_et_terms() -> tuple:
    #     """Get the indiviudal ET terms

    #     Returns:
    #         variables: variables

    #     """
    #     return (
    #         "hru_impervevap",
    #         "hru_intcpevap",
    #         "snow_evap",
    #         "dprst_evap_hru",
    #         "hru_perv_actet",
    #     )

    def set_initial_conditions(self):
        """Initialize PRMSEt variables."""

        # Derived parameter
        # "hru_percent_imperv" is actually a fraction
        self.frac_perv = one - self["hru_percent_imperv"] - self["dprst_frac"]

    def _advance_variables(self) -> None:
        # for the available_potet we probably need to check that the
        # loss components are zero at the beginning of the timestep.
        # for vv in self.get_et_terms():
        #    assert (self[vv] == zero).all() or np.isnan(self[vv]).all()
        return

    def calculate(self, simulation_time):
        """Calculate actual ET for a time step"""
        self.hru_actet[:] = self["potet"] - self.available_potet
        self.budget.calculate()
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
