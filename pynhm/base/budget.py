from warnings import warn

import numpy as np

from ..constants import zero
from .accessor import Accessor

# # For both 1) cumulative (start_date), and 2) current time step
# # * Individual input terms: name, units, & subclass
# # * Individual output terms: name, units, & subclass
# # * Sum input terms
# # * Sum output terms
# # In - out
# # ? check against a reference storage?

# budget needs to be able to restart?
# should storage unit edit its own metadata to set it's subclass on each variable?

# the input dictionarys are a strange stucture
# key: {data: np.ndarray([nhru]), **metadata}
# might consider other ways of passing the data.
# ** should carry the data and the meta separately ** since that's what we do elsewhere

# apply initial accumulations passed via e.g.
# initial_accumulations = {
#     "inputs": {"var": data, ...},
#     "outputs": {"var": data, ...},
#     "storage_changes": {"var": data, ...}, }


class Budget(Accessor):
    def __init__(
        self,
        inputs: dict,
        outputs: dict,
        storage_changes: dict,
        meta: dict = None,
        init_accumulations: dict = None,
        verbosity: int = 0,
    ):
        self.name = "Budget"
        self.inputs = inputs
        self.outputs = outputs
        self.storage_changes = storage_changes
        self.meta = meta
        self.verbosity = verbosity

        self._inputs_sum = None
        self._outputs_sum = None
        self._storage_changes_sum = None
        self._balance = None
        self._accumulations = None

        self.set_initial_accumulations(init_accumulations)
        return

    @staticmethod
    def from_storage_unit(storage_unit, **kwargs):
        # assemble the meta data, which will determine the budget component variables
        kwargs["meta"] = {}
        kwargs["meta"]["inputs"] = storage_unit.meta._get_attr_key_val(
            "input", "var_category", "mass flux"
        )
        kwargs["meta"]["outputs"] = storage_unit.meta._get_attr_key_val(
            "var", "var_category", "mass flux"
        )
        kwargs["meta"][
            "storage_changes"
        ] = storage_unit.meta._get_attr_key_val(
            "var", "var_category", "mass storage change"
        )

        # keep only the desried metadata for budget
        meta_keys = Budget.get_meta_keys()
        for component in self.components:
            for key, val in kwargs["meta"][component].items():
                kwargs["meta"][component][key] = {
                    kk: vv for kk, vv in val.items() if kk in meta_keys
                }

        for component in self.components:
            kwargs[component] = {}
            for var in kwargs["meta"][component].keys():
                kwargs[component][var] = storage_unit[var][:]

        budget = Budget(**kwargs)
        budget.storage_unit_name = storage_unit.name
        return budget

    @staticmethod
    def get_meta_keys():
        """Return a tuple of the metadata keys used by Budget."""
        return ("desc", "modules", "var_category", "units")

    @staticmethod
    def get_components():
        return ("inputs", "outputs", "storage_changes")

    @property
    def components(self):
        return self.get_components()

    def set_initial_accumulations(self, init_accumulations):
        # init to zero
        self._accumulations = {}
        for component in self.components:
            self._accumulations[component] = {}
            for var in self[component].keys():
                self._accumulations[component][var] = (
                    self[component][var] * zero
                )

        if init_accumulations is None:
            return None

        for component in self.components:
            if component not in init_accumulations.keys():
                continue
            for var in self[component].keys():
                if var in init_accumulations[component].keys():
                    self._accumulations[component][var] = init_accumulations[
                        "component"
                    ][var]

        return None

    def accumulate(self):
        for component in self.components:
            for var in self[component].keys():
                self._accumulations[component][var] += self[component][var]
        return None

    @property
    def accumulations(self):
        return self._accumulations

    def _sum(self, attr):
        key0 = list(self[attr].keys())[0]
        sum = self[attr][key0] * zero
        for kk, vv in self[attr].items():
            sum += vv
        return sum

    def _sum_inputs(self):
        return self._sum("inputs")

    def _sum_outputs(self):
        return self._sum("outputs")

    def _sum_storage_changes(self):
        return self._sum("storage_changes")

    def _calc_balance(self):
        balance = self._inputs_sum - self._outputs_sum
        close = np.isclose(balance, self._storage_changes_sum)
        # diff = balance, self._storage_changes_sum)
        if not close.all():
            msg = "The flux balance not equal to the change in storage"
            raise ValueError(msg)
        return balance

    @property
    def balance(self):
        return self._balance

    def calculate(self):
        self._inputs_sum = self._sum_inputs()
        self._outputs_sum = self._sum_outputs()
        self._storage_changes_sum = self._sum_storage_changes()
        self._balance = self._calc_balance()
        self.accumulate()
        return None
