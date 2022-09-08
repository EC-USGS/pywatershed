from typing import Union
from warnings import warn

import numpy as np

from pynhm.base.control import Control

from ..constants import zero
from ..utils.formatting import pretty_print
from .accessor import Accessor

# Todo
# * documentation
# * units (check consistent): know about metadata and just look up names?
# * get time_units / time_step from control
# * decide on the fate of "from storage unit"


class Budget(Accessor):
    """Budget class for mass and energy conservation.

    Currently no energy budget has been implmenented, todo.
    """

    def __init__(
        self,
        control: Control,
        inputs: Union[list, dict],
        outputs: Union[list, dict],
        storage_changes: Union[list, dict],
        init_accumulations: dict = None,
        accum_start_time: np.datetime64 = None,
        units: str = None,
        time_unit: str = "D",  # TODO: get this from control
        description: str = None,
        rtol: float = 1e-5,
        atol: float = 1e-5,
        basis: str = "unit",
        imbalance_fatal: bool = False,
        verbose: bool = True,
    ):

        self.name = "Budget"
        self.control = control
        self.inputs = self.init_component(inputs)
        self.outputs = self.init_component(outputs)
        self.storage_changes = self.init_component(storage_changes)
        self.time_unit = time_unit
        self.description = description
        self.rtol = rtol
        self.atol = atol
        self.imbalance_fatal = imbalance_fatal
        self.verbose = verbose
        self.basis = basis

        self._inputs_sum = None
        self._outputs_sum = None
        self._storage_changes_sum = None
        self._balance = None
        self._accumulations = None
        self._accumulations_sum = None
        self._zero_sum = None

        self._time = self.control.current_time
        self._itime_step = self.control.itime_step

        # metadata
        all_vars = [list(self[cc].keys()) for cc in self.get_components()]
        all_vars = [x for xs in all_vars for x in xs]
        self.meta = self.control.meta.get_vars(all_vars)
        if len(self.meta) == len(all_vars):
            all_units = [val["units"] for val in self.meta.values()]
            # check consistent units
            if not (np.array(all_units) == all_units[0]).all():
                msg = "Units not consistent over all terms"
                raise ValueError(msg)
            self.units = all_units[0]
        else:
            msg = f"Metadata unavailable for some Budget terms in {all_vars}"
            warn(msg)
            self.units = None

        self.set_initial_accumulations(init_accumulations, accum_start_time)
        return

    @staticmethod
    def init_component(data: Union[list, dict]) -> dict:
        if isinstance(data, dict):
            return data
        else:
            return {dd: None for dd in data}

    def set(self, data: dict):
        """Set the data on the components after initialization
        Args:
            data: a dict of dicts with top level optional keys:
            [inputs, outputs, storage_changes].
            Each of those is a dict with var: np.ndarray, eg.
            data = {'inputs': {'var': np.ndarray}}
        """
        for comp_name, comp_dict in data.items():
            for var_name, var_data in comp_dict.items():
                if self[comp_name][var_name] is not None:
                    msg = (
                        f"Component '{comp_name}' variable '{var_name}'"
                        f"has already been set and should not be reset."
                    )
                    raise ValueError(msg)
                elif var_name not in self[comp_name].keys():
                    msg = (
                        f"Component '{comp_name}' has no variable '{var_name}'"
                    )
                    raise KeyError(msg)
                else:
                    self[comp_name][var_name] = var_data

    @classmethod
    def from_storage_unit(cls, storage_unit, **kwargs):
        mass_budget_terms = storage_unit.get_mass_budget_terms()
        for component in mass_budget_terms.keys():
            kwargs[component] = {}
            for var in mass_budget_terms[component]:
                kwargs[component][var] = storage_unit[var]

        return Budget(storage_unit.control, **kwargs)

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

    def set_initial_accumulations(self, init_accumulations, accum_start_time):
        self._itime_accumulated = self._itime_step  # -1
        self._time_accumulated = self._time  # None
        self._accumulations = {}
        self._accumulations_sum = {}

        # init to zero
        self.reset_accumulations()

        if init_accumulations is None:
            self._accum_start_time = self.control.init_time
            return

        self._accum_start_time = accum_start_time
        for component in self.components:
            if component not in init_accumulations.keys():
                continue
            for var in self[component].keys():
                if var in init_accumulations[component].keys():
                    self._accumulations[component][var] = init_accumulations[
                        component
                    ][var]

        self._sum_component_accumulations()
        return

    def reset_accumulations(self):
        self._accum_start_time = self._time_accumulated
        for component in self.components:
            self._accumulations[component] = {}
            for var in self[component].keys():
                self._accumulations[component][var] = zero
        self._sum_component_accumulations()
        return

    def advance(self):
        """Advance time (taken from storageUnit)"""
        if self._itime_step >= self.control.itime_step:
            if self.verbose:
                msg = (
                    f"{self.name} did not advance because it is "
                    f"not behind control time"
                )
                print(msg)
            return

        self._itime_step = self.control.itime_step
        self._time = self.control.current_time

    def calculate(self):
        """Accumulate for the timestep."""
        if self._itime_accumulated >= self._itime_step:
            raise ValueError("Can not accumulate twice per timestep")

        self._inputs_sum = self._sum_inputs()
        self._outputs_sum = self._sum_outputs()
        self._storage_changes_sum = self._sum_storage_changes()

        # accumulate
        for component in self.components:
            for var in self[component].keys():
                self._accumulations[component][var] += self[component][
                    var
                ] * self.control.time_step.astype(
                    f"timedelta64[{self.time_unit}]"
                ).astype(
                    int
                )

        self._sum_component_accumulations()

        # check balance
        if self.basis == "unit":
            self._balance = self._calc_unit_balance()
        elif self.basis == "global":
            self._balance = self._calc_global_balance()

        self._itime_accumulated = self._itime_step
        self._time_accumulated = self._time
        return

    def _sum_component_accumulations(self):
        # sum the individual component accumulations
        for component in self.components:
            self._accumulations_sum[component] = None
            for var in self[component].keys():
                if self._accumulations_sum[component] is None:
                    self._accumulations_sum[component] = self._accumulations[
                        component
                    ][var].copy()
                else:
                    self._accumulations_sum[component] += self._accumulations[
                        component
                    ][var]
        return

    @property
    def accumulations(self):
        return self._accumulations

    def _sum(self, attr):
        """Sum over the individual terms in a budget component."""
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

    def _calc_unit_balance(self):
        unit_balance = self._inputs_sum - self._outputs_sum
        self._zero_sum = True
        if not np.allclose(
            unit_balance,
            self._storage_changes_sum,
            rtol=self.rtol,
            atol=self.atol,
        ):
            self._zero_sum = False
            msg = (
                "The flux unit balance not equal to the change in unit storage"
            )
            if self.imbalance_fatal:
                raise ValueError(msg)
            else:
                warn(msg, UserWarning)

        return unit_balance

    def _calc_global_balance(self):
        global_balance = self._inputs_sum.sum() - self._outputs_sum.sum()
        self._zero_sum = True
        if not np.allclose(
            global_balance,
            self._storage_changes_sum.sum(),
            rtol=self.rtol,
            atol=self.atol,
        ):
            self._zero_sum = False
            msg = (
                "The global flux balance not equal to the change in global "
                "storage."
            )
            if self.imbalance_fatal:
                raise ValueError(msg)
            else:
                warn(msg, UserWarning)

        return global_balance

    @property
    def balance(self):
        return self._balance

    def __repr__(self):

        if self._itime_step == -1:
            msg = (
                f"Budget (units: {self.units}) of {self.description} "
                f"only initialized. Check back later."
            )
            return msg

        n_in = len(self.inputs)
        n_out = len(self.outputs)
        n_stor = len(self.storage_changes)
        n_report = max(n_in, n_out, n_stor)
        in_keys = list(self.inputs.keys())
        out_keys = list(self.outputs.keys())
        stor_keys = list(self.storage_changes.keys())

        # the following do not copy the ndarrays
        in_vals = list(self.inputs.values())
        out_vals = list(self.outputs.values())
        stor_vals = list(self.storage_changes.values())
        acc_in_vals = list(self._accumulations["inputs"].values())
        acc_out_vals = list(self._accumulations["outputs"].values())
        acc_stor_vals = list(self._accumulations["storage_changes"].values())

        #           ': ' + value spaces
        col_extra_colon = 2
        col_extra_vals = 10
        col_extra = col_extra_colon + col_extra_vals
        in_col_key_width = max([len(kk) for kk in in_keys])
        out_col_key_width = max([len(kk) for kk in out_keys])
        stor_col_key_width = max([len(kk) for kk in stor_keys])
        in_col_width = in_col_key_width + col_extra
        out_col_width = out_col_key_width + col_extra
        stor_col_width = stor_col_key_width + col_extra

        indent_width = 9
        indent_fill = " " * indent_width
        in_col_fill = " " * in_col_width
        out_col_fill = " " * out_col_width
        stor_col_fill = " " * stor_col_width
        sep_width = 4
        col_sep = " " * sep_width
        terms_width = (
            in_col_width + out_col_width + stor_col_width + (2 * sep_width)
        )
        total_width = indent_width + terms_width

        # volume/mass or energy
        # budget name
        summary = []
        summary += ["*-" * int((total_width) / 2)]

        if self.basis == "unit":
            summary += [
                f"Individual spatial unit budget for {self.description} "
                f"(units: {self.units}).\n"
                "Budget is checked on each spatial unit. This is summary shows"
                "\nspatial sums for the entire model domain.\n"
            ]
        elif self.basis == "global":
            summary += [
                f"Global budget of {self.description} (units: {self.units}).\n"
                "Budget is checked on full domain: spatially summed fluxes\n"
                "and storage changes are checked for balance."
            ]

        # print the model time. this is
        summary += [f"@ time: {self._time} (itime_step: {self._itime_step})"]

        # Timestep rates
        summary += [""]
        # header
        summary += ["This timestep, rates:"]
        summary += [
            indent_fill
            + "input rates".ljust(in_col_width)
            + col_sep
            + "output rates".ljust(out_col_width)
            + col_sep
            + "storage change rates".ljust(stor_col_width)
        ]
        separator = [
            (indent_fill + "-" * in_col_width + col_sep)
            + ("-" * out_col_width + col_sep)
            + ("-" * stor_col_width)
        ]
        summary += separator

        # terms line by line with cols: in, out, storage
        term_data = (
            (n_in, in_keys, in_vals, in_col_width, in_col_fill),
            (n_out, out_keys, out_vals, out_col_width, out_col_fill),
            (n_stor, stor_keys, stor_vals, stor_col_width, stor_col_fill),
        )

        def table_terms_col_wise():
            "Fill terms table column wise"
            table = []
            with np.printoptions(precision=2):
                for ll in range(n_report):
                    line = indent_fill
                    for n_items, keys, vals, col_width, col_fill in term_data:
                        if ll <= (n_items - 1):
                            nk = len(keys[ll])
                            # subtract ": ", "m." and "e+ee"
                            prec_width = col_width - nk - 2 - 2 - 4
                            vals_sum = vals[ll].sum()
                            if vals_sum < 0:
                                prec_width -= 1
                            vv = np.format_float_scientific(
                                vals_sum,
                                precision=prec_width,
                            )
                            line += (
                                pretty_print(f"{keys[ll]}: {vv}", col_width)
                                + col_sep
                            )
                        else:
                            line += col_fill + col_sep

                    line += col_sep
                    table += [line]

            return table

        summary += table_terms_col_wise()

        # balance line
        summary += [indent_fill + "-" * terms_width]

        eq_op = "="
        if not self._zero_sum:
            eq_op = "!=!"

        term_data = (
            ("", self._inputs_sum, in_col_width, in_col_key_width),
            ("-", self._outputs_sum, out_col_width, out_col_key_width),
            (
                eq_op,
                self._storage_changes_sum,
                stor_col_width,
                stor_col_key_width,
            ),
        )

        def balance_line_col_wise():
            bal_line = "Balance: "
            for oper, vals_sum, col_width, col_key_width in term_data:

                bal_line += (
                    oper
                    + (" " * (col_key_width + col_extra_colon - len(oper)))
                    + pretty_print(
                        np.format_float_scientific(
                            vals_sum.sum(),
                            precision=col_width
                            - col_key_width
                            - col_extra_colon
                            - 6,
                        ),
                        col_width - col_key_width - col_extra_colon,
                    )
                    + col_sep
                )
            return bal_line

        summary += [balance_line_col_wise()]

        # Accumulated volumes
        summary += [""]
        # header
        summary += [f"Accumulated volumes (since {self._accum_start_time}):"]
        summary += [
            indent_fill
            + "input volumes".ljust(in_col_width)
            + col_sep
            + "output volumes".ljust(out_col_width)
            + col_sep
            + "storage change volumes".ljust(out_col_width)
        ]
        # separator
        summary += separator

        # terms line by line with cols: in, out, storage
        term_data = (
            (n_in, in_keys, acc_in_vals, in_col_width, in_col_fill),
            (n_out, out_keys, acc_out_vals, out_col_width, out_col_fill),
            (n_stor, stor_keys, acc_stor_vals, stor_col_width, stor_col_fill),
        )
        # fill terms column wise (abuse scope a bit)
        summary += table_terms_col_wise()

        # balance line
        summary += [indent_fill + "-" * terms_width]
        term_data = (
            (
                "",
                self._accumulations_sum["inputs"],
                in_col_width,
                in_col_key_width,
            ),
            (
                "-",
                self._accumulations_sum["outputs"],
                out_col_width,
                out_col_key_width,
            ),
            (
                eq_op,
                self._accumulations_sum["storage_changes"],
                stor_col_width,
                stor_col_key_width,
            ),
        )
        summary += [balance_line_col_wise()]

        # conclude
        summary += [""]
        return "\n".join(summary)
