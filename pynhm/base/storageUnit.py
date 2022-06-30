import os
import pathlib as pl

import numpy as np

from pynhm.base.budget import Budget

from ..base.adapter import Adapter, adapter_factory
from ..utils.netcdf_utils import NetCdfWrite
from .accessor import Accessor
from .control import Control

# These could be changed in the metadata yaml
type_translation = {
    "D": "float64",
    "F": "float64",
    "I": "int32",
    "B": "bool",  # not used despite the popularity of "flags"
}


class StorageUnit(Accessor):
    def __init__(
        self,
        control: Control,
        verbose: bool,
        subclass_name="StorageUnit",
    ):

        self.name = subclass_name
        self.control = control
        self.params = self.control.params
        self.verbose = verbose
        self._simulation_time = 0.0

        # netcdf output variables
        self._output_netcdf = False
        self._netcdf = None
        self._separate_netcdf = True
        self._itime_step = -1

        self.get_metadata()
        self.initialize_self_variables()
        self.set_initial_conditions()

        return None

    def output(self) -> None:
        """Output

        Returns:
            None

        """
        if self._output_netcdf:
            self.__output_netcdf()
        return

    def finalize(self) -> None:
        """Finalize storageUnit

        Returns:
            None

        """
        self._finalize_netcdf()
        return

    @staticmethod
    def get_parameters() -> list:
        raise Exception("This must be overridden")

    @staticmethod
    def get_inputs() -> list:
        raise Exception("This must be overridden")

    @staticmethod
    def get_variables() -> list:
        raise Exception("This must be overridden")

    @staticmethod
    def get_restart_variables() -> list:
        raise Exception("This must be overridden")

    @staticmethod
    def get_init_values() -> list:
        raise Exception("This must be overridden")

    @property
    def parameters(self) -> list:
        return self.get_parameters()

    @property
    def inputs(self) -> list:
        return self.get_inputs()

    @property
    def variables(self) -> list:
        return self.get_variables()

    @property
    def restart_variables(self) -> list:
        return self.get_restart_variables()

    @property
    def init_values(self) -> list:
        return self.get_init_values()

    def initialize_self_variables(self, restart: bool = False):
        for name in self.parameters:
            setattr(self, name, self.params.parameters[name])
        for name in self.inputs:
            setattr(self, name, np.zeros(self.nhru, dtype=float) + np.nan)
        # skip restart variables if restart (for speed) ?
        # the code is below but commented.
        # restart_variables = self.restart_variables
        for name in self.variables:
            # if restart and (name in restart_variables):
            #     continue
            self.initialize_var(name)

        # demote any self variables to single?
        # for vv in self.__dict__.keys():
        # for vv in self.parameters:
        #     if (
        #         isinstance(self[vv], np.ndarray)
        #         and self[vv].dtype == "float64"
        #     ):
        #         print(f"converting {vv} from float64 to float32")
        #         self[vv] = self[vv].astype("float32")

        return

    def initialize_var(self, var_name):
        if self.var_meta[var_name]["dimensions"][0] == "nsegment":
            nsize = self.nsegment
        else:
            nsize = self.nhru
        init_vals = self.get_init_values()
        if var_name in init_vals.keys():
            init_type = type_translation[self.var_meta[var_name]["type"]]
            setattr(
                self,
                var_name,
                np.full(nsize, init_vals[var_name], dtype=init_type),
            )
        elif self.verbose:
            print(f"{var_name} not initialized (no initial value specified)")

        return

    def set_initial_conditions(self):
        raise Exception("This must be overridden")

    def _advance_variables(self):
        raise Exception("This must be overridden")

    def _advance_inputs(self):
        for key, value in self._input_variables_dict.items():
            value.advance()
            self[key][:] = value.current

        return

    def set_inputs(self, args):
        self._input_variables_dict = {}
        for ii in self.inputs:
            self._input_variables_dict[ii] = adapter_factory(
                args[ii], ii, args["control"]
            )
            if self._input_variables_dict[ii]:
                self[ii] = self._input_variables_dict[ii].current

        return

    def set_input_to_adapter(self, input_variable_name: str, adapter: Adapter):

        self._input_variables_dict[input_variable_name] = adapter
        # can NOT use [:] on the LHS as we are relying on pointers between
        # boxes. [:] on the LHS here means it's not a pointer and then
        # requires that the calculation of the input happens before the
        # advance of this storage unit. But that gives the incorrect budget
        # for et.
        self[input_variable_name] = adapter.current

        # Using a pointer between boxes means that the same pointer has to
        # be used for the budget, so there's no way to have a preestablished
        # pointer between storageUnit and its budget. So this stuff...
        if self.budget is not None:
            for comp in self.budget.components:
                if input_variable_name in self.budget[comp].keys():
                    # can not use [:] on the LHS?
                    self.budget[comp][input_variable_name] = self[
                        input_variable_name
                    ]

        return

    def set_budget(self, budget_type):
        if budget_type is None:
            self.budget = None
        elif budget_type in ["strict", "diagnostic"]:
            self.budget = Budget.from_storage_unit(
                self,
                time_unit="D",
                description=self.name,
                imbalance_fatal=(budget_type == "strict"),
            )
        else:
            raise ValueError(f"Illegal budget_type: {budget_type}")

        return

    def advance(self):
        """
        Advance the storage unit in time.

        Returns:
            None
        """
        if self._itime_step >= self.control.itime_step:
            if self.verbose:
                msg = (
                    f"{self.name} did not advance because "
                    f"it is not behind control time"
                )
                # warn(msg)
                print(msg)  # can/howto make warn flush in real time?
                # is a warning sufficient? an error
            return

        self._advance_variables()
        self._advance_inputs()
        self._itime_step += 1
        return

    def calculate(self, time_length: float, **kwargs) -> None:
        """Calculate storageUnit terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None
        """
        # self._calculate must be implemented by the subclass
        self._calculate(time_length, *kwargs)

        # move to a timestep finalization method at some future date.
        if self.budget is not None:
            self.budget.advance()
            self.budget.calculate()

        return

    def get_metadata(self):
        self.var_meta = self.control.meta.get_vars(self.variables)
        self.input_meta = self.control.meta.get_vars(self.inputs)
        self.param_meta = self.control.meta.get_params(self.parameters)

        # This a hack as we are mushing dims into params. time dimension
        # is also not currently handled at all. Probably need to define dimensions
        # on StorageUnits
        dims = set(self.parameters).difference(set(self.param_meta.keys()))
        self.param_meta = self.control.meta.get_dims(dims)

        if self.verbose:
            from pprint import pprint

            print("Metadata print out for StorageUnit subclass {self.name} :")
            print(f"\n\nParameters ({len(self.param_meta.keys())}): {'*'* 70}")
            pprint(self.param_meta)
            print(f"\n\nInputs ({len(self.input_meta.keys())}): {'*'* 70}")
            pprint(self.input_meta)
            print(f"\n\nVariables ({len(self.var_meta.keys())}): {'*'* 70}")
            pprint(self.var_meta)

        return

    def output_to_csv(self, pth):
        """
        Save each output variable to separate csv file in specified path

        """
        output_data = self.get_output_dataframes()
        for key in output_data:
            df = output_data[key]
            fname = os.path.join(pth, f"{key}.csv")
            df.to_csv(fname)
        return

    def initialize_netcdf(
        self,
        name: str,
        separate_files: bool = True,
    ) -> None:
        """Initialize

        Args:
            name: base directory path or NetCDF file path if separate_files
                is True
            separate_files: boolean indicating if storage component output
                variables should be written to a separate file for each
                variable

        Returns:
            None

        """
        self._output_netcdf = True
        self._netcdf = {}
        if separate_files:
            self._separate_netcdf = True
            # make working directory
            working_path = pl.Path(name)
            working_path.mkdir(parents=True, exist_ok=True)
            for variable_name in self.variables:
                nc_path = pl.Path(working_path) / f"{variable_name}.nc"
                self._netcdf[variable_name] = NetCdfWrite(
                    nc_path,
                    self.params.nhm_coordinates,
                    [variable_name],
                    {variable_name: self.var_meta[variable_name]},
                )
        else:
            initial_variable = self.variables[0]
            pl.Path(name).mkdir(parents=True, exist_ok=True)
            self._netcdf[initial_variable] = NetCdfWrite(
                name,
                self.params.nhm_coordinates,
                self.variables,
                self.var_meta,
            )
            for variable in self.variables[1:]:
                self._netcdf[variable] = self._netcdf[initial_variable]
        return

    def __output_netcdf(self) -> None:
        """Output variable data for a time step

        Returns:
            None

        """
        if self._output_netcdf:
            for idx, variable in enumerate(self.variables):
                if idx == 0 or self._separate_netcdf:
                    self._netcdf[variable].add_simulation_time(
                        self._itime_step,
                        self._simulation_time,
                    )
                self._netcdf[variable].add_data(
                    variable,
                    self._itime_step,
                    getattr(self, variable),
                )
        return

    def _finalize_netcdf(self) -> None:
        if self._output_netcdf:
            for idx, variable in enumerate(self.variables):
                self._netcdf[variable].close()
                if not self._separate_netcdf:
                    break
        return
