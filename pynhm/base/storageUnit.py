import os
import pathlib as pl
from warnings import warn

import numpy as np
import pandas as pd

from ..atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from ..constants import nan, one, zero
from ..utils.netcdf_utils import NetCdfWrite
from ..utils.parameters import PrmsParameters
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
        params: PrmsParameters,
        verbose: bool,
        subclass_name="StorageUnit",
    ):

        self.name = subclass_name
        self.control = control
        self.params = params
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

    def calculate(self, simulation_time: float) -> None:
        """Calculate storageUnit terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """
        raise NotImplementedError("Override in child class.")

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
        # skip restart variables if restart (for speed) ? the code is below but commented.
        # restart_variables = self.restart_variables
        for name in self.parameters:
            setattr(self, name, self.params.parameters[name])
        for name in self.inputs:
            setattr(self, name, np.zeros(self.nhru, dtype=float) + np.nan)
        for name in self.variables:
            # if restart and (name in restart_variables):
            #     continue
            self.initialize_var(name)
        return

    def initialize_var(self, var_name):
        if self.input_meta[var_name]["dimensions"][0] == "nsegment":
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
            value.advance()  # (self.control.itime_step)
            v = getattr(self, key)
            v[:] = value.current
        return

    def advance(self):
        """
        Advance the storage unit in time.

        Returns:
            None
        """
        if self._itime_step >= self.control.itime_step:
            if self.verbose:
                msg = f"{self.name} did not advance because it is not behind control time"
                # warn(msg)
                print(msg)  # can/howto make warn flush in real time?
                # is a warning sufficient? an error
            return

        self._advance_variables()
        self._advance_inputs()
        self._itime_step += 1
        return

    def get_metadata(self):
        self.var_meta = self.control.meta.get_vars(self.variables)
        self.input_meta = self.control.meta.get_vars(self.variables)
        self.param_meta = self.control.meta.get_params(self.parameters)

        # This a hack as we are mushing dims into params. time dimension
        # is also not currently handled at all. Probably need to define dimensions
        # on StorageUnits
        dims = set(self.parameters).difference(set(self.param_meta.keys()))
        self.param_meta = self.control.meta.get_dims(dims)

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
