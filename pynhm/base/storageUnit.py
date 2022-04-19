import os
import pathlib as pl

import numpy as np
import pandas as pd

from ..atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from ..utils.netcdf_utils import NetCdfWrite
from ..utils.parameters import PrmsParameters
from .accessor import Accessor
from .control import Control


class StorageUnit(Accessor):
    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        verbose: bool,
    ):

        self.name = "Storage Unit"
        self.control = control
        self.params = params
        self.verbose = verbose
        self._simulation_time = 0.0

        # netcdf output variables
        self._output_netcdf = False
        self._netcdf = None
        self._separate_netcdf = True
        self._itime_step = -1  # JLM CHECK

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

    def finalize(self) -> None:
        """Finalize storageUnit

        Returns:
            None

        """
        self._finalize_netcdf()

    @staticmethod
    def get_parameters() -> list:
        raise Exception("This must be overridden")

    @staticmethod
    def get_inputs() -> list:
        raise Exception("This must be overridden")

    @staticmethod
    def get_variables() -> list:
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

    def initialize_self_variables(self):
        # todo: get the type from metadata
        for name in self.parameters:
            setattr(self, name, self.params.parameters[name])
        for name in self.variables:
            setattr(self, name, np.zeros(self.nhru, dtype=float))  # + np.nan)
        for name in self.inputs:
            setattr(self, name, np.zeros(self.nhru, dtype=float))  # + np.nan)

    def set_initial_conditons(self):
        raise Exception("This must be overridden")

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
            for variable in self.variables:
                nc_path = pl.Path(working_path) / f"{variable}.nc"
                self._netcdf[variable] = NetCdfWrite(
                    nc_path,
                    self.params.nhm_coordinate,
                    [variable],
                )
        else:
            initial_variable = self.variables[0]
            pl.Path(name).mkdir(parents=True, exist_ok=True)
            self._netcdf[initial_variable] = NetCdfWrite(
                name,
                self.id,
                self.variables,
            )
            for variable in self.variables[1:]:
                self._netcdf[variable] = self._netcdf[initial_variable]

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

    def _finalize_netcdf(self) -> None:
        if self._output_netcdf:
            for idx, variable in enumerate(self.variables):
                self._netcdf[variable].close()
                if not self._separate_netcdf:
                    break
