import os
import pathlib as pl

import pandas as pd

from ..atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from ..utils.netcdf_utils import NetCdfWrite
from ..utils.parameters import PrmsParameters


class StorageUnit:
    def __init__(
        self,
        storage_type,
        id: list,
        params: PrmsParameters,
        atm: NHMBoundaryLayer,
        verbose: bool,
    ):

        self.storage_type = storage_type
        self.id = id
        self.params = params
        self.atm = atm
        self.verbose = verbose
        self._simulation_time = 0.0

        # netcdf output variables
        self._output_netcdf = False
        self._netcdf = None
        self._separate_netcdf = True
        self._itime_step = -1

        # Go through list of parameters for this process and assign them
        # to self with a value of None
        for param in self.get_required_parameters():
            setattr(self, param, None)

        # Go through the parameters for this process and see if self
        # has a variable with that name.  If so, then assign the parameter
        # value to the self variable.
        for key in params.parameters:
            if hasattr(self, key):
                setattr(self, key, params.parameters[key])

        # if any of the required parameters are still none,
        # then we should terminate with an error
        for key in self.get_required_parameters():
            value = getattr(self, key)
            if value is None:
                print(
                    f"{storage_type} storage unit requires {key} but it was not found in parameters."
                )

        return

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
    def get_required_parameters() -> list:
        raise Exception("This must be overridden")

    def initialize_output_data(self):
        self.output_column_names = ["date"]
        if "nhm_id" in self.params.parameters:
            self.output_column_names += [
                f"nhru_hru_intcpstor_{nhmid}"
                for nhmid in self.params.parameters["nhm_id"]
            ]
        else:
            self.output_column_names += [
                f"intcpstor_{id}" for id in range(self.nhru)
            ]

        # Initialize the output_data dictionary for each entry with an empty list
        self.output_data = {}
        for output_name in self.output_data_names:
            self.output_data[output_name] = []
        return

    def get_output_dataframes(self):
        """
        Return a dictionary of data frames of output data from this process

        """
        output_data = {}
        for key in self.output_data:
            value = self.output_data[key]
            df = pd.DataFrame(value, columns=self.output_column_names)
            df.set_index("date", inplace=True)
            output_data[key] = df
        return output_data

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
            for variable in self.get_output_variables():
                nc_path = pl.Path(working_path) / f"{variable}.nc"
                self._netcdf[variable] = NetCdfWrite(
                    nc_path,
                    self.id,
                    [variable],
                )
        else:
            initial_variable = self.get_output_variables()[0]
            pl.Path(name).mkdir(parents=True, exist_ok=True)
            self._netcdf[initial_variable] = NetCdfWrite(
                name,
                self.id,
                self.get_output_variables(),
            )
            for variable in self.get_output_variables()[1:]:
                self._netcdf[variable] = self._netcdf[initial_variable]

    def __output_netcdf(self) -> None:
        """Output variable data for a time step

        Returns:
            None

        """
        if self._output_netcdf:
            for idx, variable in enumerate(self.get_output_variables()):
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
            for idx, variable in enumerate(self.get_output_variables()):
                self._netcdf[variable].close()
                if not self._separate_netcdf:
                    break
