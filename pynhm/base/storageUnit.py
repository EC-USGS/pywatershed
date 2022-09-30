import os
import pathlib as pl
from warnings import warn

import numpy as np

from pynhm.base import meta
from pynhm.base.budget import Budget

from ..base.adapter import Adapter, adapter_factory
from ..base.timeseries import TimeseriesArray
from ..utils.netcdf_utils import NetCdfWrite
from .accessor import Accessor
from .control import Control


class StorageUnit(Accessor):
    """StorageUnit base class

    The StorageUnit is a base class for conserving mass and energy.

    It has budgets that can optionally be established for mass an energy and
    these can be enforced or simply diagnosed with the model run.

    The  class aims to describe itself through it sstaticmethods and
    properties.

    Conventions are adopted through the use of the following
    properties/methods:

        inputs/get_inputs():
            List the names of variables required from external sources.
            Still working on conventions if these are to be modified but
            the storageUnit. For an input to be successfully inicluded,
            that variable must be defined in the metadata
            (pynhm/static/metadata/variables.yaml).
            Efforts should be made to not use diagnostic variables as input
            as much as possible.
        variables/get_variables():
            List the names of internal public variables. If not necessary,
            to be public, variables should be made private with a single,
            leading underscore and not maintained in this list. For an input
            to be successfully inicluded, that variable must be defined in the
            metadata (pynhm/static/metadata/variables.yaml).
            Efforts should be made not to track diagnostic variables in this
            public variable set, as much as possible.
        parameters/get_parameters():
            List the names of parameters used by the subclass.
        description:
            Return a dictionary of with the storage unit cubclass name and its
            metadata for all variables for each of inputs, variables, and
            parameters.
        get_init_values:
            Return a dictionary of initialization values for variables.
            Note that these may be overridden by subclass initialization
            routines (e.g. using parameters) or by restart values. So these
            are not "initial values", they are initialization values that
            are set when the variable is declared from metadata in
            _initialize_var(). Initization values should be nan as much as
            possible.
        mass_budget_terms/get_mass_budget_terms():
            These terms must all in in the same units across all components of
            the budget (inputs, outputs, storage_changes). Diagnostic variables
            should not appear in the budget terms, only prognostic variables
            should.
        _advance_variables():
            This advance should exactly specify the prognostic variables in
            setting previous values to current values. When/if necessary to
            keep previous diagnostic variables, those must not appear here but
            in _calculate().
        _calculate():
            This method is to be overridden by the subclass. Near the end of
            the method, the subclass should calculate its changes in mass and
            energy storage in an obvious way. As commented for
            mass_budget_terms, storage changes should only be tracked for
            prognostic variables. (For example is snow_water_equiv = snow_ice +
            snow_liquid, then storage changes for snow_ice and snow_liquid
            should be tracked and not for snow_water_equiv).

    """

    def __init__(
        self,
        control: Control,
        verbose: bool,
    ):

        self.name = "StorageUnit"
        self.control = control
        self.params = self.control.params
        self.verbose = verbose

        # netcdf output variables
        self._output_netcdf = False
        self._netcdf = None
        self._separate_netcdf = True
        self._itime_step = -1

        self._get_metadata()
        self._initialize_self_variables()
        self._set_initial_conditions()

        return None

    def output(self) -> None:
        """Output data to previously initialized output types.

        Writes output for initalized output types.

        Returns:
            None

        """
        if self._output_netcdf:
            if self.verbose:
                print(f"writing output for: {self.name}")
            self.__output_netcdf()

        if self.budget is not None:
            self.budget.output()
        return

    def finalize(self) -> None:
        """Finalize storageUnit

        Finalizes the object, including output methods.

        Returns:
            None

        """
        if self.verbose:
            print(f"finalizing: {self.name}")

        self._finalize_netcdf()
        if self.budget is not None:
            self.budget._finalize_netcdf()
        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get a tuple of parameter names."""
        raise Exception("This must be overridden")

    @staticmethod
    def get_inputs() -> tuple:
        """Get a tuple of input variable names."""
        raise Exception("This must be overridden")

    @classmethod
    def get_variables(cls) -> tuple:
        """Get a tuple of internal public variable names."""
        return list(cls.get_init_values().keys())

    @classmethod
    def get_mass_budget_terms(cls) -> dict:
        """Get a dictionary of variable names for mass budget terms."""
        mass_budget_terms = {
            "inputs": list(
                meta.filter_vars(
                    cls.get_inputs(), "var_category", "mass flux"
                ).keys()
            ),
            "outputs": list(
                meta.filter_vars(
                    cls.get_variables(), "var_category", "mass flux"
                ).keys()
            ),
            "storage_changes": list(
                meta.filter_vars(
                    cls.get_variables(), "var_category", "mass storage change"
                ).keys()
            ),
        }
        return mass_budget_terms

    @classmethod
    def description(cls) -> dict:
        """A description (all metadata) for all variables in inputs, variables,
        and parameters."""
        return {
            "class_name": cls.__name__,
            "mass_budget_terms": cls.get_mass_budget_terms(),
            "inputs": meta.get_vars(cls.get_inputs()),
            "variables": meta.get_vars(cls.get_variables()),
            "parameters": meta.get_params(cls.get_parameters()),
        }

    @staticmethod
    def get_restart_variables() -> list:
        """Get a list of restart varible names."""
        raise Exception("This must be overridden")

    @staticmethod
    def get_init_values() -> dict:
        """Get a dictionary of initialization values for each public
        variable."""
        raise Exception("This must be overridden")

    @property
    def parameters(self) -> tuple:
        """A tuple of parameter names."""
        return self.get_parameters()

    @property
    def inputs(self) -> tuple:
        """A tuple of input variable names."""
        return self.get_inputs()

    @property
    def variables(self) -> tuple:
        """A tuple of public variable names."""
        return self.get_variables()

    @property
    def restart_variables(self) -> tuple:
        """A tuple of restart variable names."""
        return self.get_restart_variables()

    @property
    def init_values(self) -> dict:
        """A dictionary of initial values for each public variable."""
        return self.get_init_values()

    @property
    def mass_budget_terms(self) -> dict:
        """A dictionary of variable names for the mass budget terms."""
        return self.get_mass_budget_terms()

    def _initialize_self_variables(self, restart: bool = False):
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
            self._initialize_var(name)

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

    def _initialize_var(self, var_name: str, flt_to_dbl: bool = True):
        """Initialize a variable using get_init_values and metadata.

        Initialized variables can be for single time or they can be a timeries
        array object if they have a time dimension in metadata.
        """
        init_vals = self.get_init_values()
        if var_name not in init_vals.keys():
            if self.verbose:
                warn(
                    f"{var_name} not initialized (no initial value specified)"
                )
            return

        dims = [
            self[vv] for vv in self.var_meta[var_name]["dimensions"].values()
        ]
        init_type = self.var_meta[var_name]["type"]

        if len(dims) == 1:
            self[var_name] = np.full(
                dims, init_vals[var_name], dtype=init_type
            )
        else:
            self[var_name] = TimeseriesArray(
                var_name=var_name,
                control=self.control,
                array=np.full(dims, init_vals[var_name], dtype=init_type),
                time=self._time,
            )
        return

    def _set_initial_conditions(self):
        raise Exception("This must be overridden")

    def _advance_variables(self):
        raise Exception("This must be overridden")

    def _advance_inputs(self):
        for key, value in self._input_variables_dict.items():
            value.advance()
            self[key][:] = value.current

        return

    def _set_inputs(self, args):
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

    def _set_budget(self, behavior, basis: str = "unit"):
        if behavior is None:
            self.budget = None
        elif behavior in ["error", "warn"]:
            self.budget = Budget.from_storage_unit(
                self,
                time_unit="D",
                description=self.name,
                imbalance_fatal=(behavior == "error"),
                basis=basis,
            )
        else:
            raise ValueError(f"Illegal behavior: {behavior}")

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

        if self.verbose:
            print(f"advancing: {self.name}")

        self._advance_variables()
        self._advance_inputs()
        self._itime_step += 1
        return

    def _calculate(self):
        raise Exception("This must be overridden")

    def calculate(self, time_length: float, **kwargs) -> None:
        """Calculate storageUnit terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None
        """
        if self.verbose:
            print(f"calculating: {self.name}")

        # self._calculate must be implemented by the subclass
        self._calculate(time_length, *kwargs)

        # move to a timestep finalization method at some future date.
        if self.budget is not None:
            self.budget.advance()
            self.budget.calculate()

        return

    def _get_metadata(self):
        self.var_meta = self.control.meta.get_vars(self.variables)
        self.input_meta = self.control.meta.get_vars(self.inputs)
        self.param_meta = self.control.meta.get_params(self.parameters)

        # This a hack as we are mushing dims into params. time dimension
        # is also not currently handled at all. Probably need to define
        # dimensions on StorageUnits
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
        output_dir: str,
        separate_files: bool = True,
        budget_args: dict = {},
    ) -> None:
        """Initialize NetCDF output.

        Args:
            output_dir: base directory path or NetCDF file path if separate_files
                is True
            separate_files: boolean indicating if storage component output
                variables should be written to a separate file for each
                variable
            budget_args: a dict of argument key: values to pass to
                initialize_netcdf on this storage unit's budget. see budget
                object for options.

        Returns:
            None

        """
        if self.verbose:
            print(f"initializing netcdf output for: {self.output_dir}")

        self._output_netcdf = True
        self._netcdf = {}
        if separate_files:
            self._separate_netcdf = True
            # make working directory
            output_dir = pl.Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            for variable_name in self.variables:
                nc_path = pl.Path(output_dir) / f"{variable_name}.nc"
                self._netcdf[variable_name] = NetCdfWrite(
                    nc_path,
                    self.params.nhm_coordinates,
                    [variable_name],
                    {variable_name: self.var_meta[variable_name]},
                )
        else:
            initial_variable = self.variables[0]
            pl.Path(output_dir).mkdir(parents=True, exist_ok=True)
            self._netcdf[initial_variable] = NetCdfWrite(
                output_dir,
                self.params.nhm_coordinates,
                self.variables,
                self.var_meta,
            )
            # JLM: this code seems unnecessary/overwrought
            # we are currently not testing single file output.
            for variable in self.variables[1:]:
                self._netcdf[variable] = self._netcdf[initial_variable]

        if self.budget is not None:
            if "output_dir" not in budget_args.keys():
                budget_args["output_dir"] = output_dir
            self.budget.initialize_netcdf(**budget_args)

        return

    def __output_netcdf(self) -> None:
        """Output variable data to NetCDF for a time step.

        Returns:
            None

        """
        if self._output_netcdf:
            for idx, variable in enumerate(self.variables):
                if idx == 0 or self._separate_netcdf:
                    self._netcdf[variable].add_simulation_time(
                        self.control.itime_step, self.control.current_datetime
                    )
                self._netcdf[variable].add_data(
                    variable,
                    self._itime_step,
                    getattr(self, variable),
                )

        return

    def _finalize_netcdf(self) -> None:
        """Finalize NetCDF output to disk.

        Returns:
            None
        """
        if self._output_netcdf:
            for idx, variable in enumerate(self.variables):
                self._netcdf[variable].close()
                if not self._separate_netcdf:
                    break
        return
