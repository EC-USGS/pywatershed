import inspect
import os
import pathlib as pl
from typing import Literal
from warnings import warn

import numpy as np

from ..base import meta
from ..base.adapter import Adapter, adapter_factory
from ..base.data_model import _merge_dicts
from ..base.timeseries import TimeseriesArray
from ..parameters import Parameters
from ..utils.netcdf_utils import NetCdfWrite
from .accessor import Accessor
from .control import Control


class Process(Accessor):
    """Base class for physical process representation.

    The  class aims to describe itself through its staticmethods and
    properties.

    Conventions are adopted through the use of the following
    properties/methods:

    inputs/get_inputs():
        List the names of variables required from external sources.
        Still working on conventions if these are to be modified.
        For an input to be successfully inicluded,
        that variable must be defined in the metadata
        (pywatershed/static/metadata/variables.yaml).
        Efforts should be made to not use diagnostic variables as input
        as much as possible.
    variables/get_variables():
        List the names of internal public variables. If not necessary,
        to be public, variables should be made private with a single,
        leading underscore and not maintained in this list. For an input
        to be successfully inicluded, that variable must be defined in the
        metadata (pywatershed/static/metadata/variables.yaml).
        Efforts should be made not to track diagnostic variables in this
        public variable set, as much as possible.
    parameters/get_parameters():
        List the names of parameters used by the subclass.
    description:
        Return a dictionary of with the process subclass name and its
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

    See Also
    --------
    pywatershed.base.ConservativeProcess

    Args
    ----
    control:
        A Control object
    discretization:
        A discretization object
    parameters:
        The parameters for this object
    metadata_patches:
        Override static metadata for any public parameter or variable --
        experimental.
    metadata_patch_conflicts:
        How to handle metadata_patches conflicts. Experimental.
    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        metadata_patches: dict[dict] = None,
        metadata_patch_conflicts: Literal["ignore", "warn", "error"] = "error",
    ):
        self.name = "Process"
        self.control = control

        # params from dis and params: methodize this
        missing_params = set(self.parameters).difference(
            set(parameters.variables.keys())
        )
        if missing_params:
            self.params = type(parameters).merge(parameters, discretization)
        else:
            self.params = parameters.subset(self.get_parameters())

        # netcdf output variables
        self._netcdf_initialized = False

        self._itime_step = -1

        # TODO metadata patching.
        self._set_metadata()
        if metadata_patches is not None:
            self._patch_metadata(
                metadata_patches,
                conflicts=metadata_patch_conflicts,
            )

        self._initialize_self_variables()
        self._set_initial_conditions()

        return None

    def output(self) -> None:
        """Output data to previously initialized output types.
        Returns:
            None
        """
        if self._netcdf_initialized:
            if self._verbose:
                print(f"writing output for: {self.name}")
            self._output_netcdf()
        return

    def finalize(self) -> None:
        """Finalizes the Process, including output methods.
        Returns:
            None
        """
        if self._verbose:
            print(f"finalizing: {self.name}")

        self._finalize_netcdf()
        return

    @staticmethod
    def get_dimensions() -> tuple:
        """Get a tuple of dimension names for this Process."""
        raise Exception("This must be overridden")

    @staticmethod
    def get_parameters() -> tuple:
        """Get a tuple of parameter names for this Process."""
        raise Exception("This must be overridden")

    @staticmethod
    def get_inputs() -> tuple:
        """Get a tuple of input variable names for this Process."""
        raise Exception("This must be overridden")

    @classmethod
    def get_variables(cls) -> tuple:
        """Get a tuple of (public) variable names for this Process."""
        return list(cls.get_init_values().keys())

    @classmethod
    def description(cls) -> dict:
        """A dictionary description of this Process.

        Returns:
            All metadata for all variables in inputs, variables,
            and parameters."""
        return {
            "class_name": cls.__name__,
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
    def dimensions(self) -> tuple:
        """A tuple of parameter names."""
        return self.get_dimensions()

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

    def _initialize_self_variables(self, restart: bool = False):
        # dims
        for name in self.dimensions:
            if name == "ntime":
                setattr(self, name, self.control.n_times)
            else:
                setattr(self, name, self.params.dims[name])

        # parameters
        for name in self.parameters:
            setattr(self, name, self.params.get_param_values(name))

        # inputs
        for name in self.inputs:
            # dims of internal variables never have time, so they are spatial
            spatial_dims = self.params.get_dim_values(
                list(meta.find_variables(name)[name]["dims"])
            )
            spatial_dims = tuple(spatial_dims.values())
            setattr(self, name, np.zeros(spatial_dims, dtype=float) + np.nan)

        # variables
        # skip restart variables if restart (for speed) ?
        # the code is below but commented.
        # restart_variables = self.restart_variables
        for name in self.variables:
            # if restart and (name in restart_variables):
            #     continue
            self._initialize_var(name)

        return

    def _initialize_var(self, var_name: str, flt_to_dbl: bool = True):
        """Initialize a variable using get_init_values and metadata.

        Initialized variables can be for single time or they can be a timeries
        array object if they have a time dimension in metadata.
        """
        init_vals = self.get_init_values()
        if var_name not in init_vals.keys():
            if self._verbose:
                warn(
                    f"{var_name} not initialized (no initial value specified)"
                )
            return

        dims = [self[vv] for vv in self.meta[var_name]["dims"]]
        init_type = self.meta[var_name]["type"]

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
            ii_dims = self.control.meta.get_dimensions(ii)[ii]
            # This accomodates Timeseries like objects that need to init
            # both full rank and reduced rank versions of their data
            # this is pretty adhoc
            check_list = ["time", "doy"]
            if len([mm for mm in check_list if mm in ii_dims[0]]):
                ii_dims = ii_dims[1:]

            ii_dim_sizes = tuple(self.params.get_dim_values(ii_dims).values())
            ii_type = self.control.meta.get_numpy_types(ii)[ii]

            self._input_variables_dict[ii] = adapter_factory(
                args[ii],
                variable_name=ii,
                control=args["control"],
                variable_dim_sizes=ii_dim_sizes,
                variable_type=ii_type,
            )
            if self._input_variables_dict[ii]:
                self[ii] = self._input_variables_dict[ii].current

        return

    def _set_options(self, init_locals):
        """Set options options on self if supplied on init, else take
        from control"""
        # some self and Process introspection reveals the option names
        init_arg_names = set(
            inspect.signature(self.__init__).parameters.keys()
        )
        process_init_args = set(
            inspect.signature(Process.__init__).parameters.keys()
        )
        inputs_args = self.inputs
        non_option_args = process_init_args.union(inputs_args)
        option_names = init_arg_names.difference(non_option_args)
        for opt in option_names:
            if opt in init_locals.keys() and init_locals[opt] is not None:
                setattr(self, f"_{opt}", init_locals[opt])
            elif opt in self.control.options.keys():
                setattr(self, f"_{opt}", self.control.options[opt])
            else:
                setattr(self, f"_{opt}", None)

        return

    def set_input_to_adapter(self, input_variable_name: str, adapter: Adapter):
        # Todo: document
        self._input_variables_dict[input_variable_name] = adapter
        # can NOT use [:] on the LHS as we are relying on pointers between
        # boxes. [:] on the LHS here means it's not a pointer and then
        # requires that the calculation of the input happens before the
        # advance of this process.
        self[input_variable_name] = adapter.current
        return

    def advance(self):
        """
        Advance the Process in time.

        Returns:
            None
        """

        if self._itime_step >= self.control.itime_step:
            if self._verbose:
                msg = (
                    f"{self.name} did not advance because "
                    f"it is not behind control time"
                )
                # warn(msg)
                print(msg)  # can/howto make warn flush in real time?
                # is a warning sufficient? an error
            return

        if self._verbose:
            print(f"advancing: {self.name}")

        self._advance_variables()
        self._advance_inputs()
        self._itime_step += 1
        return

    def _calculate(self):
        raise Exception("This must be overridden")

    def calculate(self, time_length: float, **kwargs) -> None:
        """Calculate Process terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None
        """
        if self._verbose:
            print(f"calculating: {self.name}")

        # self._calculate must be implemented by the subclass
        self._calculate(time_length, *kwargs)

        return

    def _set_metadata(self):
        """Set metadata on self for self's inputs, parameters, and variables"""
        meta_keys = (*self.variables, *self.inputs, *self.parameters)
        msg = (
            "Duplicate varible names amongst self's variables, "
            "inputs, and parameters"
        )
        assert len(meta_keys) == len(self.variables) + len(self.inputs) + len(
            self.parameters
        ), msg

        self.meta = self.control.meta.find_variables(meta_keys)
        if "global" not in self.meta.keys():
            self.meta["global"] = {}
        return

    def _patch_metadata(
        self, patches, conflicts: Literal["left", "warn", "error"] = "error"
    ):
        patch_meta_on_self = {
            kk: vv for kk, vv in patches.items() if kk in self.meta.keys()
        }
        self.meta = _merge_dicts(
            [self.meta, patch_meta_on_self], conflicts=conflicts
        )
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
        output_dir: [str, pl.Path] = None,
        separate_files: bool = None,
        output_vars: list = None,
    ) -> None:
        """Initialize NetCDF output.

        Args:
            output_dir: base directory path or NetCDF file path if
                separate_files is True
            separate_files: boolean indicating if storage component output
                variables should be written to a separate file for each
                variable
            output_vars: list of variable names to outuput.

        Returns:
            None

        """
        if self._netcdf_initialized:
            msg = (
                f"{self.name} class previously initialized netcdf output "
                f"in {self._netcdf_output_dir}"
            )
            warn(msg)
            return

        if self._verbose:
            print(f"initializing netcdf output for: {self.name}")

        (
            output_dir,
            output_vars,
            separate_files,
        ) = self._reconcile_nc_args_w_control_opts(
            output_dir, output_vars, separate_files
        )

        # apply defaults if necessary
        if output_dir is None:
            msg = (
                "An output directory is required to be specified for netcdf"
                "initialization."
            )
            raise ValueError(msg)

        if separate_files is None:
            separate_files = True
        self._netcdf_separate = separate_files

        self._netcdf_initialized = True
        self._netcdf_output_dir = pl.Path(output_dir)
        if output_vars is None:
            self._netcdf_output_vars = self.variables
        else:
            self._netcdf_output_vars = list(
                set(output_vars).intersection(set(self.variables))
            )
            if len(self._netcdf_output_vars) == 0:
                self._netcdf_initialized = False
                return

        self._netcdf = {}

        if self._netcdf_separate:
            # make working directory
            self._netcdf_output_dir.mkdir(parents=True, exist_ok=True)
            for variable_name in self.variables:
                if (self._netcdf_output_vars is not None) and (
                    variable_name not in self._netcdf_output_vars
                ):
                    continue
                nc_path = self._netcdf_output_dir / f"{variable_name}.nc"
                self._netcdf[variable_name] = NetCdfWrite(
                    nc_path,
                    self.params.coords,
                    [variable_name],
                    {variable_name: self.meta[variable_name]},
                    {"process class": self.name},
                )

        else:
            if self._netcdf_output_vars is None:
                the_out_vars = self.variables
            else:
                the_out_vars = self._netcdf_output_vars

            initial_variable = the_out_vars[0]
            self._netcdf_output_dir.mkdir(parents=True, exist_ok=True)
            self._netcdf[initial_variable] = NetCdfWrite(
                self._netcdf_output_dir / f"{self.name}.nc",
                self.params.coords,
                self._netcdf_output_vars,
                self.meta,
                {"process class": self.name},
            )
            for variable in the_out_vars[1:]:
                self._netcdf[variable] = self._netcdf[initial_variable]

        return

    def _output_netcdf(self) -> None:
        """Output variable data to NetCDF for a time step.

        Returns:
            None

        """
        if self._netcdf_initialized:
            time_added = False
            for variable in self.variables:
                if (self._netcdf_output_vars is not None) and (
                    variable not in self._netcdf_output_vars
                ):
                    continue
                if not time_added or self._netcdf_separate:
                    time_added = True
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
        if self._netcdf_initialized:
            for idx, variable in enumerate(self.variables):
                if (self._netcdf_output_vars is not None) and (
                    variable not in self._netcdf_output_vars
                ):
                    continue

                self._netcdf[variable].close()
                if not self._netcdf_separate:
                    break
        return

    def _reconcile_nc_args_w_control_opts(
        self, output_dir, output_vars, separate_files
    ):
        # can treat the other args but they are not yet in the available opts
        arg_opt_name_map = {
            "output_dir": "netcdf_output_dir",
            "output_vars": "netcdf_output_var_names",
            "separate_files": "netcdf_output_separate_files",
        }

        args = {
            "output_dir": output_dir,
            "output_vars": output_vars,
            "separate_files": separate_files,
        }

        for vv in args.keys():
            arg_val = args[vv]
            opt_name = arg_opt_name_map[vv]
            opts = self.control.options
            if opt_name in opts.keys():
                opt_val = opts[opt_name]
            else:
                opt_val = None

            # if options is a list (output_vars), take it's intersection with
            # self.variables

            # set value into args to return to return

            # the 4 cases:
            if opt_val is None and arg_val is None:
                pass

            elif opt_val is None:
                pass

            elif arg_val is None:
                args[vv] = opt_val

            elif opt_val is not None and arg_val is not None:
                if opt_val == arg_val:
                    pass
                else:
                    msg = (
                        f"control.option '{opt_name}' conflicts with "
                        f"initialize_netcdf() argument {vv}"
                    )
                    raise ValueError(msg)

        return args["output_dir"], args["output_vars"], args["separate_files"]
