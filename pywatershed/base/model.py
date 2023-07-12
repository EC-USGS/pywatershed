import pathlib as pl
from copy import deepcopy
from pprint import pprint
from typing import Union
from warnings import warn

from tqdm.auto import tqdm

from ..base.adapter import adapter_factory
from ..base.control import Control
from ..constants import fileish
from ..parameters import Parameters, PrmsParameters
from ..utils.path import path_rel_to_yml

# This is a convenience
process_order_nhm = [
    "PRMSSolarGeometry",
    "PRMSAtmosphere",
    "PRMSCanopy",
    "PRMSSnow",
    "PRMSRunoff",
    "PRMSSoilzone",
    "PRMSEt",
    "PRMSGroundwater",
    "PRMSChannel",
]


class Model:
    """Build a model in pywatershed.

    This is the class that helps execute sets of Processes in order.

    There are two distinct ways of instatniating the Model class described
    below: 1) PRMS-legacy instantation, 2) pywatershed-centric instatiation.

    Args:
        process_list_or_model_dict: a process list or a model dictionary,
            see ways of instatiation below.
        control: Control object or None, see below.
        parameters: A Parameters object of a dictionary of Parameters objects.
            Default is None. See ways of instatnation below.
        find_input_files: Search/find input file on __init__ or delay until run
           or advance of the model. Delaying (False) allows ModelGraph of the
           specified model without the need for input files.

    PRMS-legacy instantiation
    -----------------------------

    This method of instantiation may be deprecated in the future. It is
    intended to support PRMS files with few modifications.

    Note: Pywatershed v0.2.0 has a minor break to backwards compatibility (pre
    v0.2.0) where the list of processes is not passed as a list rather than an
    upacked list or a number of arguments.

    Old way, first three args:

    - **process_list_or_model_dict** - A process list of PRMS model components
    - **control** - A control object
    - **parameters** - A PrmsParameters object

    See first example below for more details.
    There are variations of this where

    pywatershed-centric instatiation
    ------------------------------------

    The new/pywatershed way was introduced to handle more aribtrary and subdry
    models. It is loosely based on how MF6 structures its inputs. All of the
    work is contained in the first argument, the model dictionary, and the
    control and parameters are both None.

    New way, only one of first three args:

    - **process_list_or_model_dict** - a "model dictionary", detailed below.
    - (*control* - None is default)
    - (*parameters* - None is default)

    This means that the user must understand the model dictionary supplied. A
    model dictionary has aribtrary keys but prescribed values. Values may be
    the following kinds of objects: Control, discretization dictionary, process
    dictionary, process order list, and exchanges dictionaries. (Exchanges are
    not yet supported). Each of these is described below. Model dictionaries
    can also be specfied via yaml files. Notes on yaml versus in-memory
    requirements for the values are described as they are (necessarily)
    different.

    Model dictionary values description:
    ====================================

    - **control** - The control object supplies two things for the model 1) the
      global time discretization (start and stop times, as well as time
      step), and 2) default global options for all model processes.
      Only one control object can be included in the model dictionary. Though
      the key for the control can be arbitrary, the value is either an instance
      of class Control or, in the case of a yaml model dictionary, a control
      yaml file to be loaded by Control.from_yml() (todo: link to this
      staticmethod).
    - **discretizations** - Multiple discretizations may be supplied to the
      model dictionary, each with arbitrary names. These provide spatial
      discretization information which may be shared by multiple processes
      specified later in `process_list`. Each process will refer to its
      required discretization using the key for that discretization in the
      model dict. The value of each is an instance of the class Parameters. In
      a yaml model dictionary, the value of each discretization is a path
      to a netcdf file for that discretization.
    - **processes** - Multiple processes with arbitrary keys are to be
      supplied to model dictionary (if running just one process, there's no
      need for the Model class). The values for each process are dictionaries
      and have the required keys:

      - *class* - The desired class, in the case of a yaml model dictionary,
        the string name of a class in the pywatershed namespace
      - *parameters* - Either a Parameters object or, in the case of a yaml
        model dictionary, the path of a netcdf file with the parameters
      - *dis* - The key of the dis to use as specified at the top level of
        model dictionary.

    - **model_order** - A list of the processess *keys* in the order the model
      is to be executed.
    - **exchanges** - Future.

    Examples:
    ---------

    These examples will work if you have an editable install of the repository
    (not installed from pypi).

    Construct a PRMS-legacy based model:

    ..
        import pywatershed as pws
        test_data_dir = pws.constants.__pywatershed_root__ / "../test_data"
        domain_dir = test_data_dir / "drb_2yr"
        # A PRMS-native control file
        control_file = domain_dir / "control.test"
        # PRMS-native parameter file
        parameter_file = domain_dir / "myparam.param"
        control = pws.Control.load(control_file)
        control.options['input_dir'] = domain_dir / "output"
        params = pws.parameters.PrmsParameters.load(parameter_file)
        model_procs = [pws.PRMSGroundwater, pws.PRMSChannel,]
        model = pws.Model(
            model_procs,
            control=control,
            parameters=params,
        )
        model.run()

    >>> import pywatershed as pws
    >>> test_data_dir = pws.constants.__pywatershed_root__ / "../test_data"
    >>> domain_dir = test_data_dir / "drb_2yr"
    >>> # A PRMS-native control file
    >>> control_file = domain_dir / "control.test"
    >>> # PRMS-native parameter file
    >>> parameter_file = domain_dir / "myparam.param"
    >>> control = pws.Control.load(control_file)
    >>> control.options["input_dir"] = domain_dir / "output"
    >>> params = pws.parameters.PrmsParameters.load(parameter_file)
    >>> model_procs = [
    ...     pws.PRMSGroundwater,
    ...     pws.PRMSChannel,
    ... ]
    >>> model = pws.Model(
    ...     model_procs,
    ...     control=control,
    ...     parameters=params,
    ... )
    PRMSGroundwater jit compiling with numba
    PRMSChannel jit compiling with numba
    >>> model.run()
    100%|█████████████████████████████████████████████████████████| 731/731 [00:00<00:00, 1249.26it/s]
    model.run(): finalizing


    Construct a model the pywatershed-centric way, in memory:

    ..
        import pywatershed as pws
        test_data_dir = pws.constants.__pywatershed_root__ / "../test_data"
        domain_dir = test_data_dir / "drb_2yr"
        dis_hru = pws.Parameters.from_netcdf(
            domain_dir / "parameters_dis_hru.nc", encoding=False
        )
        control_file = domain_dir / "control.test"
        control = pws.Control.load(control_file)
        control.options['input_dir'] = domain_dir
        params = {}
        for proc in ["SolarGeometry", "Atmosphere", "Canopy", "Snow"]:
            param_file = domain_dir / f"parameters_PRMS{proc}.nc"
            params[proc.lower()] = pws.Parameters.from_netcdf(param_file)
        model_dict = {
           'control': control,
           'dis_hru': dis_hru,
           'model_order': [
                'prmssolargeometry',
                'prmsatmosphere',
                'prmscanopy',
                'prmssnow',
            ],
            'prmssolargeometry': {
                'class': pws.PRMSSolarGeometry,
                'parameters': params['solargeometry'],
                'dis': 'dis_hru',
            },
            'prmsatmosphere': {
                'class': pws.PRMSAtmosphere,
                'parameters': params['atmosphere'],
                'dis': 'dis_hru',
            },
            'prmscanopy': {
                'class': pws.PRMSCanopy,
                'parameters': params['canopy'],
                'dis': 'dis_hru',
            },
            'prmssnow': {
                'class': pws.PRMSSnow,
                'parameters': params['snow'],
                'dis': 'dis_hru',
            }
        }
        model = pws.Model(model_dict)
        model.run()

    >>> import pywatershed as pws
    >>> test_data_dir = pws.constants.__pywatershed_root__ / "../test_data"
    >>> domain_dir = test_data_dir / "drb_2yr"
    >>> dis_hru = pws.Parameters.from_netcdf(
    ...     domain_dir / "parameters_dis_hru.nc", encoding=False
    ... )
    >>> control_file = domain_dir / "control.test"
    >>> control = pws.Control.load(control_file)
    >>> control.options["input_dir"] = domain_dir
    >>> params = {}
    >>> for proc in ["SolarGeometry", "Atmosphere", "Canopy", "Snow"]:
    ...     param_file = domain_dir / f"parameters_PRMS{proc}.nc"
    ...     params[proc.lower()] = pws.Parameters.from_netcdf(param_file)
    ...
    >>> model_dict = {
    ...     "control": control,
    ...     "dis_hru": dis_hru,
    ...     "model_order": [
    ...         "prmssolargeometry",
    ...         "prmsatmosphere",
    ...         "prmscanopy",
    ...         "prmssnow",
    ...     ],
    ...     "prmssolargeometry": {
    ...         "class": pws.PRMSSolarGeometry,
    ...         "parameters": params["solargeometry"],
    ...         "dis": "dis_hru",
    ...     },
    ...     "prmsatmosphere": {
    ...         "class": pws.PRMSAtmosphere,
    ...         "parameters": params["atmosphere"],
    ...         "dis": "dis_hru",
    ...     },
    ...     "prmscanopy": {
    ...         "class": pws.PRMSCanopy,
    ...         "parameters": params["canopy"],
    ...         "dis": "dis_hru",
    ...     },
    ...     "prmssnow": {
    ...         "class": pws.PRMSSnow,
    ...         "parameters": params["snow"],
    ...         "dis": "dis_hru",
    ...     },
    ... }
    >>> model = pws.Model(model_dict)
    PRMSCanopy jit compiling with numba
    PRMSSnow jit compiling with numba
    >>> model.run()
      0%|                                                                     | 0/731 [00:00<?, ?it/s]
    /Users/jmccreight/usgs/pywatershed2/pywatershed/hydrology/prms_snow.py:1086: NumbaExperimentalFeatureWarning: Use of isinstance() detected. This is an experimental feature.
        through_rain[:] = np.where(wh_through, net_rain, zero)
    100%|███████████████████████████████████████████████████████████| 731/731 [00:07<00:00, 96.76it/s]
    model.run(): finalizing

    Construct a model the pywatershed-centric way, from a yaml file definition:

    ..
        import yaml
        import pywatershed as pws
        test_data_dir = pws.constants.__pywatershed_root__ / "../test_data"
        domain_dir = test_data_dir / "drb_2yr"
        control = {
            "start_time": "1979-01-01T00:00:00",
            "end_time": "1980-12-31T00:00:00",
            "time_step": 24,
            "time_step_units": "h",
            "verbosity": 0,
            "budget_type": "warn",
            "init_vars_from_file": 0,
            "input_dir": str(domain_dir),
        }
        control_file = domain_dir / "example_control.yml"

        model_dict = {
            "control": str(control_file),
            "dis_hru": "parameters_dis_hru.nc",
            "dis_both": "parameters_dis_both.nc",
            "solargeometry": {
                "class": "PRMSSolarGeometry",
                "parameters": "parameters_PRMSSolarGeometry.nc",
                "dis": "dis_hru",
            },
            "atmosphere": {
                "class": "PRMSAtmosphere",
                "parameters": "parameters_PRMSAtmosphere.nc",
                "dis": "dis_hru",
            },
            "canopy": {
                "class": "PRMSCanopy",
                "parameters": "parameters_PRMSCanopy.nc",
                "dis": "dis_hru",
            },
            "snow": {
                "class": "PRMSSnow",
                "parameters": "parameters_PRMSSnow.nc",
                "dis": "dis_hru",
            },
            "model_order": ["solargeometry", "atmosphere", "canopy", "snow"],
        }
        model_dict_file = domain_dir / "example_model_dict.yml"
        dump_dict = {control_file: control, model_dict_file: model_dict}
        for key, val in dump_dict.items():
            with open(key, "w") as file:
                documents = yaml.dump(val, file)

        model = pws.Model.from_yml(model_dict_file)
        model.run()
        control_file.unlink()
        model_dict_file.unlink()

    >>> import yaml
    >>> import pywatershed as pws
    >>> test_data_dir = pws.constants.__pywatershed_root__ / "../test_data"
    >>> domain_dir = test_data_dir / "drb_2yr"
    >>> control = {
    ...     "start_time": "1979-01-01T00:00:00",
    ...     "end_time": "1980-12-31T00:00:00",
    ...     "time_step": 24,
    ...     "time_step_units": "h",
    ...     "verbosity": 0,
    ...     "budget_type": "warn",
    ...     "init_vars_from_file": 0,
    ...     "input_dir": str(domain_dir),
    ... }
    >>> control_file = domain_dir / "example_control.yml"
    >>> model_dict = {
    ...     "control": str(control_file),
    ...     "dis_hru": "parameters_dis_hru.nc",
    ...     "dis_both": "parameters_dis_both.nc",
    ...     "solargeometry": {
    ...         "class": "PRMSSolarGeometry",
    ...         "parameters": "parameters_PRMSSolarGeometry.nc",
    ...         "dis": "dis_hru",
    ...     },
    ...     "atmosphere": {
    ...         "class": "PRMSAtmosphere",
    ...         "parameters": "parameters_PRMSAtmosphere.nc",
    ...         "dis": "dis_hru",
    ...     },
    ...     "canopy": {
    ...         "class": "PRMSCanopy",
    ...         "parameters": "parameters_PRMSCanopy.nc",
    ...         "dis": "dis_hru",
    ...     },
    ...     "snow": {
    ...         "class": "PRMSSnow",
    ...         "parameters": "parameters_PRMSSnow.nc",
    ...         "dis": "dis_hru",
    ...     },
    ...     "model_order": ["solargeometry", "atmosphere", "canopy", "snow"],
    ... }
    >>> model_dict_file = domain_dir / "example_model_dict.yml"
    >>> dump_dict = {control_file: control, model_dict_file: model_dict}
    >>> for key, val in dump_dict.items():
    ...     with open(key, "w") as file:
    ...         documents = yaml.dump(val, file)
    ...
    >>> model = pws.Model.from_yml(model_dict_file)
    PRMSCanopy jit compiling with numba
    PRMSSnow jit compiling with numba
    >>> model.run()
      0%|                                                                     | 0/731 [00:00<?, ?it/s]
    /Users/jmccreight/usgs/pywatershed2/pywatershed/hydrology/prms_snow.py:1086: NumbaExperimentalFeatureWarning: Use of isinstance() detected. This is an experimental feature.
      through_rain[:] = np.where(wh_through, net_rain, zero)
      0%|                                                           | 1/731 [00:05<1:08:27,  5.63s/it]
    /Users/jmccreight/usgs/pywatershed2/pywatershed/base/budget.py:317: UserWarning: The flux unit balance not equal to the change in unit storage: PRMSSnow
      warn(msg, UserWarning)
    100%|██████████████████████████████████████████████████████████| 731/731 [00:06<00:00, 119.95it/s]
    model.run(): finalizing
    >>> control_file.unlink()
    >>> model_dict_file.unlink()

    """

    def __init__(
        self,
        process_list_or_model_dict,
        control: Control = None,
        parameters: Union[Parameters, dict[Parameters]] = None,
        find_input_files: bool = True,
    ):
        self.control = control
        self.parameters = parameters

        # This is for backwards compatibility
        msg = "Inputs are inconsistent"
        if isinstance(process_list_or_model_dict, (list, tuple)):
            # take the old-school-style inputs and convert to new-school inputs
            # may be deprecated in the future.
            assert control is not None, msg
            assert isinstance(parameters, PrmsParameters), msg

            # eventually handle a file for parameters?
            # proc_param_file = domain["dir"] / f"parameters_{proc_name}.nc"
            # proc["parameters"] = PrmsParameters.from_netcdf(
            #     proc_param_file
            # )

            self.model_dict = {}
            model_dict = self.model_dict
            model_dict["control"] = self.control
            for process in process_list_or_model_dict:
                proc_name = process.__name__
                model_dict[proc_name] = {}
                model_dict[proc_name]["class"] = process
                model_dict[proc_name]["parameters"] = parameters

            model_dict["model_order"] = [
                key for key in process_order_nhm if key in model_dict.keys()
            ]

        elif isinstance(process_list_or_model_dict, dict):
            assert control is None, msg
            assert parameters is None, msg
            self.model_dict = process_list_or_model_dict

        else:
            raise ValueError("Invalid type of process_list_or_model_dict")

        self._categorize_model_dict()
        self._validate_model_dict()
        self._set_input_dir()
        self._solve_inputs()
        self._init_procs()
        self._connect_procs()

        self._found_input_files = False
        if find_input_files:
            self._find_input_files()

        # methodize this netcdf section
        self._netcdf_initialized = False
        if "netcdf_output_dir" in self.control.options.keys():
            if "netcdf_output_var_names" in self.control.options.keys():
                output_vars = self.control.options["netcdf_output_var_names"]
            else:
                output_vars = None
            if "netcdf_output_separate_files" in self.control.options.keys():
                separate_files = self.control.options[
                    "netcdf_output_separate_files"
                ]
            else:
                separate_files = True

            self.initialize_netcdf(
                self.control.options["netcdf_output_dir"],
                output_vars=output_vars,
                separate_files=separate_files,
            )

        return

    def _categorize_model_dict(self):
        """Categorize model_dict entries

        Categorize entries into [control, dis, order, process, exchange]
        categories based on types and attributes
        """
        # TODO: add exchanges
        category_type_dict = {
            "dis": Parameters,
            "control": Control,
            "order": list,
            "process": dict,
            # "exchange": exchange class?
        }
        category_key_dict = {}
        for cat, typ in category_type_dict.items():
            category_key_dict[cat] = [
                kk for kk, vv in self.model_dict.items() if isinstance(vv, typ)
            ]

        self._category_key_dict = category_key_dict

        return

    def _validate_model_dict(self):
        cat_key_dict = self._category_key_dict
        # * a single control object
        assert len(cat_key_dict["control"]) == 1
        cat_key_dict["control"] = cat_key_dict["control"][0]
        self._control_name = cat_key_dict["control"]
        if self.control is None:
            self.control = self.model_dict[self._control_name]
        # * a single order list
        assert len(cat_key_dict["order"]) == 1
        cat_key_dict["order"] = cat_key_dict["order"][0]
        self._order_name = cat_key_dict["order"]
        self.process_order = self.model_dict[self._order_name]
        # * all processes are in order and vice-versa
        assert set(cat_key_dict["process"]) == set(
            self.model_dict[cat_key_dict["order"]]
        )
        # * all required proc dis are present
        # * all required proc in exchanges are present

        return

    def _solve_inputs(self):
        """What processes supply inputs to others, what files are needed?

        TODO: this does not currently take order into account, which could
              be important and is now possible with order specified.
        """
        proc_dict = {
            proc_name: self.model_dict[proc_name]["class"]
            for proc_name in self._category_key_dict["process"]
        }

        proc_inputs = {kk: vv.get_inputs() for kk, vv in proc_dict.items()}
        proc_vars = {kk: vv.get_variables() for kk, vv in proc_dict.items()}

        # Solve where inputs come from
        inputs_from = {}
        for comp in proc_dict.keys():
            inputs = deepcopy(proc_inputs)
            vars = deepcopy(proc_vars)
            c_inputs = inputs.pop(comp)
            _ = vars.pop(comp)

            inputs_from[comp] = {}
            for input in c_inputs:
                inputs_ptr = inputs_from
                inputs_ptr[comp][input] = []  # could use None
                for other in inputs.keys():
                    if input in vars[other]:
                        inputs_ptr[comp][input] += [other]
                        # this should be a list of length one
                        # check?

        # If inputs dont come from other processes, assume they come from
        # file in input_dir. Exception is that PRMSAtmosphere requires its
        # files on init, so dont adapt these
        file_input_names = set([])
        for k0, v0 in inputs_from.items():
            for k1, v1 in v0.items():
                if not v1:
                    file_input_names = file_input_names.union([k1])

        # initiate the file inputs here rather than in the processes
        file_inputs = {}
        # Use dummy names for now
        for name in file_input_names:
            file_inputs[name] = pl.Path(name)

        self._proc_dict = proc_dict
        self._inputs_from = inputs_from
        self._file_inputs = file_inputs
        self._file_input_names = file_input_names
        return

    def _init_procs(self):
        # instantiate processes: instance dict
        self.processes = {}
        for proc_name in self.model_dict[self._order_name]:
            # establish the args for init
            proc_specs = self.model_dict[proc_name]

            dis = None
            if "dis" in proc_specs:
                dis = self.model_dict[proc_specs["dis"]]

            process_inputs = {
                input: None
                for input in self._proc_dict[proc_name].get_inputs()
            }

            args = {
                "control": self.control,
                "discretization": dis,
                "parameters": proc_specs["parameters"],
                **process_inputs,
            }

            self.processes[proc_name] = self._proc_dict[proc_name](**args)

        # <
        return

    def _connect_procs(self):
        self.process_input_from = {}
        for process in self.process_order:
            self.process_input_from[process] = {}
            for input, frm in self._inputs_from[process].items():
                if not frm:
                    self.process_input_from[process][
                        input
                    ] = self._file_inputs[input]
                else:
                    self.process_input_from[process][input] = frm[0]
                    self.processes[process].set_input_to_adapter(
                        input,
                        adapter_factory(
                            self.processes[frm[0]][input],
                            control=self.control,
                        ),  # drop list above
                    )
        #   <   <   <
        return

    def _set_input_dir(self):
        if "input_dir" not in self.control.options.keys():
            msg = "Required control option 'input_dir' not found"
            raise ValueError(msg)
        else:
            self._input_dir = pl.Path(
                self.control.options["input_dir"]
            ).resolve()

        return

    def _find_input_files(self) -> None:
        file_inputs = {}
        for name in self._file_input_names:
            nc_path = self._input_dir / f"{name}.nc"
            file_inputs[name] = adapter_factory(
                nc_path,
                name,
                control=self.control,
            )
        for process in self.process_order:
            for input, frm in self._inputs_from[process].items():
                if not frm:
                    fname = file_inputs[input]._fname
                    self.process_input_from[process][input] = fname
                    self.processes[process].set_input_to_adapter(
                        input, file_inputs[input]
                    )

        self._found_input_files = True
        return

    @staticmethod
    def model_dict_from_yml(yml_file: Union[str, pl.Path]) -> dict:
        """Generate a model dictionary from a yaml file.

        Instead of Model.from_yml() it can be useful to get the model
        dictionary before passing it to Model.

        Args:
            yml_file: a yml file

        Returns:
            A model dictionary.
        """
        import yaml

        import pywatershed

        with pl.Path(yml_file).open("r") as file_stream:
            model_dict = yaml.load(file_stream, Loader=yaml.Loader)

        for key, val in model_dict.items():
            if isinstance(val, str):
                val_pl = path_rel_to_yml(val, yml_file)
                if val.endswith(".yml"):
                    model_dict[key] = Control.from_yml(val_pl)
                elif val.endswith(".nc"):
                    model_dict[key] = Parameters.from_netcdf(val_pl)
                else:
                    msg = (
                        "Unsupported file extension for control (.yml)"
                        "and parameter (.nc) file paths in model yaml file"
                    )
                    raise ValueError(msg)

            elif isinstance(val, dict):
                # Until exchanges are introduced processes are the only
                # accepted dictionaries
                if "class" in val.keys():
                    # processes
                    cls = val["class"]
                    val["class"] = getattr(pywatershed, cls)
                    par = val["parameters"]
                    par_pl = path_rel_to_yml(par, yml_file)
                    val["parameters"] = Parameters.from_netcdf(
                        par_pl, encoding=False
                    )
                    # dis = val["dis"]
                    # val["dis"] = model_dict[dis]

            elif isinstance(val, list):
                pass

        return model_dict

    @staticmethod
    def from_yml(yml_file: Union[str, pl.Path]):
        """Instantiate a Model from a yaml file

        A yaml file that specifies a model_dict as the first argument of Model.

        Args:
           yml_file: str or pathlib.Path

        Returns:
           An instance of Model.

        Yaml file structure (strict order not required, but suggested):

        Control object: Any name can be used but the value must be a control
            yaml file specified with the suffix ".yml". E.g "name: control.yml"
            would appear in the passed yaml file. Only one control
            specification is allowed in the yml_file. For details on the
            requirements of the control.yml file see `Control.from_yml`
        Discretization objects: Any number of discretization objects can be
            supplied with arbitrary (though unique) names. The values supplied
            for each discretization must be a valid netcdf file with suffix
            ".nc". These can generally be obtained by calling
            `parameters.to_netcdf()` on subclasses of Parameters.
        Process objects: Any number of processes may be specified with
            arbitrary (though unique) names. Processes are specified as
            dictionaries in the yaml file, they have additonal key:value pairs.
            Required key:value pairs are:
                class: a class that can be `from pywatershed import class`
                parameters: a netcdf file that specifies the parameters for
                    the class
                dis: the name of the discretization for the class as given by
                    a discretization object speficication above.
            Optional key:value pairs:
                TBD
        Model order list: a list supplying the order in which the processes are
            to be executed.

        Note: To get a model_dict specfied by the yml_file, call
        `model_dict_from_yml` instead.

        """
        return Model(Model.model_dict_from_yml(yml_file))

    def initialize_netcdf(
        self,
        output_dir: str,
        separate_files: bool = True,
        budget_args: dict = None,
        output_vars: list = None,
    ):
        """Initialize NetCDF output files for model (all processes)."""

        if self._netcdf_initialized:
            msg = (
                "Model class previously initialized netcdf output "
                f"in {self._netcdf_dir}"
            )
            warn(msg)
            return

        if not self._found_input_files:
            self._find_input_files()

        self._netcdf_dir = pl.Path(output_dir)
        for cls in self.process_order:
            self.processes[cls].initialize_netcdf(
                output_dir=self._netcdf_dir,
                separate_files=separate_files,
                budget_args=budget_args,
                output_vars=output_vars,
            )
        self._netcdf_initialized = True
        return

    def run(
        self,
        netcdf_dir: fileish = None,
        finalize: bool = True,
        n_time_steps: int = None,
        output_vars: list = None,
    ):
        """Run the model.

        Running the model wraps:
            * initializing netcdf output (optional on netcdf_dir argument)
            * full time loop (from control object) with progress reports
                * advance
                * calculate
                * output (optional on output methods requested)
            * finalize (optional)

        Args:
            netcdf_dir: optional directory to netcdf output files (initializes
               netcdf and outputs at each timestep).
            finalize: option to not finalize at the end of the time loop.
               Default is to finalize.
            n_time_steps: the number of timesteps to run
            output_vars: the vars to output to the netcdf_dir
        """
        if not self._found_input_files:
            self._find_input_files()

        if netcdf_dir:
            print("model.run(): initializing NetCDF output")
            self.initialize_netcdf(netcdf_dir, output_vars=output_vars)

        if not n_time_steps:
            n_time_steps = self.control.n_times

        for istep in tqdm(range(n_time_steps)):
            self.advance()
            self.calculate()
            self.output()

        if finalize:
            print("model.run(): finalizing")
            self.finalize()

        return

    def advance(self):
        """Advance the model in time."""
        if not self._found_input_files:
            self._find_input_files()

        self.control.advance()
        for cls in self.process_order:
            self.processes[cls].advance()
        return

    def calculate(self):
        """Calculate the model."""
        for cls in self.process_order:
            self.processes[cls].calculate(1.0)
        return

    def output(self):
        """Output the model at the current time."""
        for cls in self.process_order:
            self.processes[cls].output()
        return

    def finalize(self):
        """Finalize the model."""
        for cls in self.process_order:
            self.processes[cls].finalize()
        return
