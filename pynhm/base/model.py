import math
import pathlib as pl
from copy import deepcopy

from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.constants import fileish


class Model:
    """pynhmn model builder class.

    Build a model from process classes. Model wires the inputs and outputs,
    searching for unavailable inputs from file (input_dir).

    Args:
        process_classes: N process classes composing which compose the model.
            Currently, unknown classes will break the model builder with
            an informative error.
        control: Control object.
        input_dir: A directory to search for input files.
        budget_type: None, "warn", or "error".
        verbose: Boolean.
        find_input_files: Search/find input file on __init__ or delay until run
           or advance of the model. Delaying (False) allows ModelGraph of the
           specified model without the need for input files.
    """

    def __init__(
        self,
        *process_classes,
        control: Control,
        input_dir: str = None,
        budget_type: str = "error",  # todo: also pass dict
        netcdf_output_dir: str = None,
        verbose: bool = False,
        find_input_files: bool = True,
    ):

        self.control = control
        self.input_dir = input_dir
        self.verbose = verbose

        class_dict = {comp.__name__: comp for comp in process_classes}
        class_inputs = {kk: vv.get_inputs() for kk, vv in class_dict.items()}
        class_vars = {kk: vv.get_variables() for kk, vv in class_dict.items()}

        # Solve where inputs come from
        self._inputs_from = {}
        # inputs_from_prev = {}
        for comp in class_dict.keys():
            inputs = deepcopy(class_inputs)
            vars = deepcopy(class_vars)
            c_inputs = inputs.pop(comp)
            _ = vars.pop(comp)
            self._inputs_from[comp] = {}
            # inputs_from_prev[comp] = {}
            # inputs
            for input in c_inputs:
                inputs_ptr = self._inputs_from
                # if re.compile(".*_(ante|prev|old)").match(input):
                #    inputs_ptr = inputs_from_prev
                inputs_ptr[comp][input] = []  # could use None
                for other in inputs.keys():
                    if input in vars[other]:
                        inputs_ptr[comp][input] += [other]
                        # this should be a list of length one
                        # check?

        # determine process order
        # for now this is sufficient. when more classes exist, we can analyze
        # dependencies, might inputs_from_prev above.
        process_order_prms = [
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
        self.process_order = []
        for comp in process_order_prms:
            if comp in class_dict.keys():
                self.process_order += [comp]

        missing_processes = set(self.process_order).difference(
            set(class_dict.keys())
        )
        if missing_processes:
            raise ValueError(
                f"Process {missing_processes} not currently "
                f"handled by Model class."
            )

        # If inputs dont come from other processes, assume they come from
        # file in input_dir. Exception is that PRMSAtmosphere requires its
        # files on init, so dont adapt these
        self.file_input_names = set([])
        for k0, v0 in self._inputs_from.items():
            for k1, v1 in v0.items():
                if not v1:
                    self.file_input_names = self.file_input_names.union([k1])

        # initiate the file inputs here rather than in the processes
        file_inputs = {}
        # Use dummy names for now
        for name in self.file_input_names:
            file_inputs[name] = pl.Path(name)

        # instantiate processes: instance dict
        self.processes = {}
        for process in self.process_order:
            process_inputs = {
                input: None for input in class_dict[process].get_inputs()
            }

            # This is a hack. Need to make StorageUnit a subclass of a more
            # general model/element/process class
            args = {
                "control": control,
                **process_inputs,
            }
            if process not in ["PRMSSolarGeometry"]:
                args["budget_type"] = budget_type
            self.processes[process] = class_dict[process](**args)

        # Wire it up
        self.process_input_from = {}
        for process in self.process_order:
            self.process_input_from[process] = {}
            for input, frm in self._inputs_from[process].items():
                if not frm:
                    self.process_input_from[process][input] = file_inputs[
                        input
                    ]
                else:
                    self.process_input_from[process][input] = frm[0]
                    self.processes[process].set_input_to_adapter(
                        input,
                        adapter_factory(
                            self.processes[frm[0]][input], control=control
                        ),  # drop list above
                    )

        self._found_input_files = False
        if find_input_files:
            self._find_input_files()

        if netcdf_output_dir:
            self.initialize_netcdf(output_dir=netcdf_output_dir)

        return

    def _find_input_files(self):
        file_inputs = {}
        for name in self.file_input_names:
            nc_path = self.input_dir / f"{name}.nc"
            file_inputs[name] = adapter_factory(
                nc_path, name, control=self.control
            )
        for process in self.process_order:
            self.process_input_from[process] = {}
            for input, frm in self._inputs_from[process].items():
                if not frm:
                    fname = file_inputs[input]._fname
                    self.process_input_from[process][input] = fname
                    self.processes[process].set_input_to_adapter(
                        input, file_inputs[input]
                    )

        self._found_input_files = True
        return

    def initialize_netcdf(self, **kwargs):
        """Initialize NetCDF output files for model (all processes).

        The kwargs are meant to be passed to storageUnit.initialize_netcdf

        kwrgs:
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
        if not self._found_input_files:
            self._find_input_files()

        self._netcdf_dir = kwargs["output_dir"]
        for cls in self.process_order:
            self.processes[cls].initialize_netcdf(**kwargs)
        return

    def run(self, netcdf_dir: fileish = None, finalize: bool = True):
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
        """
        if not self._found_input_files:
            self._find_input_files()

        if netcdf_dir:
            print("model.run(): initializing NetCDF output")
            self.initialize_netcdf(output_dir=netcdf_dir)

        last_pct_comp = 0
        print(f"model.run(): {last_pct_comp} % complete")

        n_time_steps = self.control.n_times
        for istep in range(n_time_steps):

            # progress for the impatient
            pct_complete = math.floor((istep + 1) / n_time_steps * 100)
            if not (pct_complete % 10) and pct_complete != last_pct_comp:
                last_pct_comp = pct_complete
                print(f"model.run(): {pct_complete} % complete")

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
