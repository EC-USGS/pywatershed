import pathlib as pl
import shutil
from typing import Union, Literal

from . import _is_pws, parameterized, test_data_dir

if _is_pws:
    import pywatershed as pws
else:
    import pynhm as pws

# TODO remove backwards compatiability with pynhm once
#      reported slows downs are sorted output
# TODO remove backwards compatability with <0.2.0 once
#      it is released.

domains = ["hru_1", "drb_2yr", "ucb_2yr"]
outputs = [None, "separate", "together"]
n_time_steps = 183

model_tests = {
    "solar": (pws.PRMSSolarGeometry,),
    "atm": (pws.PRMSAtmosphere,),
    "canopy": (pws.PRMSCanopy,),
    "snow": (pws.PRMSSnow,),
    "runoff": (pws.PRMSRunoff,),
    "soil": (pws.PRMSSoilzone,),
    "gw": (pws.PRMSGroundwater,),
    "channel": (pws.PRMSChannel,),
    "nhm": (
        pws.PRMSSolarGeometry,
        pws.PRMSAtmosphere,
        pws.PRMSCanopy,
        pws.PRMSSnow,
        pws.PRMSRunoff,
        pws.PRMSSoilzone,
        pws.PRMSGroundwater,
        pws.PRMSChannel,
    ),
}
model_tests_inv = {v: k for k, v in model_tests.items()}


class PRMSBasics:
    """Benchmark simple components object of PRMS models"""

    @parameterized(
        ["domain"],
        (domains),
    )
    def time_prms_parameter_read(self, domain):
        parameter_file = test_data_dir / f"{domain}/myparam.param"
        if _is_pws:
            _ = pws.parameters.PrmsParameters.load(parameter_file)
        else:
            _ = pws.PrmsParameters.load(parameter_file)
        return

    @parameterized(
        ["domain"],
        (domains),
    )
    def time_prms_control_read(self, domain):
        control_file = test_data_dir / f"{domain}/control.test"
        if pws.__version__ == "0.2.0":
            _ = pws.Control.load(control_file)
        else:
            _ = pws.Control.load_prms(control_file, warn_unused_options=False)

        return


class PRMSModels:
    """Benchmark simple components object of PRMS models"""

    def setup(self, *args):
        self.domain = args[0]
        self.tag = args[1]
        self.processes = model_tests[self.tag]

        self.control_file = test_data_dir / f"{self.domain}/control.test"
        self.parameter_file = test_data_dir / f"{self.domain}/myparam.param"

        # backwards compatability pre pywatershed
        if _is_pws:
            self.params = pws.parameters.PrmsParameters.load(
                self.parameter_file
            )
        else:
            self.params = pws.PrmsParameters.load(self.parameter_file)

        # backwards compatability pre 0.2.0
        try:
            self.lt_v0_2_0 = True
            self.control = pws.Control.load(
                self.control_file, params=self.params
            )
        except:
            self.lt_v0_2_0 = False

        if hasattr(self, "control"):
            del self.control

        # setup input_dir with symlinked prms inputs and outputs
        self.domain_dir = pl.Path(f"PRMSModels_{self.domain}")
        self.domain_dir.mkdir(exist_ok=True)
        self.tag_dir = self.domain_dir / self.tag
        self.tag_input_dir = self.tag_dir / "input"
        self.tag_input_dir.mkdir(parents=True)
        for ff in (test_data_dir / self.domain).glob("*.nc"):
            shutil.copy(ff, self.tag_input_dir / ff.name)
            for ff in (test_data_dir / f"{self.domain}/output").glob("*.nc"):
                shutil.copy(ff, self.tag_input_dir / ff.name)

    def teardown(self, *args):
        shutil.rmtree(self.domain_dir)
        pass

    # Helper function
    def model_setup_run(
        self,
        domain: str = None,
        processes: tuple = None,
        write_output: Union[bool, Literal["separate", "together"]] = None,
    ):
        # seem to need to load control inside the model setup run bc
        # results are strange/inconsistent

        if not self.lt_v0_2_0:
            if pws.__version__ == "0.2.0":
                self.control = pws.Control.load(self.control_file)
            else:
                self.control = pws.Control.load_prms(
                    self.control_file, warn_unused_options=False
                )

            for oo in ["netcdf_output_dir", "netcdf_output_var_names"]:
                if oo in self.control.options.keys():
                    del self.control.options[oo]

            self.control.options["input_dir"] = self.tag_input_dir
            self.control.options["budget_type"] = "warn"
            self.control.options["calc_method"] = "numba"
            self.control.edit_n_time_steps(n_time_steps)

            model = pws.Model(
                self.processes,
                control=self.control,
                parameters=self.params,
            )

        else:
            self.control = pws.Control.load(
                self.control_file, params=self.params
            )
            self.control.edit_n_time_steps(n_time_steps)

            model = pws.Model(
                *self.processes,
                control=self.control,
                input_dir=self.tag_input_dir,
                budget_type="warn",
                calc_method="numba",
            )

        if write_output is not None:
            model.initialize_netcdf(
                output_dir=self.tag_dir,
                separate_files=(write_output == "separate"),
            )
        model.run(finalize=True)
        del model
        del self.control

    @parameterized(
        ["domain", "procs", "output"],
        (
            domains,
            list(model_tests.keys()),
            outputs,
        ),
    )
    def time_prms_run(
        self,
        domain: str,
        procs: str,
        output: Union[None, Literal["separate", "together"]],
    ):
        print(
            "\nPRMSModels args: \n",
            f"domain: {domain}\n",
            f"procs:{procs}\n",
            f"output: {output}\n",
        )

        _ = self.model_setup_run(
            domain=domain, processes=model_tests[procs], write_output=output
        )
