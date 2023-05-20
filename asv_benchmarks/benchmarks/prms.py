import pathlib as pl
import shutil

import pynhm as pws

from . import parameterized, test_data_dir

domains = ["hru_1", "drb_2yr", "ucb_2yr"]


model_tests = {
    # "solar": (pws.PRMSSolarGeometry,),
    # "atm": (pws.PRMSAtmosphere),
    "canopy": (pws.PRMSCanopy,),
    "solar-atm": (pws.PRMSSolarGeometry, pws.PRMSAtmosphere),
}
model_tests_inv = {v: k for k, v in model_tests.items()}


class PRMSBasics:
    """Benchmark simple components object of PRMS models"""

    @parameterized(
        ["domain"],
        (domains),
    )
    def time_parameter_read_prms(self, domain):
        parameter_file = test_data_dir / f"{domain}/myparam.param"
        _ = pws.PrmsParameters.load(parameter_file)
        return

    @parameterized(
        ["domain"],
        (domains),
    )
    def time_control_read(self, domain):
        control_file = test_data_dir / f"{domain}/control.test"
        _ = pws.Control.load(control_file)
        return


class PRMSModels:
    """Benchmark simple components object of PRMS models"""

    def setup(self, *args):
        self.domain = args[0]
        self.processes = args[1]
        self.tag = model_tests_inv[self.processes]

        self.control_file = test_data_dir / f"{self.domain}/control.test"
        self.parameter_file = test_data_dir / f"{self.domain}/myparam.param"

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

    # Helper function
    def model_setup_run(
        self,
        domain: str = None,
        processes: tuple = None,
        write_output: bool = None,
    ):
        print(f"model_setup_run tag: {self.tag}")
        params = pws.PrmsParameters.load(self.parameter_file)
        control = pws.Control.load(self.control_file, params=params)

        model = pws.Model(
            *self.processes,
            control=control,
            input_dir=self.tag_input_dir,
            budget_type="warn",
            calc_method="numba",
        )
        if write_output:
            model.initialize_netcdf(self.tag_dir)
        model.run(finalize=True)

    @parameterized(
        ["domain", "processes"],
        (
            domains,
            list(model_tests.values()),
        ),
    )
    def time_model_setup_run_no_output(
        self,
        domain: str,
        processes: tuple,
    ):
        _ = self.model_setup_run(
            domain=domain, processes=processes, write_output=False
        )
