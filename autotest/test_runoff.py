import pathlib as pl
from datetime import datetime

import numpy as np
import pytest

from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.hydrology.PRMSRunoff import PRMSRunoff
from pynhm.preprocess import CsvFile
from pynhm.utils import ControlVariables
from pynhm.utils.netcdf_utils import NetCdfCompare
from pynhm.utils.parameters import PrmsParameters


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


class TestPRMSRunoffDomain:
    def test_init(self, domain, control, params, tmp_path):
        tmp_path = pl.Path(tmp_path)

        # get the answer data

        comparison_var_names = [
            "contrib_fraction",
        ]
        output_dir = domain["prms_output_dir"]

        # Read PRMS output into ans for comparison with pynhm results
        ans = {}
        for key in comparison_var_names:
            nc_pth = output_dir / f"{key}.nc"
            ans[key] = adapter_factory(nc_pth, variable_name=key)

        # instantiate runoff
        input_variables = {}
        for key in PRMSRunoff.get_inputs():
            nc_pth = output_dir / f"{key}.nc"
            if "soil_moist" in str(nc_pth):
                nc_pth = output_dir / "soil_moist_prev.nc"
            input_variables[key] = nc_pth

        runoff = PRMSRunoff(control=control, params=params, **input_variables)

        all_success = True
        for istep in range(control.n_times):
            # control.advance()
            runoff.advance()
            runoff.calculate(1.0)

            # advance the answer, which is being read from a netcdf file
            for key, val in ans.items():
                val.advance()

            # make a comparison check with answer
            check = False
            if check:
                atol = 1.0e-5
                success = self.check_timestep_results(runoff, istep, ans, atol)
                if not success:
                    all_success = False

        runoff.finalize()

        # check at the end and error if one or more steps didn't pass
        if not all_success:
            raise Exception("pynhm results do not match prms results")

        return

    @staticmethod
    def check_timestep_results(storageunit, istep, ans, atol):
        for key in ans.keys():
            a1 = ans[key].current
            a2 = storageunit[key]
            success = np.isclose(a1, a2, atol=atol).all()
            if not success:
                diff = a1 - a2
                diffmin = diff.min()
                diffmax = diff.max()
                if True:
                    print(f"time step {istep}")
                    print(f"output variable {key}")
                    print(f"prms   {a1.min()}    {a1.max()}")
                    print(f"pynhm  {a2.min()}    {a2.max()}")
                    print(f"diff   {diffmin}  {diffmax}")
        return success
