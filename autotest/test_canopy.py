import pathlib as pl
from datetime import datetime

import numpy as np
import pytest

from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.hydrology.PRMSCanopy import PRMSCanopy
from pynhm.preprocess import CsvFile
from pynhm.utils.netcdf_utils import NetCdfCompare
from pynhm.utils.parameters import PrmsParameters


class TestPRMSCanopySimple:
    def test_init(self):

        time_dict = {
            "start_time": np.datetime64("1979-01-03T00:00:00.00"),
            "end_time": np.datetime64("1979-01-04T00:00:00.00"),
            "time_step": np.timedelta64(1, "D"),
        }
        control = Control(**time_dict)

        nhru = 2
        prms_params = {
            "nhru": nhru,
            "hru_area": np.array(nhru * [1.0]),
            "covden_sum": np.array(nhru * [0.5]),
            "covden_win": np.array(nhru * [0.5]),
            "srain_intcp": np.array(nhru * [1.0]),
            "wrain_intcp": np.array(nhru * [1.0]),
            "snow_intcp": np.array(nhru * [1.0]),
            "epan_coef": np.array(nhru * [1.0]),
            "potet_sublim": np.array(nhru * [1.0]),
            "cov_type": np.array(nhru * [1]),
        }
        prms_params = PrmsParameters(prms_params)

        input_variables = {}
        for key in PRMSCanopy.get_inputs():
            input_variables[key] = np.ones([nhru])

        # todo: this is testing instantiation, but not physics
        cnp = PRMSCanopy(
            control=control, params=prms_params, **input_variables
        )
        control.advance()
        cnp.advance()
        cnp.calculate(time_length=1.0)

        return


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


class TestPRMSCanopyDomain:
    def test_init(self, domain, control, params, tmp_path):
        tmp_path = pl.Path(tmp_path)
        output_dir = domain["prms_output_dir"]

        # get the answer data
        comparison_var_names = [
            "net_rain",
            "net_snow",
            "net_ppt",
            "intcp_stor",
            "intcp_evap",
            "hru_intcpstor",
            "hru_intcpevap",
        ]
        ans = {}
        for key in comparison_var_names:
            nc_pth = output_dir / f"{key}.nc"
            ans[key] = adapter_factory(nc_pth, variable_name=key)

        # setup the canopy
        input_variables = {}
        for key in PRMSCanopy.get_inputs():
            nc_pth = output_dir / f"{key}.nc"
            input_variables[key] = nc_pth

        cnp = PRMSCanopy(control=control, params=params, **input_variables)

        all_success = True
        for istep in range(control.n_times):
            control.advance()
            cnp.advance()
            cnp.calculate(1.0)

            # compare along the way
            atol = 1.0e-5
            for key, val in ans.items():
                val.advance()
            for key in ans.keys():
                a1 = ans[key].current
                a2 = cnp[key]
                success = np.isclose(a1, a2, atol=atol).all()
                if not success:
                    all_success = False
                    diff = a1 - a2
                    diffmin = diff.min()
                    diffmax = diff.max()
                    print(f"time step {istep}")
                    print(f"output variable {key}")
                    print(f"prms   {a1.min()}    {a1.max()}")
                    print(f"pynhm  {a2.min()}    {a2.max()}")
                    print(f"diff   {diffmin}  {diffmax}")

        cnp.finalize()

        if not all_success:
            raise Exception("pynhm results do not match prms results")
        # is comparing along the way slower or faster than comparing netcdf?

        # prms_output_dataframes = {}
        # for cv in comparison_variables:
        #     fname = prms_output_files[cv]
        #     print(f"loading {fname}")
        #     csvobj = CsvFile(fname)
        #     df = csvobj.to_dataframe()
        #     prms_output_dataframes[cv] = df

        # # get a dictionary of dataframes for process model output
        # pynhm_output_dataframes = cnp.get_output_dataframes()

        # # compare prms and pynhm data
        # for cv in comparison_variables:
        #     prms_data = prms_output_dataframes[cv]
        #     pynhm_data = pynhm_output_dataframes[cv]

        #     print(f"\n{50*'*'}")
        #     print(f"{cv}  min  max")
        #     a1 = prms_data.to_numpy()
        #     a2 = pynhm_data.to_numpy()
        #     diff = a1 - a2
        #     diffmin = diff.min()
        #     diffmax = diff.max()
        #     print(f"prms   {a1.min()}    {a1.max()}")
        #     print(f"pynhm  {a2.min()}    {a2.max()}")
        #     print(f"diff   {diffmin}  {diffmax}")

        #     atol = 1.0e-5
        #     errmsg = f"Canopy variable {cv} does not match to within {atol}"
        #     assert np.allclose(diffmin, 0.0, atol=atol), errmsg
        #     assert np.allclose(diffmax, 0.0, atol=atol), errmsg

        # # save cnp output as dataframes in temp/domain_name
        # saveoutput = False
        # if saveoutput:
        #     pth = pathlib.Path(".", "temp", domain["domain_name"])
        #     pth.mkdir(parents=True, exist_ok=True)
        #     cnp.output_to_csv(pth)

        return
