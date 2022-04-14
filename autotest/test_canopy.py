import pathlib
from datetime import datetime

import numpy as np
import pytest

from pynhm.atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from pynhm.base.control import Control
from pynhm.hydrology.PRMSCanopy import PRMSCanopy
from pynhm.preprocess import CsvFile
from pynhm.utils import ControlVariables
from pynhm.utils.parameters import PrmsParameters

forcings_dict = {
    "datetime": np.array(
        [
            datetime(1979, 1, 3, 0, 0),
            datetime(1979, 1, 4, 0, 0),
        ]
    ),
    "spatial_id": np.array([5307, 5308]),
    "tmin": np.array(
        [
            [46.1209, 45.76805],
            [37.7609, 37.4881],
        ]
    ),
    "rhavg": np.array(
        [
            [82.45999908447266, 82.5999984741211],
            [81.98999786376953, 82.3499984741211],
        ]
    ),
    "tmax": np.array(
        [
            [57.41188049316406, 56.47270965576172],
            [55.511878967285156, 55.032711029052734],
        ]
    ),
    "snowfall": np.array([[0.0, 0.0], [0.0, 0.0]]),
    "prcp": np.array(
        [
            [0.31392958760261536, 0.24780480563640594],
            [0.6605601906776428, 0.5214226245880127],
        ]
    ),
    "rainfall": np.array(
        [
            [0.31392958760261536, 0.24780480563640594],
            [0.6605601906776428, 0.5214226245880127],
        ]
    ),
    "potet": np.array(
        [
            [0.25, 0.26],
            [0.26, 0.27],
        ]
    ),
}


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


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
        for key in PRMSCanopy.get_input_variables():
            input_variables[key] = np.ones([nhru])

        # todo: this is testing instantiation, but not physics
        ntimes = control.n_times
        self.cnp = PRMSCanopy(
            control=control, params=prms_params, **input_variables
        )
        self.cnp.advance(itime_step=0)
        self.cnp.calculate(time_length=1.0)

        return


class TestPRMSCanopyDomain:
    def test_init(self, domain):
        prms_params = PrmsParameters.load(domain["param_file"])

        # Set information from the control file
        control_file = domain["control_file"]
        control = ControlVariables.load(control_file)
        start_time = control.control.start_time
        end_time = control.control.end_time
        initial_deltat = control.control.initial_deltat

        atm_information_dict = {
            "start_time": start_time,
            "end_time": end_time,
            "time_step": initial_deltat,
            "verbosity": 3,
            "height_m": 5,
        }
        var_translate = {
            "prcp_adj": "prcp",
            "rainfall_adj": "rainfall",
            "snowfall_adj": "snowfall",
            "tmax_adj": "tmax",
            "tmin_adj": "tmin",
            "swrad": "swrad",
            "potet": "potet",
        }
        var_file_dict = {
            var_translate[var]: file
            for var, file in domain["prms_outputs"].items()
            if var in var_translate.keys()
        }
        atm = NHMBoundaryLayer.load_prms_output(
            **var_file_dict, parameters=prms_params, **atm_information_dict
        )
        atm.calculate_sw_rad_degree_day()
        atm.calculate_potential_et_jh()

        ## NEWNEWNEW
        output_files = domain["prms_outputs"]
        input_variables = {}
        for key in PRMSCanopy.get_input_variables():
            nc_pth = output_files[key].with_suffix(".nc")
            input_variables[key] = nc_pth

        # pkwater_equiv comes from snowpack; it is lagged by a time step
        prms_output_files = domain["prms_outputs"]
        pkwater_equiv = None
        if "pkwater_equiv" in prms_output_files:
            fname = prms_output_files["pkwater_equiv"]
            df = CsvFile(fname).to_dataframe()
            pkwater_equiv = df.to_numpy()

        transp_on = None
        if "transp_on" in prms_output_files:
            fname = prms_output_files["transp_on"]
            df = CsvFile(fname).to_dataframe()
            transp_on = df.to_numpy()

        self.cnp = PRMSCanopy(prms_params, atm, pkwater_equiv, transp_on)
        self.cnp.advance(itime_step=0)

        for istep in range(atm.n_time):
            if istep > 0:
                atm.advance()
            self.cnp.advance(istep)
            self.cnp.calculate(1.0)

        # build a dictionary of the prms output dataframes for comparison
        prms_output_files = domain["prms_outputs"]
        comparison_variables = [
            "rainfall_adj",
            "snowfall_adj",
            "potet",
            "intcpstor",
            "net_rain",
            "net_snow",
            "intcp_evap",
        ]
        prms_output_dataframes = {}
        for cv in comparison_variables:
            fname = prms_output_files[cv]
            print(f"loading {fname}")
            csvobj = CsvFile(fname)
            df = csvobj.to_dataframe()
            prms_output_dataframes[cv] = df

        # get a dictionary of dataframes for process model output
        pynhm_output_dataframes = self.cnp.get_output_dataframes()

        # compare prms and pynhm data
        for cv in comparison_variables:
            prms_data = prms_output_dataframes[cv]
            pynhm_data = pynhm_output_dataframes[cv]

            print(f"\n{50*'*'}")
            print(f"{cv}  min  max")
            a1 = prms_data.to_numpy()
            a2 = pynhm_data.to_numpy()
            diff = a1 - a2
            diffmin = diff.min()
            diffmax = diff.max()
            print(f"prms   {a1.min()}    {a1.max()}")
            print(f"pynhm  {a2.min()}    {a2.max()}")
            print(f"diff   {diffmin}  {diffmax}")

            atol = 1.0e-5
            errmsg = f"Canopy variable {cv} does not match to within {atol}"
            assert np.allclose(diffmin, 0.0, atol=atol), errmsg
            assert np.allclose(diffmax, 0.0, atol=atol), errmsg

        # save cnp output as dataframes in temp/domain_name
        saveoutput = False
        if saveoutput:
            pth = pathlib.Path(".", "temp", domain["domain_name"])
            pth.mkdir(parents=True, exist_ok=True)
            self.cnp.output_to_csv(pth)

        return
