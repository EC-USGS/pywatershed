
import sys
sys.path.append("..")
import numpy as np
from datetime import datetime, timedelta


from pynhm.utils import ControlVariables
from pynhm.utils.parameters import PrmsParameters
from pynhm.atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from pynhm.canopy.PRMSCanopy import PRMSCanopy

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


class TestPRMSCanopySimple():

    def test_init(self):

        nhru = 2
        prms_params = {
            "nhru": nhru,
            "hru_area": np.array(nhru * [1.0]),
            "covden_sum": np.array(nhru * [0.5]),
            "covden_win": np.array(nhru * [0.5]),
            "srain_intcp": np.array(nhru * [1.0]),
            "wrain_intcp": np.array(nhru * [1.0]),
            "snow_intcp": np.array(nhru * [1.0]),
        }
        prms_params = PrmsParameters(prms_params)
        atm = NHMBoundaryLayer(
            forcings_dict,
            start_time=np.datetime64("1979-01-03T00:00:00.00"),
            time_step=np.timedelta64(1, "D"),
            verbosity=3,
            height_m=5,
        )

        # todo: this is testing instantiation, but not physics
        self.cnp = PRMSCanopy(prms_params, atm)
        self.cnp.advance(itime_step=0)
        self.cnp.calculate(time_length=1.0)

        return


class TestPRMSCanopyDomain():

    def test_init(self, domain):
        prms_params = PrmsParameters.load(domain["param_file"])

        control_file = domain["control_file"]
        control = ControlVariables.load(control_file)

        # this doesn't work.  cannot pass in cbh files
        #input_files_dict = domain["input_files_dict"]
        #atm = NHMBoundaryLayer(input_files_dict)
        start_time = control.control.start_time
        end_time = control.control.end_time
        initial_deltat = control.control.initial_deltat

        atm_init_test_dict = {
            "start_time": start_time,
            "end_time": end_time,
            "time_step": initial_deltat,
            "verbosity": 3,
            "height_m": 5,
        }
        print(domain["cbh_nc"])
        atm = NHMBoundaryLayer(domain["cbh_nc"], **atm_init_test_dict)
        atm.calculate_sw_rad_degree_day(prms_params)
        atm.calculate_potential_et_jh(prms_params)

        self.cnp = PRMSCanopy(prms_params, atm)
        self.cnp.advance(itime_step=0)

        while atm.current_time < atm.end_time:
            print(f"Running canopy for time: {atm.current_time}")
            self.cnp.advance(0)
            self.cnp.calculate(1.0)
            atm.advance()

        return
