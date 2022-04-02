from datetime import datetime, timedelta

import numpy as np
import pandas as pd

import pynhm.preprocess
from pynhm.atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from pynhm.canopy.PRMSCanopy import PRMSCanopy
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


class TestPRMSCanopySimple:
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
            "epan_coef": np.array(nhru * [1.0]),
            "potet_sublim": np.array(nhru * [1.0]),
        }
        prms_params = PrmsParameters(prms_params)
        atm = NHMBoundaryLayer(
            forcings_dict,
            prms_params,
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
        print(domain["cbh_nc"])
        atm = NHMBoundaryLayer.load_netcdf(
            domain["cbh_nc"], prms_params, **atm_information_dict
        )
        atm.calculate_sw_rad_degree_day()
        atm.calculate_potential_et_jh()

        self.cnp = PRMSCanopy(prms_params, atm)
        self.cnp.advance(itime_step=0)

        for istep in range(atm.n_time):
            if istep > 0:
                atm.advance()
            # print(f"Running canopy for step {istep} and day: {atm.current_time}")
            self.cnp.advance(0)
            self.cnp.calculate(1.0)

        # create data frame of prms interception storage (from nhru_hru_intcpstor.csv)
        output_files = domain["prms_outputs"]
        fname = output_files["intcpstor"]
        print(f"loading {fname}")
        intcpstor_prms = pynhm.preprocess.CsvFile(fname)
        intcpstor_prms = intcpstor_prms.to_dataframe()

        # create data frame of pynhm interception storage
        intcpstor_pynhm = pd.DataFrame(
            self.cnp.output_data, columns=self.cnp.output_column_names
        )
        intcpstor_pynhm.set_index("date", inplace=True)

        print("intcp_stor  min  max")
        a1 = intcpstor_prms.to_numpy()
        a2 = intcpstor_pynhm.to_numpy()
        print(f"prms  {a1.min()}    {a1.max()}")
        print(f"pynhm  {a2.min()}    {a2.max()}")

        makeplot = False
        if makeplot:
            import matplotlib.pyplot as plt

            f = plt.figure()
            ax = plt.subplot(1, 1, 1, aspect="equal")
            xmin = 1e30
            xmax = -1e30
            for colname in intcpstor_prms.columns:
                x1 = intcpstor_prms.loc[:, colname].values
                x2 = intcpstor_pynhm.loc[:, colname].values
                ax.scatter(
                    x1, x2, facecolors="none", edgecolors="k", linewidth=0.25
                )
                xmin = min(xmin, x1.min(), x2.min())
                xmax = max(xmax, x1.max(), x2.max())
            ax.set_title("Interception Storage")
            ax.set_xlabel("PRMS Interception Storage, in inches")
            ax.set_ylabel("PYNHM Interception Storage, in inches")
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(xmin, xmax)
            pname = domain["domain_name"] + "_1-1.png"
            print(f"Creating plot {pname}")
            plt.savefig(pname, dpi=300)

            f = plt.figure()
            ax = plt.subplot(1, 1, 1)
            pc = plt.imshow(a1 - a2, cmap="jet")
            plt.colorbar(pc, shrink=0.5)
            ax.set_title("Difference Between PRMS and PYNHM")
            ax.set_xlabel("HRU number")
            ax.set_ylabel("Time step")
            pname = domain["domain_name"] + "_2d.png"
            print(f"Creating plot {pname}")
            plt.savefig(pname, dpi=300)

        return
