import numpy as np

from ..atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from ..base.storageUnit import StorageUnit
from ..utils.parameters import PrmsParameters

RAIN = 0
SNOW = 1


class PRMSCanopy(StorageUnit):
    def __init__(
        self,
        params: PrmsParameters,
        atm: NHMBoundaryLayer,
    ):

        verbose = True
        super().__init__("cnp", id, params, atm, verbose)

        # Do we somehow name variables so that we can create
        # budget table?
        # inflows, outflows, dS, residual

        # define self variables that will be used for the calculation
        # todo: may need way to initialize interception storage to something non-zero
        self.intcp_stor_old = np.array(self.nhru * [0.0])
        self.intcp_stor = np.array(self.nhru * [0.0])
        self.tranpiration_on = True  # prms variable Transp_on
        self.covden = None  # will be set to covden_sum or covden_win
        self.interception_form = np.array(self.nhru * [RAIN], dtype=int)

        # define information on the output data that will be created
        self.output_column_names = ["date"]
        self.output_data = []
        if "nhm_id" in params.parameters:
            self.output_column_names += [
                f"nhru_hru_intcpstor_{nhmid}"
                for nhmid in params.parameters["nhm_id"]
            ]
        else:
            self.output_column_names += [
                f"intcpstor_{id}" for id in range(self.nhru)
            ]
        return

    @staticmethod
    def get_required_parameters() -> list:
        """
        Returt a list of the paramaters required for this process

        """
        return [
            "nhru",
            "hru_area",
            "covden_sum",
            "covden_win",
            "srain_intcp",
            "wrain_intcp",
            "snow_intcp",
            "epan_coef",
            "potet_sublim",
        ]

    def advance(self, itime_step):
        self.intcp_stor_old = self.intcp_stor

        # set variables that depend on transpiration on/off setting
        # todo: this is currently hardwired to be on
        if self.tranpiration_on:
            self.covden = self.covden_sum
            self.stor_max_rain = self.srain_intcp
        else:
            self.covden = self.covden_win
            self.stor_max_rain = self.wrain_intcp

        self.interception_form[:] = RAIN
        snowfall = self.atm.get_current_state("snowfall")
        idx = np.where(snowfall > 0)
        self.interception_form[idx] = SNOW

        return

    def calculate(self, time_length):

        # Retrieve atmospheric forcings
        hru_rain = self.atm.get_current_state("rainfall")
        hru_snow = self.atm.get_current_state("snowfall")
        potet = self.atm.get_current_state("potet")
        prcp = self.atm.get_current_state("prcp")

        # initialize calculation variables
        intcp_stor = self.intcp_stor_old.copy()
        net_rain = np.array(self.nhru * [0.0])
        net_snow = np.array(self.nhru * [0.0])
        net_ppt = np.array(self.nhru * [0.0])

        # todo: add these as needed
        # intcp_evap = np.array(self.nhru * [0.0])
        # change_over = np.array(self.nhru * [0.0])
        # extra_water = np.array(self.nhru * [0.0])

        # todo: Lakes not handled; but lakes are not in NHM so probably okay

        # todo: Handle changeover water going from summer to winter

        # todo: Handle changeover water going from winter to summer

        # Update intcp_stor and net_rain with rainfall
        self.update_net_precip(
            hru_rain, self.stor_max_rain, self.covden, intcp_stor, net_rain
        )

        # Update intcp_stor and net_snow with snowfall
        self.update_net_precip(
            hru_snow, self.snow_intcp, self.covden, intcp_stor, net_snow
        )

        # todo: Handle irrigation water?  Depends on whether or not this is part of NHM

        # todo: Handle evaporation and sublimination
        # todo: epan_coef is supposed to be specified by month; here it is fixed to 1.0
        epan_coef = 1.0

        # if there precip, then shut off potet and sublimation
        evrn = np.where(prcp > 0.0, 0.0, potet / epan_coef)
        evsn = np.where(prcp > 0.0, 0.0, potet * self.potet_sublim)
        intcp_stor_save = intcp_stor.copy()
        depth = np.where(
            self.interception_form == SNOW,
            intcp_stor - evsn,
            intcp_stor - evrn,
        )
        intcp_stor = np.maximum(depth, 0.0)
        intcp_evap = intcp_stor_save - intcp_stor

        # accumulate into net_ppt
        net_ppt[:] = net_rain + net_snow

        self.intcp_stor[:] = intcp_stor[:]
        hru_intcpstor = intcp_stor * self.covden
        self.output_data.append([self.atm.current_time] + list(hru_intcpstor))

        return

    @staticmethod
    def update_net_precip(precip, stor_max, covden, intcp_stor, net_precip):
        net_precip[:] = precip * (1.0 - covden)
        intcp_stor[:] += precip[:]
        idx = np.where(intcp_stor > stor_max)
        net_precip[idx] += (intcp_stor[idx] - stor_max[idx]) * covden[idx]
        intcp_stor[idx] = stor_max[idx]
        return
