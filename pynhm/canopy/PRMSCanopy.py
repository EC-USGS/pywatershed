import numpy as np

from ..atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from ..base.storageUnit import StorageUnit
from ..utils.parameters import PrmsParameters

RAIN = 0
SNOW = 1


class PRMSCanopy(StorageUnit):
    @staticmethod
    def get_required_parameters():
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

    def __init__(
        self,
        params: PrmsParameters,
        atm: NHMBoundaryLayer,
    ):

        verbose = True
        super().__init__("cnp", id, params, atm, verbose)

        # potential variables are managed by a parent StateAccess class.  This cannot be
        # set until after super().__init__() called to set it up
        self._potential_variables = ["intcp_stor"]

        # Do we somehow name variables so that we can create
        # budget table?
        # inflows, outflows, dS, residual

        # volume_start = 0.  # do we need to set this from input?   intcp_stor_start * covden
        # self.intcp_stor_max = intcp_stor_max
        # self.intcp_stor_new = intcp_stor_start

        # todo: may need way to interception storage to something non-zero
        self.intcp_stor_old = np.array(self.nhru * [0.0])
        self["intcp_stor"] = np.array(self.nhru * [0.0])
        self.tranpiration_on = True  # prms variable Transp_on
        self.covden = (
            None  # prms variable will be set to covden_sum or covden_win
        )
        self.interception_form = np.array(self.nhru * [RAIN], dtype=int)

        self.output_column_names = [
            "date",
            "precip",
            "aet",
            "intcp_ds",
            "net_precip",
            "residual",
            "intcp_stor_new",
            "intcp_stor_old",
        ]
        return

    def advance(self, itime_step):
        super().advance(itime_step)
        self.intcp_stor_old = self["intcp_stor"]

        # set variables that depend on transpiration on/off setting
        # todo: this is hardwired to be on
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

        # initialize calculation variables
        intcp_stor = self["intcp_stor"]
        net_rain = np.array(self.nhru * [0.0])
        net_snow = np.array(self.nhru * [0.0])

        # todo: add these as needed
        # intcp_evap = np.array(self.nhru * [0.0])
        # change_over = np.array(self.nhru * [0.0])
        # extra_water = np.array(self.nhru * [0.0])

        # todo: Lakes not handled; but not in NHM so probably okay

        # todo: Handle changeover water going from summer to winter

        # todo: Handle changeover water going from winter to summer

        # Interception from rain
        self.update_net_precip(
            hru_rain, self.stor_max_rain, self.covden, intcp_stor, net_rain
        )

        # todo: Interception from snow
        self.update_net_precip(
            hru_snow, self.snow_intcp, self.covden, intcp_stor, net_snow
        )

        # todo: Handle irrigation water?  Depends on whether or not this is part of NHM

        # todo: Handle evaporation and sublimination

        return

    @staticmethod
    def update_net_precip(precip, stor_max, covden, intcp_stor, net_precip):
        net_precip[:] = precip * (1.0 - covden)
        intcp_stor[:] += precip[:]
        idx = np.where(intcp_stor > stor_max)
        net_precip[idx] += (intcp_stor[idx] - stor_max[idx]) * covden[idx]
        return

    def calculate_old(self, time_length):

        # Retrieve forcings
        # precip = self.forcing.precip_current
        # pot_et = self.forcing.pot_et_current
        precip = self.atm.get_current_state("prcp")
        potet = self.atm.get_current_state("potet")

        # Calculate canopy ET as the minimum of the canopy storage volume
        # and the potential ET.  Then let the forcing module know that
        intcp_stor = self.intcp_stor_old
        requested_et = min(intcp_stor, pot_et)
        canopy_et = self.forcing.consume_pot_et(requested_et)
        intcp_stor -= canopy_et

        # calculate change in storage in canopy and precipitation throughfall
        available_storage = self.intcp_stor_max - intcp_stor
        if precip > available_storage:
            precipitation_throughfall = precip - available_storage
            intcp_stor = self.intcp_stor_max
        else:
            precipitation_throughfall = 0.0
            intcp_stor += precip
        self.intcp_stor_new = intcp_stor
        used_storage = self.intcp_stor_new - self.intcp_stor_old
        used_storage *= self.area * self.covden

        # net precipitation
        net_precipitation = (
            precip * (1 - self.covden)
            + precipitation_throughfall * self.covden
        )
        net_precipitation *= self.area

        # add and remove volumetric flow rates
        self.add_water("precipitation", precip * self.area, time_length)
        self.remove_water(
            "aet", -canopy_et * self.area * self.covden, time_length
        )
        self.remove_water("intcp_ds", -used_storage, time_length)
        self.remove_water("net_precipitation", -net_precipitation, time_length)

        # Send net_precipitation to any registered recipients
        for recipient_info in self.recipients:
            recipient, process_name = recipient_info
            if process_name == "net_precipitation":
                recipient.add_water(
                    "net_precipitation", net_precipitation, time_length
                )

        # Calculate residual and save results
        super().calculate(time_length)
        output = [self.forcing.current_date]
        output.append(precip * self.area)
        output.append(-canopy_et * self.area * self.covden)
        output.append(-used_storage)
        output.append(-net_precipitation)
        output.append(self.residual_new)
        output.append(self.intcp_stor_new)
        output.append(self.intcp_stor_old)
        self.output_data.append(output)
        return
