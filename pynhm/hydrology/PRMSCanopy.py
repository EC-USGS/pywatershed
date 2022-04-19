from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.control import Control
from ..base.adapter import Adapter, adapter_factory

adaptable = Union[str, np.ndarray, Adapter]

NEARZERO = 1.0e-6
DNEARZERO = np.finfo(float).eps  # EPSILON(0.0D0)

RAIN = 0
SNOW = 1

BARESOIL = 0
GRASSES = 1

OFF = 0
ACTIVE = 1

LAND = 1
LAKE = 2


class PRMSCanopy(StorageUnit):
    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        pkwater_equiv: adaptable,
        transp_on: adaptable,
        hru_ppt: adaptable,
        hru_rain: adaptable,
        hru_snow: adaptable,
        potet: adaptable,
        verbose: bool = False,
    ):

        verbose = True
        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
        )

        self.name = "PRMSCanopy"
        # self._current_time = self.control.current_time

        # store dependencies
        self._input_variables_dict = {}
        self._input_variables_dict["pkwater_equiv"] = adapter_factory(
            pkwater_equiv, "pkwater_equiv"
        )
        self._input_variables_dict["transp_on"] = adapter_factory(
            transp_on, "transp_on"
        )
        self._input_variables_dict["hru_ppt"] = adapter_factory(
            hru_ppt, "hru_ppt"
        )
        self._input_variables_dict["hru_rain"] = adapter_factory(
            hru_rain, "hru_rain"
        )
        self._input_variables_dict["hru_snow"] = adapter_factory(
            hru_snow, "hru_snow"
        )
        self._input_variables_dict["potet"] = adapter_factory(potet, "potet")

        # Where does the initial storage come from? Document here.
        # apparently it's just zero?
        # self.inctp_stor = self.intcp_stor_init.copy()
        self.intcp_stor[:] = np.zeros([1])[0]
        self.intcp_stor_old = None

        return

    @staticmethod
    def get_parameters() -> tuple:
        """
        Return a list of the parameters required for this process

        """
        return (
            "nhru",
            "hru_area",
            "cov_type",
            "covden_sum",
            "covden_win",
            "srain_intcp",
            "wrain_intcp",
            "snow_intcp",
            "epan_coef",
            "potet_sublim",
            "snow_intcp",
        )

    @staticmethod
    def get_input_variables() -> tuple:
        """Get canopy reservoir input variables

        Returns:
            variables: input variables

        """
        return (
            "pkwater_equiv",
            "transp_on",
            "hru_ppt",
            "hru_rain",
            "hru_snow",
            "potet",
        )

    @staticmethod
    def get_variables() -> tuple:
        """Get canopy self variables

        Returns:
            variables: self variables

        """
        return (
            "intcp_stor",
            "net_rain",
            "net_snow",
            "intcp_evap",
            "rainfall_adj",
            "snowfall_adj",
            "potet",
            "interception_form",
            "intcp_transp_on",
        )

    def advance(self):

        # if self._current_time >= self.control.current_time:
        #    return None
        # self._current_time = self.control.current_time

        self.intcp_stor_old = self.intcp_stor
        self._itime_step += 1

        # set variables that depend on transpiration on/off setting
        # todo: this is currently hardwired to be on
        # if self.tranpiration_on:
        #    self.covden = self.covden_sum
        #    self.stor_max_rain = self.srain_intcp
        # else:
        #    self.covden = self.covden_win
        #    self.stor_max_rain = self.wrain_intcp

        for key, value in self._input_variables_dict.items():
            value.advance()
            v = getattr(self, key)
            v[:] = value.current

        self.interception_form[:] = RAIN
        idx = np.where(self.hru_snow > 0)
        self.interception_form[idx] = SNOW

        assert self.pkwater_equiv.shape == (self.nhru,)

        return

    def calculate(self, time_length, vectorized=False):
        if vectorized:
            self.calculate_vectorized(time_length)
        else:
            self.calculate_procedural(time_length)
        return

    def calculate_procedural(self, time_length):

        # todo: verify that hru_ppt is prcp
        hru_ppt = self.hru_ppt
        potet = self.potet
        hru_rain = self.hru_rain
        hru_snow = self.hru_snow

        net_rain = np.array(self.nhru * [0.0])
        net_snow = np.array(self.nhru * [0.0])
        net_ppt = np.array(self.nhru * [0.0])
        intcp_form = np.array(self.nhru * [RAIN])

        hru_type = np.array(self.nhru * [LAND])
        transp_on = self.transp_on
        intcp_evap = np.array(self.nhru * [0.0])
        hru_intcpstor = np.array(self.nhru * [0.0])

        for i in range(self.nhru):
            harea = self.hru_area[i]
            netrain = hru_rain[i]
            netsnow = hru_snow[i]

            if transp_on[i] == ACTIVE:
                cov = self.covden_sum[i]
                stor_max_rain = self.srain_intcp[i]
            else:
                cov = self.covden_win[i]
                stor_max_rain = self.wrain_intcp[i]

            intcp_form[i] = RAIN
            if hru_snow[i] > 0.0:
                intcp_form[i] = SNOW

            intcpstor = self.intcp_stor[i]
            intcpevap = 0.0
            changeover = 0.0
            extra_water = 0.0

            # lake or bare ground hrus
            if hru_type[i] == LAKE or self.cov_type[i] == BARESOIL:
                if self.cov_type[i] == BARESOIL and intcpstor > 0.0:
                    extra_water = self.intcp_stor[i]
                intcpstor = 0.0

            # ***** go from summer to winter cover density
            if transp_on[i] == OFF and self.intcp_transp_on[i] == ACTIVE:
                self.intcp_transp_on[i] = OFF
                if intcpstor > 0.0:
                    diff = self.covden_sum[i] - cov
                    changeover = intcpstor * diff
                    if cov > 0.0:
                        if changeover < 0.0:
                            intcpstor = intcpstor * self.covden_sum[i] / cov
                            changeover = 0.0
                    else:
                        intcpstor = 0.0

            # ****** go from winter to summer cover density, excess = throughfall
            elif transp_on[i] == ACTIVE and self.intcp_transp_on[i] == OFF:
                self.intcp_transp_on[i] = ACTIVE
                if intcpstor > 0.0:
                    diff = self.covden_win[i] - cov
                    changeover = intcpstor * diff
                    if cov > 0.0:
                        if changeover < 0.0:
                            intcpstor = intcpstor * self.covden_win[i] / cov
                            changeover = 0.0
                    else:
                        intcpstor = 0.0

            # *****Determine the amount of interception from rain
            # IF ( Hru_type(i)/=LAKE .AND. Cov_type(i)/=BARESOIL ) THEN
            if hru_type[i] != LAKE and self.cov_type[i] != BARESOIL:
                if hru_rain[i] > 0.0:
                    if cov > 0.0:
                        # IF ( Cov_type(i)>GRASSES ) THEN
                        if self.cov_type[i] > GRASSES:
                            # intercept(Hru_rain(i), stor_max_rain, cov, intcpstor, netrain)
                            intcpstor, netrain = self.intercept(
                                hru_rain[i],
                                stor_max_rain,
                                cov,
                                intcpstor,
                                netrain,
                            )
                        elif self.cov_type[i] == GRASSES:
                            # if there is no snowpack and no snowfall, then apparently, grasses
                            # can intercept rain.
                            # IF ( Pkwater_equiv(i)<DNEARZERO .AND. netsnow<NEARZERO ) THEN
                            if (
                                self.pkwater_equiv[i] < DNEARZERO
                                and netsnow < NEARZERO
                            ):
                                intcpstor, netrain = self.intercept(
                                    hru_rain[i],
                                    stor_max_rain,
                                    cov,
                                    intcpstor,
                                    netrain,
                                )

            # ******Determine amount of interception from snow
            if hru_snow[i] > 0.0:
                if cov > 0.0:
                    if self.cov_type[i] > GRASSES:
                        intcpstor, netsnow = self.intercept(
                            hru_snow[i],
                            self.snow_intcp[i],
                            cov,
                            intcpstor,
                            netsnow,
                        )
                        if netsnow < NEARZERO:
                            netrain = netrain + netsnow
                            netsnow = 0.0
                            # todo: deal with newsnow and pptmix?
                            # Newsnow(i) = OFF
                            # Pptmix(i) = OFF

            # todo: canopy application of irrigation water based on irr_type

            # ******compute evaporation or sublimation of interception

            # if precipitation assume no evaporation or sublimation
            if intcpstor > 0.0:
                if hru_ppt[i] < NEARZERO:
                    epan_coef = 1.0
                    evrn = potet[i] / epan_coef
                    evsn = potet[i] * self.potet_sublim[i]

                    # todo: pan et
                    # IF ( Use_pandata==ACTIVE ) THEN
                    # evrn = Pan_evap(Hru_pansta(i))
                    # IF ( evrn<0.0 ) evrn = 0.0
                    # ENDIF

                    if intcp_form[i] == SNOW:
                        z = intcpstor - evsn
                        if z > 0:
                            intcpstor = z
                            intcpevap = evsn
                        else:
                            intcpevap = intcpstor
                            intcpstor = 0.0
                    else:
                        d = intcpstor - evrn
                        if d > 0.0:
                            intcpstor = d
                            intcpevap = evrn
                        else:
                            intcpevap = intcpstor
                            intcpstor = 0.0

            if intcpevap * cov > potet[i]:
                last = intcpevap
                if cov > 0.0:
                    intcpevap = potet[i] / cov
                else:
                    intcpevap = 0.0
                intcpstor = intcpstor + last - intcpevap

            self.intcp_evap[i] = intcpevap
            self.intcp_stor[i] = intcpstor
            hru_intcpstor[i] = intcpstor * cov
            self.net_rain[i] = netrain
            self.net_snow[i] = netsnow
            net_ppt[i] = netrain + netsnow

        return

    def calculate_vectorized(self, time_length):

        # Retrieve atmospheric forcings - rename?
        rainfall_adj = rainfall
        snowfall_adj = snowfall
        potet = potet
        prcp = prcp

        # initialize calculation variables
        net_rain = rainfall_adj.copy()
        net_snow = snowfall_adj.copy()
        intcp_stor = self.intcp_stor_old.copy()
        intcp_evap = np.array(self.nhru * [0.0])
        net_ppt = np.array(self.nhru * [0.0])

        # todo: Lakes not handled; but lakes are not in NHM so probably okay

        # todo: Handle changeover water going from summer to winter

        # todo: Handle changeover water going from winter to summer

        # update interception and net_rain
        idx = np.where(
            (self.cov_type != BARESOIL)
            & (self.covden > 0)
            & (self.cov_type > GRASSES)
        )
        self.update_net_precip(
            rainfall_adj,
            self.stor_max_rain,
            self.covden,
            intcp_stor,
            net_rain,
            idx,
        )
        idx = np.where(
            (self.cov_type == GRASSES)
            & (self.covden > 0.0)
            & (self.pkwater_equiv < NEARZERO)
            & (snowfall_adj < NEARZERO)
        )
        self.update_net_precip(
            rainfall_adj,
            self.stor_max_rain,
            self.covden,
            intcp_stor,
            self.net_rain,
            idx,
        )

        # Update intcp_stor and net_snow with snowfall for anything greater than grass
        idx = np.where((self.cov_type > GRASSES) & (self.covden > 0.0))
        self.update_net_precip(
            snowfall_adj,
            self.stor_max_rain,
            self.covden,
            intcp_stor,
            self.net_snow,
            idx,
        )

        # todo: Handle irrigation water?  Depends on whether or not this is part of NHM

        # todo: epan_coef is supposed to be specified by month; here it is fixed to 1.0
        epan_coef = 1.0

        # # if there is precip, then shut off potet and sublimation
        # evrn = np.where(prcp < NEARZERO, potet / epan_coef, 0.)
        # evsn = np.where(prcp < NEARZERO, potet * self.potet_sublim, 0.)
        # intcp_stor_save = intcp_stor.copy()
        # depth = np.where(
        #     self.interception_form == SNOW,
        #     intcp_stor - evsn,
        #     intcp_stor - evrn,
        # )
        # intcp_stor = np.maximum(depth, 0.0)
        # intcp_evap = intcp_stor_save - intcp_stor

        for i in range(self.nhru):
            intcpstor = intcp_stor[i]
            intcpevap = 0.0
            if intcpstor > 0.0:
                if prcp[i] < NEARZERO:
                    evrn = potet[i] / epan_coef
                    evsn = potet[i] * self.potet_sublim[i]

                    if self.interception_form[i] == SNOW:
                        z = intcpstor - evsn
                        if z > 0:
                            intcpstor = z
                            intcpevap = evsn
                        else:
                            intcpevap = intcpstor
                            intcpstor = 0.0
                    else:
                        d = intcpstor - evrn
                        if d > 0.0:
                            intcpstor = d
                            intcpevap = evrn
                        else:
                            intcpevap = intcpstor
                            intcpstor = 0.0
            intcp_stor[i] = intcpstor
            intcp_evap[i] = intcpevap

        # todo: adjust intcp_evap for cover density
        # todo: but this doesn't seem to make any sense
        # IF ( intcpevap*cov>Potet(i) ) THEN
        #  last = intcpevap
        #  IF ( cov>0.0 ) THEN
        #    intcpevap = Potet(i)/cov
        #  ELSE
        #    intcpevap = 0.0
        #  ENDIF
        #  intcpstor = intcpstor + last - intcpevap
        # ENDIF

        # accumulate into net_ppt
        net_ppt[:] = self.net_rain + self.net_snow

        self.intcp_stor[:] = intcp_stor[:]
        hru_intcpstor = intcp_stor * self.covden

        return

    @staticmethod
    def intercept(precip, stor_max, cov, intcp_stor, net_precip):
        net_precip = precip * (1.0 - cov)
        intcp_stor = intcp_stor + precip
        if intcp_stor > stor_max:
            net_precip = net_precip + (intcp_stor - stor_max) * cov
            intcp_stor = stor_max
        return intcp_stor, net_precip

    @staticmethod
    def update_net_precip(
        precip, stor_max, covden, intcp_stor, net_precip, idx
    ):
        net_precip[idx] = precip[idx] * (1.0 - covden[idx])
        intcp_stor[idx] += precip[idx]
        for i in idx[0]:
            if intcp_stor[i] > stor_max[i]:
                net_precip[i] += (intcp_stor[i] - stor_max[i]) * covden[i]
                intcp_stor[i] = stor_max[i]
        return


def canopy_vectorized(
    intcpstor,
    hru_rain,
    hru_snow,
    transp_on,
    covden_sum,
    covden_win,
    srain_intcp,
    wrain_intcp,
):
    netrain = hru_rain
    netsnow = hru_snow
    if transp_on == ACTIVE:
        cov = covden_sum
        stor_max_rain = srain_intcp
    else:
        cov = covden_win
        stor_max_rain = wrain_intcp
    intcpform = RAIN
    if hru_snow > 0.0:
        intcpform = SNOW
    intcpevap = 0.0
    changeover = 0.0
    extra_water = 0.0

    return intcpstor, netrain, netsnow
