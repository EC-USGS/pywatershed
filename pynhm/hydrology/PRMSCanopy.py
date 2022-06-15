from typing import Union

import numpy as np

from pynhm.base.budget import Budget
from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control
from ..constants import CovType, HruType, zero

adaptable = Union[str, np.ndarray, Adapter]

# set constants (may need .value for enum to be used in > comparisons)
NEARZERO = 1.0e-6
DNEARZERO = np.finfo(float).eps  # EPSILON(0.0D0)
BARESOIL = CovType.BARESOIL.value
GRASSES = CovType.GRASSES.value
LAND = HruType.LAND
LAKE = HruType.LAKE
RAIN = 0
SNOW = 1
OFF = 0
ACTIVE = 1


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
        budget_type: str = None,
        verbose: bool = False,
    ):

        self.name = "PRMSCanopy"
        verbose = True
        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
            subclass_name=self.name,
        )

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

        if budget_type is None:
            self.budget = None
        else:
            budget_terms = {
                "inputs": {
                    "hru_rain": self.hru_rain,
                    "hru_snow": self.hru_snow,
                },
                "outputs": {
                    "net_rain": self.net_rain,
                    "net_snow": self.net_snow,
                    "hru_intcpevap": self.hru_intcpevap,
                },
                "storage_changes": {
                    "hru_intcpstor_change": self.hru_intcpstor_change,
                },
            }
            self.budget = Budget(
                self.control,
                **budget_terms,
                time_unit="D",
                description=self.name,
                imbalance_fatal=(budget_type == "strict"),
            )

        return

    def set_initial_conditions(self):
        # self.inctp_stor = self.intcp_stor_init.copy()
        # self.intcp_stor_old = None
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
    def get_inputs() -> tuple:
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
        """Get canopy output variables

        Returns:
            variables: output variables

        """
        return (
            "net_ppt",
            "net_rain",
            "net_snow",
            "intcp_stor",
            "intcp_evap",
            "hru_intcpstor",
            "hru_intcpstor_old",
            "hru_intcpstor_change",
            "hru_intcpevap",
            "intcp_form",
            "intcp_transp_on",  # this is private in prms6 and is not in the metadata
            # i defined metadata for it in a very adhoc way
        )

    @staticmethod
    def get_init_values() -> dict:
        """Get canopy initial values

        Returns:
            dict: initial values for named variables
        """

        return {
            "net_rain": zero,
            "net_snow": zero,
            "net_ppt": zero,
            "intcp_stor": zero,  # JLM: this is the only one that was explicitly defined
            "intcp_evap": zero,
            "hru_intcpstor": zero,
            "hru_intcpstor_old": zero,
            "hru_intcpstor_change": zero,
            "hru_intcpevap": zero,
            "intcp_form": 0,  # could make boolean but would have to make the RAIN/SNOW match
            "intcp_transp_on": 0,  # could make boolean
        }

    def _advance_variables(self):
        """Advance canopy
        Returns:
            None

        """
        # self.intcp_stor_old = self.intcp_stor
        self.hru_intcpstor_old[:] = self.hru_intcpstor
        return

    def calculate(self, time_length, vectorized=False):
        """Calculate canopy terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """
        if vectorized:
            # todo: not implemented yet
            self.calculate_vectorized(time_length)
        else:
            self.calculate_procedural(time_length)

        self.hru_intcpstor_change[:] = (
            self.hru_intcpstor - self.hru_intcpstor_old
        )

        if self.budget is not None:
            self.budget.advance()
            self.budget.calculate()

        return

    def calculate_procedural(self, time_length):

        # Assign short variables for input variables
        hru_ppt = self.hru_ppt
        potet = self.potet
        hru_rain = self.hru_rain
        hru_snow = self.hru_snow
        intcp_form = self.intcp_form
        transp_on = self.transp_on

        # set hrutype to LAND as this is only type supported in NHM
        hru_type = np.array(self.nhru * [LAND])

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

            # Store calculated values in output variables
            self.intcp_evap[i] = intcpevap
            self.intcp_stor[i] = intcpstor
            self.net_rain[i] = netrain
            self.net_snow[i] = netsnow
            self.net_ppt[i] = netrain + netsnow
            self.hru_intcpstor[i] = intcpstor * cov
            self.hru_intcpevap[i] = intcpevap * cov

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
