from typing import Literal
from warnings import warn

import numpy as np
from numba import prange

from ..base.adapter import adaptable
from ..base.conservative_process import ConservativeProcess
from ..base.control import Control
from ..constants import (
    CovType,
    HruType,
    dnearzero,
    nan,
    nearzero,
    numba_num_threads,
    zero,
)
from ..parameters import Parameters

try:
    from ..prms_canopy_f import canopy

    _calculate_fortran = canopy.calc_canopy
    has_prmscanopy_f = True
except ImportError:
    has_prmscanopy_f = False

# set constants (may need .value for enum to be used in > comparisons)
BARESOIL = CovType.BARESOIL.value
GRASSES = CovType.GRASSES.value
LAND = HruType.LAND.value
LAKE = HruType.LAKE.value
RAIN = 0
SNOW = 1
OFF = 0
ACTIVE = 1


class PRMSCanopy(ConservativeProcess):
    """PRMS canopy class.

    A canopy or vegetation representation from PRMS.

    Implementation based on PRMS 5.2.1 with theoretical documentation given in
    the PRMS-IV documentation:

    `Markstrom, S. L., Regan, R. S., Hay, L. E., Viger, R. J., Webb, R. M.,
    Payn, R. A., & LaFontaine, J. H. (2015). PRMS-IV, the
    precipitation-runoff modeling system, version 4. US Geological Survey
    Techniques and Methods, 6, B7.
    <https://pubs.usgs.gov/tm/6b7/pdf/tm6-b7.pdf>`__


    Args:
        control: a Control object
        discretization: a discretization of class Parameters
        parameters: a parameter object of class Parameters
        pk_ice_prev: Previous snowpack ice on each HRU
        freeh2o_prev: Previous snowpack free water on each HRU
        transp_on: Flag indicating whether transpiration is occurring
            (0=no;1=yes)
        hru_ppt: Precipitation on each HRU
        hru_rain: Rain on each HRU
        hru_snow: Snow on each HRU
        budget_type: one of [None, "warn", "error"]
        calc_method: one of ["fortran", "numba", "numpy"]. None defaults to
            "numba".
        verbose: Print extra information or not?
        load_n_time_batches: not-implemented
    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        pk_ice_prev: adaptable,
        freeh2o_prev: adaptable,
        transp_on: adaptable,
        hru_ppt: adaptable,
        hru_rain: adaptable,
        hru_snow: adaptable,
        potet: adaptable,
        budget_type: Literal[None, "warn", "error"] = None,
        calc_method: Literal["fortran", "numba", "numpy"] = None,
        verbose: bool = None,
    ):
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "PRMSCanopy"

        # set hrutype to LAND as this is only type supported in NHM
        self._hru_type = np.array(self.nhru * [LAND])

        self._set_inputs(locals())
        self._set_options(locals())

        self._set_budget()
        self._init_calc_method()

        return

    @staticmethod
    def get_dimensions() -> tuple:
        return ("nhru",)

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "cov_type",
            "covden_sum",
            "covden_win",
            "srain_intcp",
            "wrain_intcp",
            "snow_intcp",
            "potet_sublim",
        )

    @staticmethod
    def get_inputs() -> tuple:
        return (
            "pk_ice_prev",
            "freeh2o_prev",
            "transp_on",
            "hru_ppt",
            "hru_rain",
            "hru_snow",
            "potet",
        )

    @staticmethod
    def get_init_values() -> dict:
        return {
            "net_ppt": zero,
            "net_rain": zero,
            "net_snow": zero,
            "intcp_changeover": zero,
            "intcp_evap": zero,
            "intcp_form": nan,  # = RAIN bool would have to match RAIN/SNOW
            "intcp_stor": zero,
            "intcp_transp_on": 0,  # could make boolean
            "hru_intcpevap": zero,
            "hru_intcpstor": zero,
            "hru_intcpstor_change": zero,
            "hru_intcpstor_old": zero,
        }

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": ["hru_rain", "hru_snow"],
            "outputs": [
                "net_rain",
                "net_snow",
                "hru_intcpevap",
                "intcp_changeover",  # ?? always zero
            ],
            "storage_changes": ["hru_intcpstor_change"],
        }

    def _set_initial_conditions(self):
        return

    def _init_calc_method(self):
        if self._calc_method is None:
            self._calc_method = "numba"

        avail_methods = ["numpy", "numba", "fortran"]
        fortran_msg = ""
        if self._calc_method == "fortran" and not has_prmscanopy_f:
            _ = avail_methods.remove("fortran")
            fortran_msg = "\n(Fortran not available as installed)\n"

        if self._calc_method.lower() not in avail_methods:
            msg = (
                f"Invalid calc_method={self._calc_method} for {self.name}. "
                f"{fortran_msg}"
                f"Setting calc_method to 'numba' for {self.name}"
            )
            warn(msg)
            self._calc_method = "numba"

        if self._calc_method.lower() in ["numba"]:
            import numba as nb

            numba_msg = f"{self.name} jit compiling with numba "
            nb_parallel = (numba_num_threads is not None) and (
                numba_num_threads > 1
            )
            if nb_parallel:
                numba_msg += f"and using {numba_num_threads} threads"
            print(numba_msg, flush=True)

            # JLM: note. I gave up on specifying signatures because it
            #      appears impossible/undocumented how to specify the type
            #      of a passed function. The work is not in vain as then
            #      types will be required for f90. The signatures remain
            #      here commented for that reason and in case we can get it
            #      to work in the future.

            # self._intercept_numba = nb.njit(
            #     nb.types.Tuple(
            #         (
            #             nb.float64[:],  # intcp_stor
            #             nb.float64[:],  # net_precip
            #         )
            #     )(
            #         nb.float64[:],  # precip
            #         nb.float64[:],  # stor_max
            #         nb.float64[:],  # cov
            #         nb.float64[:],  # intcp_stor
            #         nb.float64[:],  # net_precip
            #     ),
            #     fastmath=True,
            # )(self._intercept)
            self._intercept = nb.njit(fastmath=True)(self._intercept)

            # self._calculate_numba = nb.njit(
            #     nb.types.Tuple(
            #         (
            #             nb.float64[:],  # intcp_evap
            #             nb.float64[:],  # intcp_stor
            #             nb.float64[:],  # net_rain
            #             nb.float64[:],  # net_snow
            #             nb.float64[:],  # net_ppt
            #             nb.float64[:],  # hru_intcpstor
            #             nb.float64[:],  # hru_intcpevap
            #             nb.float64[:],  # intcp_changeover
            #             nb.int32[:],  # intcp_transp_on
            #         )
            #     )(
            #         nb.int32,  # np.int32(self.nhru)
            #         nb.int64[:],  # cov_type
            #         nb.float64[:],  # covden_sum
            #         nb.float64[:],  # covden_win
            #         nb.float64[:],  # hru_intcpstor
            #         nb.float64[:],  # hru_intcpevap
            #         nb.float64[:],  # hru_ppt
            #         nb.float64[:],  # hru_rain
            #         nb.float64[:],  # hru_snow
            #         nb.float64[:],  # intcp_changeover
            #         nb.float64[:],  # intcp_evap
            #         nb.int32[:],  # intcp_form
            #         nb.float64[:],  # intcp_stor
            #         nb.int32[:],  # intcp_transp_on
            #         nb.float64[:],  # net_ppt
            #         nb.float64[:],  # net_rain
            #         nb.float64[:],  # net_snow
            #         nb.float64[:],  # self.pk_ice_prev
            #         nb.float64[:],  # self.freeh2o_prev
            #         nb.float64[:],  # potet
            #         nb.float64[:],  # potet_sublim
            #         nb.float64[:],  # snow_intcp
            #         nb.float64[:],  # srain_intcp
            #         nb.float64[:],  # transp_on
            #         nb.float64[:],  # wrain_intcp
            #         nb.float64,  # np.float64(time_length),
            #         nb.int64[:],  # self._hru_type
            #         nb.float64,  # nearzero
            #         nb.float64,  # dnearzero,
            #         nb.int32,  # BARESOIL,
            #         nb.int32,  # GRASSES,
            #         nb.int32,  # LAND,
            #         nb.int32,  # LAKE,
            #         nb.int32,  # RAIN,
            #         nb.int32,  # SNOW,
            #         nb.int32,  # OFF,
            #         nb.int32,  # ACTIVE,
            #         nb.typeof(self._intercept_numba),  # function.
            #     ),
            #     fastmath=True,
            # )(self._calculate_procedural)
            self._calculate_canopy = nb.njit(
                fastmath=True, parallel=nb_parallel
            )(self._calculate_numpy)

        elif self._calc_method.lower() == "fortran":
            pass
            # fortran has a different call signature in the last agument
            # because the intercept function is not passed.
            # so it is handled with an if statement at call time.
            # self._calculate_gw = _calculate_fortran

        else:
            self._calculate_canopy = self._calculate_numpy

        return

    def _advance_variables(self):
        """Advance canopy
        Returns:
            None

        """
        self.hru_intcpstor_old[:] = self.hru_intcpstor
        return

    def _calculate(self, time_length):
        """Calculate canopy terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """
        if self._calc_method.lower() != "fortran":
            (
                self.intcp_evap[:],
                self.intcp_form[:],
                self.intcp_stor[:],
                self.net_rain[:],
                self.net_snow[:],
                self.net_ppt[:],
                self.hru_intcpstor[:],
                self.hru_intcpevap[:],
                self.intcp_changeover[:],
                self.intcp_transp_on[:],
            ) = self._calculate_canopy(
                nhru=np.int32(self.nhru),
                cov_type=self.cov_type,
                covden_sum=self.covden_sum,
                covden_win=self.covden_win,
                hru_intcpstor=self.hru_intcpstor,
                hru_intcpevap=self.hru_intcpevap,
                hru_ppt=self.hru_ppt,
                hru_rain=self.hru_rain,
                hru_snow=self.hru_snow,
                intcp_changeover=self.intcp_changeover,
                intcp_evap=self.intcp_evap,
                intcp_stor=self.intcp_stor,
                intcp_transp_on=self.intcp_transp_on,
                net_ppt=self.net_ppt,
                net_rain=self.net_rain,
                net_snow=self.net_snow,
                pk_ice_prev=self.pk_ice_prev,
                freeh2o_prev=self.freeh2o_prev,
                potet=self.potet,
                potet_sublim=self.potet_sublim,
                snow_intcp=self.snow_intcp,
                srain_intcp=self.srain_intcp,
                transp_on=self.transp_on,
                wrain_intcp=self.wrain_intcp,
                time_length=time_length,
                hru_type=self._hru_type,
                nearzero=nearzero,
                dnearzero=dnearzero,
                BARESOIL=np.int32(BARESOIL),
                GRASSES=np.int32(GRASSES),
                LAND=np.int32(LAND),
                LAKE=np.int32(LAKE),
                RAIN=np.int32(RAIN),
                SNOW=np.int32(SNOW),
                OFF=np.int32(OFF),
                ACTIVE=np.int32(ACTIVE),
                intercept=self._intercept,
            )

        else:
            (
                self.intcp_evap[:],
                self.intcp_form[:],
                self.intcp_stor[:],
                self.net_rain[:],
                self.net_snow[:],
                self.net_ppt[:],
                self.hru_intcpstor[:],
                self.hru_intcpevap[:],
                self.intcp_changeover[:],
                self.intcp_transp_on[:],
            ) = _calculate_fortran(
                self.cov_type,
                self.covden_sum,
                self.covden_win,
                self.hru_intcpstor,
                self.hru_intcpevap,
                self.hru_ppt,
                self.hru_rain,
                self.hru_snow,
                self.intcp_changeover,
                self.intcp_evap,
                self.intcp_stor,
                self.intcp_transp_on,
                self.net_ppt,
                self.net_rain,
                self.net_snow,
                self.pk_ice_prev,
                self.freeh2o_prev,
                self.potet,
                self.potet_sublim,
                self.snow_intcp,
                self.srain_intcp,
                self.transp_on,
                self.wrain_intcp,
                time_length,
                self._hru_type,
                nearzero,
                dnearzero,
                np.int32(BARESOIL),
                np.int32(GRASSES),
                np.int32(LAND),
                np.int32(LAKE),
                np.int32(RAIN),
                np.int32(SNOW),
                np.int32(OFF),
                np.int32(ACTIVE),
            )

        # <
        self.hru_intcpstor_change[:] = (
            self.hru_intcpstor - self.hru_intcpstor_old
        )

        return

    @staticmethod
    def _calculate_numpy(
        nhru,
        cov_type,
        covden_sum,
        covden_win,
        hru_intcpstor,
        hru_intcpevap,
        hru_ppt,
        hru_rain,
        hru_snow,
        intcp_changeover,
        intcp_evap,
        intcp_stor,
        intcp_transp_on,
        net_ppt,
        net_rain,
        net_snow,
        pk_ice_prev,
        freeh2o_prev,
        potet,
        potet_sublim,
        snow_intcp,
        srain_intcp,
        transp_on,
        wrain_intcp,
        time_length,
        hru_type,
        nearzero,
        dnearzero,
        BARESOIL,
        GRASSES,
        LAND,
        LAKE,
        RAIN,
        SNOW,
        OFF,
        ACTIVE,
        intercept,
    ):
        # TODO: would be nice to alphabetize the arguments
        #       probably while keeping constants at the end.
        #       intcp_form and intcp_transp_on also do not appear to be
        #       actual inputs
        #       Keep the f90 call signature consistent with the args in
        #       python/numba.

        intcp_form = np.full_like(hru_rain, np.nan, dtype="int32")
        for i in prange(nhru):
            netrain = hru_rain[i]
            netsnow = hru_snow[i]

            if transp_on[i] == ACTIVE:
                cov = covden_sum[i]
                stor_max_rain = srain_intcp[i]
            else:
                cov = covden_win[i]
                stor_max_rain = wrain_intcp[i]

            intcp_form[i] = RAIN
            if hru_snow[i] > 0.0:
                intcp_form[i] = SNOW

            intcpstor = intcp_stor[i]
            intcpevap = 0.0
            changeover = 0.0
            extra_water = 0.0

            # lake or bare ground hrus
            if hru_type[i] == LAKE or cov_type[i] == BARESOIL:
                if cov_type[i] == BARESOIL and intcpstor > 0.0:
                    extra_water = intcp_stor[i]
                intcpstor = 0.0

            # ***** go from summer to winter cover density
            if transp_on[i] == OFF and intcp_transp_on[i] == ACTIVE:
                intcp_transp_on[i] = OFF
                if intcpstor > 0.0:
                    diff = covden_sum[i] - cov
                    changeover = intcpstor * diff
                    if cov > 0.0:
                        if changeover < 0.0:
                            intcpstor = intcpstor * covden_sum[i] / cov
                            changeover = 0.0
                    else:
                        intcpstor = 0.0

            # **** go from winter to summer cover density, excess = throughfall
            elif transp_on[i] == ACTIVE and intcp_transp_on[i] == OFF:
                intcp_transp_on[i] = ACTIVE
                if intcpstor > 0.0:
                    diff = covden_win[i] - cov
                    changeover = intcpstor * diff
                    if cov > 0.0:
                        if changeover < 0.0:
                            intcpstor = intcpstor * covden_win[i] / cov
                            changeover = 0.0
                    else:
                        intcpstor = 0.0

            # *****Determine the amount of interception from rain
            # IF ( Hru_type(i)/=LAKE .AND. Cov_type(i)/=BARESOIL ) THEN
            if hru_type[i] != LAKE and cov_type[i] != BARESOIL:
                if hru_rain[i] > 0.0:
                    if cov > 0.0:
                        # IF ( Cov_type(i)>GRASSES ) THEN
                        if cov_type[i] > GRASSES:
                            # intercept(
                            #     Hru_rain(i), stor_max_rain, cov, intcpstor,
                            #     netrain)
                            intcpstor, netrain = intercept(
                                hru_rain[i],
                                stor_max_rain,
                                cov,
                                intcpstor,
                                netrain,
                            )
                        elif cov_type[i] == GRASSES:
                            # if there is no snowpack and no snowfall, then
                            # apparently, grasses can intercept rain.
                            # IF (
                            #     pkwater_ante(i)<dnearzero
                            #     .AND. netsnow<nearzero ) THEN
                            if (
                                pk_ice_prev[i] + freeh2o_prev[i]
                            ) < dnearzero and netsnow < nearzero:
                                intcpstor, netrain = intercept(
                                    hru_rain[i],
                                    stor_max_rain,
                                    cov,
                                    intcpstor,
                                    netrain,
                                )

            # ******Determine amount of interception from snow
            if hru_snow[i] > 0.0:
                if cov > 0.0:
                    if cov_type[i] > GRASSES:
                        intcpstor, netsnow = intercept(
                            hru_snow[i],
                            snow_intcp[i],
                            cov,
                            intcpstor,
                            netsnow,
                        )
                        if netsnow < nearzero:
                            netrain = netrain + netsnow
                            netsnow = 0.0
                            # todo: deal with newsnow and pptmix?
                            # Newsnow(i) = OFF
                            # Pptmix(i) = OFF

            # todo: canopy application of irrigation water based on irr_type

            # ******compute evaporation or sublimation of interception

            # if precipitation assume no evaporation or sublimation
            if intcpstor > 0.0:
                if hru_ppt[i] < nearzero:
                    epan_coef = 1.0
                    evrn = potet[i] / epan_coef
                    evsn = potet[i] * potet_sublim[i]

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
            intcp_evap[i] = intcpevap
            intcp_stor[i] = intcpstor
            net_rain[i] = netrain
            net_snow[i] = netsnow
            net_ppt[i] = netrain + netsnow
            hru_intcpstor[i] = intcpstor * cov
            hru_intcpevap[i] = intcpevap * cov

            intcp_changeover[i] = changeover + extra_water

        return (
            intcp_evap,
            intcp_form,
            intcp_stor,
            net_rain,
            net_snow,
            net_ppt,
            hru_intcpstor,
            hru_intcpevap,
            intcp_changeover,
            intcp_transp_on,
        )

    @staticmethod
    def _intercept(precip, stor_max, cov, intcp_stor, net_precip):
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
