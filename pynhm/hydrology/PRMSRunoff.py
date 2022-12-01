import numba as nb
import numpy as np

from pynhm.base.storageUnit import StorageUnit

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import HruType, nan, zero

NEARZERO = 1.0e-6
DNEARZERO = np.finfo(float).eps  # EPSILON(0.0D0)

RAIN = 0
SNOW = 1

BARESOIL = 0
GRASSES = 1

OFF = 0
ACTIVE = 1

LAND = HruType.LAND.value
LAKE = HruType.LAKE.value


class PRMSRunoff(StorageUnit):
    """PRMS surface runoff."""

    def __init__(
        self,
        control: Control,
        soil_moist_prev: adaptable,
        net_rain: adaptable,
        net_ppt: adaptable,
        net_snow: adaptable,
        potet: adaptable,
        snowmelt: adaptable,
        snow_evap: adaptable,
        pkwater_equiv: adaptable,
        pptmix_nopack: adaptable,
        snowcov_area: adaptable,
        through_rain: adaptable,
        hru_intcpevap: adaptable,
        intcp_changeover: adaptable,
        budget_type: str = None,
        calc_method: str = None,
        verbose: bool = False,
        load_n_time_batches: int = 1,
    ) -> "PRMSRunoff":

        super().__init__(
            control=control,
            verbose=verbose,
            load_n_time_batches=load_n_time_batches,
        )
        self.name = "PRMSRunoff"

        self._calc_method = str(calc_method)

        self._set_inputs(locals())

        self._set_budget(budget_type)

        return

    def _set_initial_conditions(self):
        # Where does the initial storage come from? Document here.
        # apparently it's just zero?
        # self.var1_stor[:] = np.zeros([1])[0]
        # self.var1_stor_old = None

        # cdl -- todo:
        # this variable is calculated and stored by PRMS but does not seem
        # to be used widely
        self.dprst_in = np.zeros(self.nhru, dtype=float)
        self.dprst_vol_open_max = np.zeros(self.nhru, dtype=float)
        self.dprst_vol_clos_max = np.zeros(self.nhru, dtype=float)
        self.dprst_frac_clos = np.zeros(self.nhru, dtype=float)
        self.dprst_vol_thres_open = np.zeros(self.nhru, dtype=float)

        # call the basin_init hack to calculate basin
        # variables
        self.basin_init()

        # call the depression storage init
        self.dprst_init()

        return

    @staticmethod
    def get_parameters() -> tuple:
        """
        Return a list of the parameters required for this process

        """
        return (
            "nhru",
            "hru_type",
            "hru_area",
            "hru_percent_imperv",
            "imperv_stor_max",
            "dprst_frac",
            "carea_max",
            "smidx_coef",
            "smidx_exp",
            "soil_moist_max",
            "snowinfil_max",
            "dprst_depth_avg",
            "dprst_et_coef",
            "dprst_flow_coef",
            "dprst_frac",
            "dprst_frac_init",
            "dprst_frac_open",
            "dprst_seep_rate_clos",
            "dprst_seep_rate_open",
            "sro_to_dprst_imperv",
            "sro_to_dprst_perv",
            "va_open_exp",
            "va_clos_exp",
            "op_flow_thres",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get input variables

        Returns:
            variables: input variables

        """
        return (
            "soil_moist_prev",
            "net_rain",
            "net_ppt",
            "net_snow",
            "potet",
            "snowmelt",
            "snow_evap",
            "pkwater_equiv",
            "pptmix_nopack",
            "snowcov_area",
            "through_rain",
            "hru_intcpevap",
            "intcp_changeover",
        )

    @staticmethod
    def get_init_values() -> dict:
        """Get initial values

        Returns:
            dict: initial values for named variables
        """
        return {
            "contrib_fraction": zero,
            "infil": zero,
            "infil_hru": zero,
            "sroff": zero,  # todo: privatize and only make vol public
            "sroff_vol": zero,
            "hru_sroffp": zero,
            "hru_sroffi": zero,
            "imperv_stor": zero,
            "imperv_evap": zero,
            "hru_impervevap": zero,
            "hru_impervstor": zero,
            "hru_impervstor_old": zero,
            "hru_impervstor_change": zero,
            "dprst_vol_frac": zero,
            "dprst_vol_clos": zero,
            "dprst_vol_open": zero,
            "dprst_vol_clos_frac": zero,
            "dprst_vol_open_frac": zero,
            "dprst_area_clos": zero,
            "dprst_area_open": zero,
            "dprst_area_clos_max": zero,
            "dprst_area_open_max": zero,
            "dprst_sroff_hru": zero,
            "dprst_seep_hru": nan,
            "dprst_evap_hru": zero,
            "dprst_insroff_hru": zero,
            "dprst_stor_hru": zero,
            "dprst_stor_hru_old": zero,
            "dprst_stor_hru_change": zero,
        }

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                # "net_rain",
                "through_rain",
                "snowmelt",
                "intcp_changeover",
            ],
            "outputs": [
                # sroff = hru_sroffi + hru_sroffp + dprst_sroff_hru
                "hru_sroffi",
                "hru_sroffp",
                "dprst_sroff_hru",
                "infil_hru",
                "hru_impervevap",
                "dprst_seep_hru",
                "dprst_evap_hru",
            ],
            "storage_changes": [
                "hru_impervstor_change",
                "dprst_stor_hru_change",
            ],
        }

    def basin_init(self):
        """
        This is trying to replicate the prms basin_init function that
        calculates some of the variables needed here by runoff.  This should
        probably go somewhere else at some point as I suspect other components
        may need similar information.
        """
        dprst_flag = ACTIVE
        self.hru_perv = np.zeros(self.nhru, float)
        self.hru_frac_perv = np.zeros(self.nhru, float)
        self.hru_imperv = np.zeros(self.nhru, float)
        self.dprst_area_max = np.zeros(self.nhru, float)
        for k in range(self.nhru):
            i = k
            harea = self.hru_area[k]
            perv_area = harea
            if self.hru_percent_imperv[k] > 0.0:
                self.hru_imperv[i] = self.hru_percent_imperv[i] * harea
                perv_area = perv_area - self.hru_imperv[i]

            if dprst_flag == ACTIVE:
                self.dprst_area_max[i] = self.dprst_frac[i] * harea
                if self.dprst_area_max[i] > 0.0:
                    self.dprst_area_open_max[i] = (
                        self.dprst_area_max[i] * self.dprst_frac_open[i]
                    )
                    self.dprst_frac_clos[i] = 1.0 - self.dprst_frac_open[i]
                    self.dprst_area_clos_max[i] = (
                        self.dprst_area_max[i] - self.dprst_area_open_max[i]
                    )
                    if self.dprst_area_clos_max[i] > 0.0:
                        dprst_clos_flag = ACTIVE
                        # cdl -- todo: above variable should be stored to self?
                    if self.dprst_area_open_max[i] > 0.0:
                        dprst_open_flag = ACTIVE
                        # cdl -- todo: above variable should be stored to self?
                    perv_area = perv_area - self.dprst_area_max[i]

            self.hru_perv[i] = perv_area
            self.hru_frac_perv[i] = perv_area / harea
        return

    def dprst_init(self):
        for j in range(self.nhru):
            i = j
            if self.dprst_frac[i] > 0.0:

                if self.dprst_depth_avg[i] == 0.0:
                    raise Exception(
                        f"dprst_fac > and dprst_depth_avg == 0 for HRU {i}"
                    )

                # calculate open and closed volumes (acre-inches) of depression
                # storage by HRU
                # Dprst_area_open_max is the maximum open depression area
                # (acres) that can generate surface runoff:
                dprst_clos_flag = ACTIVE
                if dprst_clos_flag == ACTIVE:
                    self.dprst_vol_clos_max[i] = (
                        self.dprst_area_clos_max[i] * self.dprst_depth_avg[i]
                    )
                dprst_open_flag = ACTIVE
                if dprst_open_flag == ACTIVE:
                    self.dprst_vol_open_max[i] = (
                        self.dprst_area_open_max[i] * self.dprst_depth_avg[i]
                    )

                # calculate the initial open and closed depression storage
                # volume:
                dprst_open_flag = ACTIVE
                if dprst_open_flag == ACTIVE:
                    self.dprst_vol_open[i] = (
                        self.dprst_frac_init[i] * self.dprst_vol_open_max[i]
                    )
                dprst_clos_flag = ACTIVE
                if dprst_clos_flag == ACTIVE:
                    self.dprst_vol_clos[i] = (
                        self.dprst_frac_init[i] * self.dprst_vol_clos_max[i]
                    )

                # threshold volume is calculated as the % of maximum open
                # depression storage above which flow occurs *  total open
                # depression storage volume
                self.dprst_vol_thres_open[i] = (
                    self.op_flow_thres[i] * self.dprst_vol_open_max[i]
                )
                if self.dprst_vol_open[i] > 0.0:
                    open_vol_r = (
                        self.dprst_vol_open[i] / self.dprst_vol_open_max[i]
                    )
                    if open_vol_r < NEARZERO:
                        frac_op_ar = 0.0
                    elif open_vol_r > 1.0:
                        frac_op_ar = 1.0
                    else:
                        frac_op_ar = np.exp(
                            self.va_open_exp[i] * np.log(open_vol_r)
                        )
                    self.dprst_area_open[i] = (
                        self.dprst_area_open_max[i] * frac_op_ar
                    )
                    if self.dprst_area_open[i] > self.dprst_area_open_max[i]:
                        self.dprst_area_open[i] = self.dprst_area_open_max[i]

                # Closed depression surface area for each HRU:
                if self.dprst_vol_clos[i] > 0.0:
                    clos_vol_r = (
                        self.dprst_vol_clos[i] / self.dprst_vol_clos_max[i]
                    )
                    if clos_vol_r < NEARZERO:
                        frac_cl_ar = 0.0
                    elif clos_vol_r > 1.0:
                        frac_cl_ar = 1.0
                    else:
                        frac_cl_ar = np.exp(
                            self.va_clos_exp[i] * np.log(clos_vol_r)
                        )
                    self.dprst_area_clos[i] = (
                        self.dprst_area_clos_max[i] * frac_cl_ar
                    )
                    if self.dprst_area_clos[i] > self.dprst_area_clos_max[i]:
                        self.dprst_area_clos[i] = self.dprst_area_clos_max[i]

                #
                # calculate basin open and closed depression storage volumes
                self.dprst_stor_hru[i] = (
                    self.dprst_vol_open[i] + self.dprst_vol_clos[i]
                ) / self.hru_area[i]
                self.dprst_stor_hru_old[i] = self.dprst_stor_hru[i]

                if (
                    self.dprst_vol_open_max[i] + self.dprst_vol_clos_max[i]
                    > 0.0
                ):
                    self.dprst_vol_frac[i] = (
                        self.dprst_vol_open[i] + self.dprst_vol_clos[i]
                    ) / (
                        self.dprst_vol_open_max[i] + self.dprst_vol_clos_max[i]
                    )
        return

    def _advance_variables(self) -> None:
        """Advance the variables
        Returns:
            None
        """
        self.hru_impervstor_old[:] = self.hru_impervstor
        self.dprst_stor_hru_old[:] = self.dprst_stor_hru
        return None

    def _calculate(self, time_length, vectorized=False):
        """Calculate terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """
        if self._calc_method.lower() == "numba":
            import numba as nb

            # if not hasattr(self, "_calculate_numba"):

        elif self._calc_method.lower() in ["none", "numpy"]:

            self.calculate_prms_style()

        # <
        self.infil_hru[:] = self.infil * self.hru_frac_perv

        self.hru_impervstor_change[:] = (
            self.hru_impervstor - self.hru_impervstor_old
        )
        self.dprst_stor_hru_change[:] = (
            self.dprst_stor_hru - self.dprst_stor_hru_old
        )

        self.sroff_vol = self.sroff * self.control.params.hru_in_to_cf

        return

    def calculate_prms_style(self):

        # move towards replacing net_rain and other variables.
        # self.net_rain[:] = self.through_rain

        dprst_chk = 0
        self.infil[:] = 0.0
        for k in range(self.nhru):

            # TODO: remove duplicated vars
            # TODO: move setting constants outside the loop.

            # cdl i = Hru_route_order(k)
            i = k
            ihru = i

            runoff = 0.0

            hruarea = self.hru_area[i]
            perv_area = self.hru_perv[i]
            perv_frac = self.hru_frac_perv[i]
            srp = 0.0
            sri = 0.0
            self.hru_sroffp[i] = 0.0
            self.contrib_fraction[i] = 0.0
            hruarea_imperv = self.hru_imperv[i]
            imperv_frac = 0.0
            if hruarea_imperv > 0.0:
                imperv_frac = self.hru_percent_imperv[i]
                self.hru_sroffi[i] = 0.0
                self.imperv_evap[i] = 0.0
                self.hru_impervevap[i] = 0.0

            avail_et = (
                self.potet[i] - self.snow_evap[i] - self.hru_intcpevap[i]
            )
            availh2o = self.intcp_changeover[i] + self.net_rain[i]

            (
                sri,
                srp,
                self.imperv_stor[i],
                self.infil[i],
                self.contrib_fraction[i],
            ) = self.compute_infil(
                self.contrib_fraction[i],
                self.soil_moist_prev[i],
                self.soil_moist_max[i],
                self.carea_max[i],
                self.smidx_coef[i],
                self.smidx_exp[i],
                self.pptmix_nopack[i],
                self.net_rain[i],
                self.net_ppt[i],
                self.imperv_stor[i],
                self.imperv_stor_max[i],
                self.snowmelt[i],
                self.snowinfil_max[i],
                self.net_snow[i],
                self.pkwater_equiv[i],
                self.infil[i],
                self.hru_type[i],
                self.intcp_changeover[i],
                hruarea_imperv,
                sri,
                srp,
                self.check_capacity,
                self.perv_comp,
            )

            dprst_flag = ACTIVE  # cdl todo: hardwired
            frzen = OFF  # cdl todo: hardwired

            if dprst_flag == ACTIVE:
                self.dprst_in[i] = 0.0
                self.dprst_seep_hru[i] = zero
                dprst_chk = OFF
                # JLM: can this logic be moved inside dprst_comp?
                if self.dprst_area_max[i] > 0.0:
                    dprst_chk = ACTIVE
                    if frzen == OFF:
                        (
                            self.dprst_in[i],
                            self.dprst_vol_open[i],
                            avail_et,
                            self.dprst_vol_clos[i],
                            self.dprst_sroff_hru[i],
                            srp,
                            sri,
                            self.dprst_evap_hru[i],
                            self.dprst_seep_hru[i],
                        ) = self.dprst_comp(
                            self.dprst_vol_clos[i],
                            self.dprst_area_clos_max[i],
                            self.dprst_area_clos[i],
                            self.dprst_vol_open_max[i],
                            self.dprst_vol_open[i],
                            self.dprst_area_open_max[i],
                            self.dprst_area_open[i],
                            self.dprst_sroff_hru[i],
                            self.dprst_seep_hru[i],
                            self.sro_to_dprst_perv[i],
                            self.sro_to_dprst_imperv[i],
                            self.dprst_evap_hru[i],
                            avail_et,
                            availh2o,
                            self.dprst_in[i],
                            ihru,
                            srp,
                            sri,
                            imperv_frac,
                            perv_frac,
                        )
                        runoff = runoff + self.dprst_sroff_hru[i] * hruarea

            # cdl -- the upper part of this block needs to be done to calculate
            #        runoff and srunoff
            # Compute runoff for pervious and impervious area, and depression
            # storage area
            srunoff = 0.0
            if self.hru_type[i] == LAND:
                runoff = runoff + srp * perv_area + sri * hruarea_imperv
                srunoff = runoff / hruarea
                self.hru_sroffp[i] = srp * perv_frac

            # <
            # cdl -- the guts of this was implemented in the python code below
            # Compute evaporation from impervious area
            if hruarea_imperv > 0.0:
                if self.imperv_stor[i] > 0.0:
                    self.imperv_stor[i], self.imperv_evap[i] = self.imperv_et(
                        self.imperv_stor[i],
                        self.potet[i],
                        self.imperv_evap[i],
                        self.snowcov_area[i],
                        avail_et,
                        imperv_frac,
                    )
                    self.hru_impervevap[i] = self.imperv_evap[i] * imperv_frac
                    avail_et = avail_et - self.hru_impervevap[i]
                    if avail_et < 0.0:
                        self.hru_impervevap[i] = (
                            self.hru_impervevap[i] + avail_et
                        )
                        if self.hru_impervevap[i] < 0.0:
                            self.hru_impervevap[i] = 0.0
                        self.imperv_evap[i] = self.imperv_evap[i] / imperv_frac
                        self.imperv_stor[i] = (
                            self.imperv_stor[i] - avail_et / imperv_frac
                        )
                        avail_et = 0.0
                    self.hru_impervstor[i] = self.imperv_stor[i] * imperv_frac
                self.hru_sroffi[i] = sri * imperv_frac

            # <
            # cdl -- saving sroff here
            if dprst_chk == ACTIVE:
                self.dprst_stor_hru[i] = (
                    self.dprst_vol_open[i] + self.dprst_vol_clos[i]
                ) / hruarea
            self.sroff[i] = srunoff

        # <
        return

    @staticmethod
    def compute_infil(
        contrib_fraction,
        soil_moist_prev,
        soil_moist_max,
        carea_max,
        smidx_coef,
        smidx_exp,
        pptmix_nopack,
        net_rain,
        net_ppt,
        imperv_stor,
        imperv_stor_max,
        snowmelt,
        snowinfil_max,
        net_snow,
        pkwater_equiv,
        infil,
        hru_type,
        intcp_changeover,
        hruarea_imperv,
        sri,
        srp,
        check_capacity,
        perv_comp,
    ):

        isglacier = False  # todo -- hardwired
        cascade_active = False
        hru_flag = 0
        if hru_type == LAND or isglacier:
            hru_flag = 1
        if cascade_active:
            raise Exception("bad bad bad")
        else:
            avail_water = 0.0

        # compute runoff from canopy changeover water
        if intcp_changeover > 0.0:
            avail_water = avail_water + intcp_changeover
            infil = infil + intcp_changeover
            if hru_flag == 1:
                infil, srp, contrib_fraction = perv_comp(
                    soil_moist_prev,
                    carea_max,
                    smidx_coef,
                    smidx_exp,
                    intcp_changeover,
                    intcp_changeover,
                    infil,
                    srp,
                )

        # if rain/snow event with no antecedent snowpack,
        # compute the runoff from the rain first and then proceed with the
        # snowmelt computations
        if pptmix_nopack == ACTIVE:
            avail_water = avail_water + net_rain
            infil = infil + net_rain
            if hru_flag == 1:
                infil, srp, contrib_fraction = perv_comp(
                    soil_moist_prev,
                    carea_max,
                    smidx_coef,
                    smidx_exp,
                    net_rain,
                    net_rain,
                    infil,
                    srp,
                )

        # If precipitation on snowpack, all water available to the surface is
        # considered to be snowmelt, and the snowmelt infiltration
        # procedure is used.  If there is no snowpack and no precip,
        # then check for melt from last of snowpack.  If rain/snow mix
        # with no antecedent snowpack, compute snowmelt portion of runoff.
        if snowmelt > 0.0:
            # There was no snowmelt but a snowpack may exist.  If there is
            # no snowpack then check for rain on a snowfree HRU.
            avail_water = avail_water + snowmelt
            infil = infil + snowmelt
            if hru_flag == 1:
                if pkwater_equiv > 0.0 or net_ppt - net_snow < NEARZERO:
                    # Pervious area computations
                    infil, srp = check_capacity(
                        soil_moist_prev,
                        soil_moist_max,
                        snowinfil_max,
                        infil,
                        srp,
                    )
                else:
                    # Snowmelt occurred and depleted the snowpack
                    infil, srp, contrib_fraction = perv_comp(
                        soil_moist_prev,
                        carea_max,
                        smidx_coef,
                        smidx_exp,
                        snowmelt,
                        net_ppt,
                        infil,
                        srp,
                    )

        elif pkwater_equiv < DNEARZERO:
            # If no snowmelt and no snowpack but there was net snow then
            # snowpack was small and was lost to sublimation.
            if net_snow < NEARZERO and net_rain > 0.0:
                avail_water = avail_water + net_rain
                infil = infil + net_rain
                if hru_flag == 1:
                    infil, srp, contrib_fraction = perv_comp(
                        soil_moist_prev,
                        carea_max,
                        smidx_coef,
                        smidx_exp,
                        net_rain,
                        net_rain,
                        infil,
                        srp,
                    )

        # Snowpack exists, check to see if infil exceeds maximum daily
        # snowmelt infiltration rate. Infil results from rain snow mix
        # on a snowfree surface.
        elif infil > 0.0:
            if hru_flag == 1:
                infil, srp = check_capacity(
                    soil_moist_prev,
                    soil_moist_max,
                    snowinfil_max,
                    infil,
                    srp,
                )

        if hruarea_imperv > 0.0:
            imperv_stor = imperv_stor + avail_water
            if hru_flag == 1:
                if imperv_stor > imperv_stor_max:
                    sri = imperv_stor - imperv_stor_max
                    imperv_stor = imperv_stor_max

        return sri, srp, imperv_stor, infil, contrib_fraction

    def dprst_comp(
        self,
        dprst_vol_clos,
        dprst_area_clos_max,
        dprst_area_clos,
        dprst_vol_open_max,
        dprst_vol_open,
        dprst_area_open_max,
        dprst_area_open,
        dprst_sroff_hru,
        dprst_seep_hru,
        sro_to_dprst_perv,
        sro_to_dprst_imperv,
        dprst_evap_hru,
        avail_et,
        net_rain,
        dprst_in,
        ihru,
        srp,
        sri,
        imperv_frac,
        perv_frac,
    ):

        cascade_flag = OFF  # cdl -- todo: hardwired
        if cascade_flag > OFF:
            raise Exception("i am brokin")
        else:
            inflow = 0.0

        if self.pptmix_nopack[ihru]:
            inflow = inflow + net_rain

        # I f precipitation on snowpack all water available to the surface is
        # considered to be snowmelt
        # If there is no snowpack and no precip,then check for melt from last
        # of snowpack.
        # If rain/snow mix with no antecedent snowpack, compute snowmelt
        # portion of runoff.
        if self.snowmelt[ihru]:
            inflow = inflow + self.snowmelt[ihru]

        # !******There was no snowmelt but a snowpack may exist.  If there is
        # !******no snowpack then check for rain on a snowfree HRU.
        elif self.pkwater_equiv[ihru] < DNEARZERO:
            # !      If no snowmelt and no snowpack but there was net snow then
            # !      snowpack was small and was lost to sublimation.
            if self.net_snow[ihru] < NEARZERO and net_rain > 0.0:
                inflow = inflow + net_rain

        dprst_add_water_use = OFF
        if dprst_add_water_use == ACTIVE:
            raise Exception("not supported")

        dprst_in = 0.0
        if dprst_area_open_max > 0.0:
            dprst_in = inflow * dprst_area_open_max
            dprst_vol_open = dprst_vol_open + dprst_in

        if dprst_area_clos_max > 0.0:
            tmp1 = inflow * dprst_area_clos_max
            dprst_vol_clos = dprst_vol_clos + tmp1
            dprst_in = dprst_in + tmp1
        dprst_in = dprst_in / self.hru_area[ihru]

        # add any pervious surface runoff fraction to depressions
        dprst_srp = 0.0
        dprst_sri = 0.0
        if srp > 0.0:
            tmp = srp * perv_frac * sro_to_dprst_perv * self.hru_area[ihru]
            if dprst_area_open_max > 0.0:
                dprst_srp_open = tmp * self.dprst_frac_open[ihru]
                dprst_srp = dprst_srp_open / self.hru_area[ihru]
                dprst_vol_open = dprst_vol_open + dprst_srp_open
            if dprst_area_clos_max > 0.0:
                dprst_srp_clos = tmp * self.dprst_frac_clos[ihru]
                dprst_srp = dprst_srp + dprst_srp_clos / self.hru_area[ihru]
                dprst_vol_clos = dprst_vol_clos + dprst_srp_clos
            srp = srp - dprst_srp / perv_frac
            if srp < 0.0:
                if srp < -NEARZERO:
                    srp = 0.0

        if sri > 0.0:
            tmp = sri * imperv_frac * sro_to_dprst_imperv * self.hru_area[ihru]
            if dprst_area_open_max > 0.0:
                dprst_sri_open = tmp * self.dprst_frac_open[ihru]
                dprst_sri = dprst_sri_open / self.hru_area[ihru]
                dprst_vol_open = dprst_vol_open + dprst_sri_open
            if dprst_area_clos_max > 0.0:
                dprst_sri_clos = tmp * self.dprst_frac_clos[ihru]
                dprst_sri = dprst_sri + dprst_sri_clos / self.hru_area[ihru]
                dprst_vol_clos = dprst_vol_clos + dprst_sri_clos
            sri = sri - dprst_sri / imperv_frac
            if sri < 0.0:
                if sri < -NEARZERO:
                    sri = 0.0
            self.dprst_insroff_hru[ihru] = dprst_srp + dprst_sri

        dprst_area_open = 0.0
        if dprst_vol_open > 0.0:
            open_vol_r = dprst_vol_open / dprst_vol_open_max
            if open_vol_r < NEARZERO:
                frac_op_ar = 0.0
            elif open_vol_r > 1.0:
                frac_op_ar = 1.0
            else:
                frac_op_ar = np.exp(
                    self.va_open_exp[ihru] * np.log(open_vol_r)
                )
            dprst_area_open = dprst_area_open_max * frac_op_ar
            if dprst_area_open > dprst_area_open_max:
                dprst_area_open = dprst_area_open_max

        if dprst_area_clos_max > 0.0:
            dprst_area_clos = 0.0
            if dprst_vol_clos > 0.0:
                clos_vol_r = dprst_vol_clos / self.dprst_vol_clos_max[ihru]
                if clos_vol_r < NEARZERO:
                    frac_cl_ar = 0.0
                elif clos_vol_r > 1.0:
                    frac_cl_ar = 1.0
                else:
                    frac_cl_ar = np.exp(
                        self.va_clos_exp[ihru] * np.log(clos_vol_r)
                    )
                dprst_area_clos = dprst_area_clos_max * frac_cl_ar
                if dprst_area_clos > dprst_area_clos_max:
                    dprst_area_clos = dprst_area_clos_max
                if dprst_area_clos < NEARZERO:
                    dprst_area_clos = 0.0

        #
        # evaporate water from depressions based on snowcov_area
        # dprst_evap_open & dprst_evap_clos = inches-acres on the HRU
        unsatisfied_et = avail_et
        dprst_avail_et = (
            self.potet[ihru]
            * (1.0 - self.snowcov_area[ihru])
            * self.dprst_et_coef[ihru]
        )
        dprst_evap_hru = 0.0
        if dprst_avail_et > 0.0:
            dprst_evap_open = 0.0
            dprst_evap_clos = 0.0
            if dprst_area_open > 0.0:
                dprst_evap_open = min(
                    dprst_area_open * dprst_avail_et, dprst_vol_open
                )
                if dprst_evap_open / self.hru_area[ihru] > unsatisfied_et:
                    dprst_evap_open = unsatisfied_et * self.hru_area[ihru]
                if dprst_evap_open > dprst_vol_open:
                    dprst_evap_open = dprst_vol_open
                unsatisfied_et = (
                    unsatisfied_et - dprst_evap_open / self.hru_area[ihru]
                )
                dprst_vol_open = dprst_vol_open - dprst_evap_open

            if dprst_area_clos > 0.0:
                dprst_evap_clos = min(
                    dprst_area_clos * dprst_avail_et, dprst_vol_clos
                )
                if dprst_evap_clos / self.hru_area[ihru] > unsatisfied_et:
                    dprst_evap_clos = unsatisfied_et * self.hru_area[ihru]
                if dprst_evap_clos > dprst_vol_clos:
                    dprst_evap_clos = dprst_vol_clos
                dprst_vol_clos = dprst_vol_clos - dprst_evap_clos

            dprst_evap_hru = (
                dprst_evap_open + dprst_evap_clos
            ) / self.hru_area[ihru]

        # compute seepage
        dprst_seep_hru = 0.0
        if dprst_vol_open > 0.0:
            seep_open = dprst_vol_open * self.dprst_seep_rate_open[ihru]
            dprst_vol_open = dprst_vol_open - seep_open
            if dprst_vol_open < 0.0:
                seep_open = seep_open + dprst_vol_open
                dprst_vol_open = 0.0
            dprst_seep_hru = seep_open / self.hru_area[ihru]

        # compute open surface runoff
        dprst_sroff_hru = 0.0
        if dprst_vol_open > 0.0:
            dprst_sroff_hru = max(0.0, dprst_vol_open - dprst_vol_open_max)
            dprst_sroff_hru = dprst_sroff_hru + max(
                0.0,
                (
                    dprst_vol_open
                    - dprst_sroff_hru
                    - self.dprst_vol_thres_open[ihru]
                )
                * self.dprst_flow_coef[ihru],
            )
            dprst_vol_open = dprst_vol_open - dprst_sroff_hru
            dprst_sroff_hru = dprst_sroff_hru / self.hru_area[ihru]
            if dprst_vol_open < 0.0:
                dprst_vol_open = 0.0

        if dprst_area_clos_max > 0.0:
            if dprst_area_clos > NEARZERO:
                seep_clos = dprst_vol_clos * self.dprst_seep_rate_clos[ihru]
                dprst_vol_clos = dprst_vol_clos - seep_clos
                if dprst_vol_clos < 0.0:
                    seep_clos = seep_clos + dprst_vol_clos
                    dprst_vol_clos = 0.0
                dprst_seep_hru = (
                    dprst_seep_hru + seep_clos / self.hru_area[ihru]
                )

        avail_et = avail_et - dprst_evap_hru
        if dprst_vol_open_max > 0.0:
            self.dprst_vol_open_frac[ihru] = (
                dprst_vol_open / dprst_vol_open_max
            )
        if self.dprst_vol_clos_max[ihru] > 0.0:
            self.dprst_vol_clos_frac[ihru] = (
                dprst_vol_clos / self.dprst_vol_clos_max[ihru]
            )
        if dprst_vol_open_max + self.dprst_vol_clos_max[ihru] > 0.0:
            self.dprst_vol_frac[ihru] = (dprst_vol_open + dprst_vol_clos) / (
                dprst_vol_open_max + self.dprst_vol_clos_max[ihru]
            )
        self.dprst_stor_hru[ihru] = (
            dprst_vol_open + dprst_vol_clos
        ) / self.hru_area[ihru]

        return (
            dprst_in,
            dprst_vol_open,
            avail_et,
            dprst_vol_clos,
            dprst_sroff_hru,
            srp,
            sri,
            dprst_evap_hru,
            dprst_seep_hru,
        )

    @staticmethod
    def perv_comp(
        soil_moist_prev,
        carea_max,
        smidx_coef,
        smidx_exp,
        pptp,
        ptc,
        infil,
        srp,
    ):
        """Pervious area computations."""
        smidx_module = True
        if smidx_module:
            smidx = soil_moist_prev + 0.5 * ptc
            if smidx > 25.0:
                ca_fraction = carea_max
            else:
                ca_fraction = smidx_coef * 10.0 ** (smidx_exp * smidx)
        else:
            raise Exception("you did a bad thing...")

        if ca_fraction > carea_max:
            ca_fraction = carea_max
        srpp = ca_fraction * pptp
        infil = infil - srpp
        srp = srp + srpp

        return infil, srp, ca_fraction

    @staticmethod
    def check_capacity(
        soil_moist_prev, soil_moist_max, snowinfil_max, infil, srp
    ):
        """
        Fill soil to soil_moist_max, if more than capacity restrict
        infiltration by snowinfil_max, with excess added to runoff
        """
        capacity = soil_moist_max - soil_moist_prev
        excess = infil - capacity
        if excess > snowinfil_max:
            srp = srp + excess - snowinfil_max
            infil = snowinfil_max + capacity
        return infil, srp

    @staticmethod
    def imperv_et(imperv_stor, potet, imperv_evap, sca, avail_et, imperv_frac):
        if sca < 1.0:
            if potet < imperv_stor:
                imperv_evap = potet * (1.0 - sca)
            else:
                imperv_evap = imperv_stor * (1.0 - sca)
            if imperv_evap * imperv_frac > avail_et:
                imperv_evap = avail_et / imperv_frac
            imperv_stor = imperv_stor - imperv_evap
        return imperv_stor, imperv_evap
