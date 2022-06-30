import numpy as np

from pynhm.base.storageUnit import StorageUnit

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import HruType, zero

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
        hru_intcpevap: adaptable,
        intcp_changeover: adaptable,
        budget_type: str = None,
        verbose: bool = False,
    ):
        self.name = "PRMSRunoff"
        super().__init__(
            control=control,
            verbose=verbose,
            subclass_name=self.name,
        )

        self.set_inputs(locals())

        self.set_budget(budget_type)
        # if budget_type is None:
        #     self.budget = None
        # else:
        #     budget_terms = {
        #         "inputs": {
        #             "net_rain": self.net_rain,
        #             "net_snow": self.net_snow,
        #             "snowmelt": self.snowmelt,
        #         },
        #         "outputs": {
        #             "hru_sroffi": self.hru_sroffi,
        #             "hru_sroffp": self.hru_sroffp,
        #             "infil": self.infil,
        #             "hru_impervevap": self.hru_impervevap,
        #             "dprst_sroff_hru": self.dprst_sroff_hru,
        #             "dprst_seep_hru": self.dprst_seep_hru,
        #             "dprst_evap_hru": self.dprst_evap_hru,
        #             # "dprst_insroff_hru": self.dprst_insroff_hru,
        #         },
        #         "storage_changes": {
        #             "hru_impervstor_change": self.hru_impervstor_change,
        #             "dprst_stor_hru_change": self.dprst_stor_hru_change,
        #         },
        #     }

        #     self.budget = Budget(
        #         self.control,
        #         **budget_terms,
        #         time_unit="D",
        #         description=self.name,
        #         imbalance_fatal=(budget_type == "strict"),
        #     )

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

    def set_initial_conditions(self):
        # Where does the initial storage come from? Document here.
        # apparently it's just zero?
        # self.var1_stor[:] = np.zeros([1])[0]
        # self.var1_stor_old = None
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
            "hru_intcpevap",
            "intcp_changeover",
        )

    @staticmethod
    def get_variables() -> tuple:
        """Get output variables

        Returns:
            variables: output variables

        """
        return (
            "contrib_fraction",
            "infil",
            "sroff",
            "hru_sroffp",
            "hru_sroffi",
            "imperv_stor",
            "imperv_evap",
            "hru_impervevap",
            "hru_impervstor",
            "hru_impervstor_old",
            "hru_impervstor_change",
            "dprst_vol_frac",
            "dprst_vol_clos",
            "dprst_vol_open",
            "dprst_vol_clos_frac",
            "dprst_vol_open_frac",
            "dprst_area_clos",
            "dprst_area_open",
            "dprst_area_clos_max",
            "dprst_area_open_max",
            "dprst_sroff_hru",
            "dprst_seep_hru",
            "dprst_evap_hru",
            "dprst_insroff_hru",
            "dprst_stor_hru",
            "dprst_stor_hru_old",
            "dprst_stor_hru_change",
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
            "sroff": zero,
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
            "dprst_seep_hru": zero,
            "dprst_evap_hru": zero,
            "dprst_insroff_hru": zero,
            "dprst_stor_hru": zero,
            "dprst_stor_hru_old": zero,
            "dprst_stor_hru_change": zero,
        }

    def _advance_variables(self) -> None:
        """Advance the variables
        Returns:
            None
        """
        self.hru_impervstor_old[:] = self.hru_impervstor
        self.dprst_stor_hru_old[:] = self.dprst_stor_hru
        return None

    def calculate(self, time_length, vectorized=False):
        """Calculate terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """
        self.calculate_prms_style()

        self.hru_impervstor_change[:] = (
            self.hru_impervstor_old - self.hru_impervstor
        )
        self.dprst_stor_hru_change[:] = (
            self.dprst_stor_hru_old - self.dprst_stor_hru
        )

        return

    def calculate_prms_style(self):

        dprst_chk = 0
        self.infil[:] = 0.0
        # DO k = 1, Active_hrus
        for k in range(self.nhru):

            # cdl i = Hru_route_order(k)
            i = k
            ihru = i

            hruarea = self.hru_area[i]
            hruarea_dble = hruarea
            upslope = 0.0

            # cdl if cascade_flag > CASCADE_OFF:
            if False:
                upslope = None  # Upslope_hortonian(i)
            ihru = i
            runoff = 0.0
            glcrmltb = 0.0
            isglacier = 0
            active_glacier = -1

            # cdl -- following glacier code is not implemented
            # IF ( Glacier_flag==ACTIVE ) THEN
            #   IF ( Hru_type(i)==GLACIER ) THEN
            #     IF ( Glacier_flag==ACTIVE ) THEN ! glacier
            #       Isglacier = 1
            #       glcrmltb = Glacrb_melt(i)
            #       IF ( Glacier_frac(i)>0.0 ) THEN
            #         active_glacier = ACTIVE
            #       ELSE
            #         active_glacier = OFF  ! glacier capable HRU, but not glaciated
            #       ENDIF
            #     ENDIF
            #   ENDIF
            # ENDIF

            # cdl -- following lake code is not implemented
            #         IF ( Hru_type(i)==LAKE ) THEN
            # ! HRU is a lake
            # !     eventually add code for lake area less than hru_area
            # !     that includes soil_moist for fraction of hru_area that is dry bank
            #           Hortonian_lakes(i) = upslope
            #           Basin_hortonian_lakes = Basin_hortonian_lakes + Hortonian_lakes(i)*Hruarea_dble
            #           CYCLE
            #         ENDIF

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

            # cdl -- the following frozen code is not implemented
            # frzen = OFF
            # IF ( Frozen_flag==ACTIVE ) THEN
            #   IF ( Tavgc(i)>0.0 ) THEN
            #     cfgi_k = 0.5
            #   ELSE
            #     cfgi_k = 0.08
            #   ENDIF
            #   depth_cm = SNGL(Pk_depth(i))*2.54 !depth of snow cover averaged over HRU
            #   Cfgi(i) = Cfgi_decay*Cfgi_prev(i) - Tavgc(i)*( 2.71828**(-0.4*cfgi_k*depth_cm) )
            #   IF ( active_glacier==1 ) THEN
            #     Cfgi(i) = 0.0 !if glacier over, want ground completely unfrozen, or below threshold, infiltration
            #     IF ( Glacier_frac(i)<1.0 ) Cfgi(i) = Cfgi_thrshld ! glacier with some open fraction
            #   ENDIF
            #   IF ( Cfgi(i)<0.0 ) Cfgi(i) = 0.0
            #   Cfgi_prev(i) = Cfgi(i)
            #   IF ( Cfgi(i)>=Cfgi_thrshld ) THEN
            #     frzen = 1
            #     ! depression storage states are not changed if frozen
            #     cfgi_sroff = (Snowmelt(i) + availh2o + SNGL(upslope) + glcrmltb)*Hruarea
            #     IF ( Use_transfer_intcp==ACTIVE ) cfgi_sroff = cfgi_sroff + Net_apply(i)*Hruarea
            #     runoff = runoff + cfgi_sroff
            #     Basin_cfgi_sroff = Basin_cfgi_sroff + cfgi_sroff
            #   ENDIF
            #   Frozen(i) = frzen
            # ENDIF

            # cdl -- much of the following code is not implemented, except for the compute_infil call
            # !******Compute runoff for pervious, impervious, and depression storage area, only if not frozen ground
            #         IF ( frzen==OFF ) THEN
            # ! DO IRRIGATION APPLICATION, ONLY DONE HERE, ASSUMES NO SNOW and
            # ! only for pervious areas (just like infiltration)
            #           IF ( Use_transfer_intcp==ACTIVE ) THEN
            #             IF ( Net_apply(i)>0.0 ) THEN
            #               sra = 0.0
            #               Infil(i) = Infil(i) + Net_apply(i)
            #               CALL perv_comp(Net_apply(i), Net_apply(i), Infil(i), sra)
            # ! ** ADD in water from irrigation application and water-use transfer for pervious portion - sra (if any)
            #               apply_sroff = DBLE( sra*perv_area )
            #               Basin_apply_sroff = Basin_apply_sroff + apply_sroff
            #               runoff = runoff + apply_sroff
            #             ENDIF
            #           ENDIF
            #
            #           IF ( Isglacier==OFF ) THEN
            #             CALL compute_infil(Net_rain(i), Net_ppt(i), Imperv_stor(i), Imperv_stor_max(i), Snowmelt(i), &
            #      &                         Snowinfil_max(i), Net_snow(i), Pkwater_equiv(i), Infil(i), Hru_type(i), Intcp_changeover(i))
            sri, srp, self.imperv_stor[i], self.infil[i] = self.compute_infil(
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
                ihru,
            )
            #           ELSE ! glacier
            #             temp = Snowmelt(i) + glcrmltb !Snowmelt or 0.0
            #             temp2 = availh2o*(1.0-Glacier_frac(i))
            #             CALL compute_infil(temp2, Net_ppt(i), Imperv_stor(i), Imperv_stor_max(i), temp, &
            #      &                         Snowinfil_max(i), Net_snow(i), Pkwater_equiv(i), Infil(i), Hru_type(i), Intcp_changeover(i))
            #           ENDIF
            #         ENDIF

            #         IF ( Dprst_flag==ACTIVE ) THEN
            #           Dprst_in(i) = 0.0D0
            #           dprst_chk = OFF
            #           IF ( Dprst_area_max(i)>0.0 ) THEN
            #             dprst_chk = ACTIVE
            # !           ******Compute the depression storage component
            # !           only call if total depression surface area for each HRU is > 0.0
            #             IF ( frzen==OFF ) THEN
            #               CALL dprst_comp(Dprst_vol_clos(i), Dprst_area_clos_max(i), Dprst_area_clos(i), &
            #      &                        Dprst_vol_open_max(i), Dprst_vol_open(i), Dprst_area_open_max(i), Dprst_area_open(i), &
            #      &                        Dprst_sroff_hru(i), Dprst_seep_hru(i), Sro_to_dprst_perv(i), Sro_to_dprst_imperv(i), &
            #      &                        Dprst_evap_hru(i), avail_et, availh2o, Dprst_in(i))
            #               runoff = runoff + Dprst_sroff_hru(i)*Hruarea_dble
            #             ENDIF
            #           ENDIF
            #         ENDIF
            dprst_flag = ACTIVE  # cdl todo: hardwired
            frzen = OFF  # cdl todo: hardwired
            if dprst_flag == ACTIVE:
                self.dprst_in[i] = 0.0
                dprst_chk = OFF
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

            # cdl -- the upper part of this block needs to be done to calculate runoff and srunoff
            #         srunoff = 0.0
            #         IF ( Hru_type(i)==LAND .OR. active_glacier==OFF ) THEN ! could be an glacier-capable HRU with no ice
            # !******Compute runoff for pervious and impervious area, and depression storage area
            #           runoff = runoff + DBLE( Srp*perv_area + Sri*Hruarea_imperv )
            #           srunoff = SNGL( runoff/Hruarea_dble )
            srunoff = 0.0
            if self.hru_type[i] == LAND:
                runoff = runoff + srp * perv_area + sri * hruarea_imperv
                srunoff = runoff / hruarea
                #
                # !******Compute HRU weighted average (to units of inches/dt)
                #           IF ( Cascade_flag>CASCADE_OFF ) THEN
                #             hru_sroff_down = 0.0D0
                #             IF ( srunoff>0.0 ) THEN
                #               IF ( Ncascade_hru(i)>0 ) CALL run_cascade_sroff(Ncascade_hru(i), srunoff, hru_sroff_down)
                #               Hru_hortn_cascflow(i) = hru_sroff_down
                #               !IF ( Hru_hortn_cascflow(i)<0.0D0 ) Hru_hortn_cascflow(i) = 0.0D0
                #               !IF ( Upslope_hortonian(i)<0.0D0 ) Upslope_hortonian(i) = 0.0D0
                #               Basin_sroff_upslope = Basin_sroff_upslope + Upslope_hortonian(i)*Hruarea_dble
                #               Basin_sroff_down = Basin_sroff_down + hru_sroff_down*Hruarea_dble
                #             ELSE
                #               Hru_hortn_cascflow(i) = 0.0D0
                #             ENDIF
                #           ENDIF
                self.hru_sroffp[i] = srp * perv_frac
            #           Hru_sroffp(i) = Srp*Perv_frac
            #           Basin_sroffp = Basin_sroffp + Srp*perv_area
            #         ENDIF

            # Basin_infil = Basin_infil + DBLE(Infil(i) * perv_area)
            # Basin_contrib_fraction = Basin_contrib_fraction + DBLE(Contrib_fraction(i) * perv_area)

            # cdl -- the guts of this was implemented in the python code below
            # !******Compute evaporation from impervious area
            #         IF ( frzen==OFF ) THEN
            #         IF ( Hruarea_imperv>0.0 ) THEN
            #           IF ( Imperv_stor(i)>0.0 ) THEN
            #             CALL imperv_et(Imperv_stor(i), Potet(i), Imperv_evap(i), Snowcov_area(i), avail_et)
            #             Hru_impervevap(i) = Imperv_evap(i)*Imperv_frac
            #             !IF ( Hru_impervevap(i)<0.0 ) Hru_impervevap(i) = 0.0
            #             avail_et = avail_et - Hru_impervevap(i)
            #             IF ( avail_et<0.0 ) THEN
            #                ! sanity check
            # !              IF ( avail_et<-NEARZERO ) PRINT*, 'avail_et<0 in srunoff imperv', i, Nowmonth, Nowday, avail_et
            #               Hru_impervevap(i) = Hru_impervevap(i) + avail_et
            #               IF ( Hru_impervevap(i)<0.0 ) Hru_impervevap(i) = 0.0
            #               Imperv_evap(i) = Hru_impervevap(i)/Imperv_frac
            #               Imperv_stor(i) = Imperv_stor(i) - avail_et/Imperv_frac
            #               avail_et = 0.0
            #             ENDIF
            #             Basin_imperv_evap = Basin_imperv_evap + DBLE( Hru_impervevap(i)*Hruarea )
            #             Hru_impervstor(i) = Imperv_stor(i)*Imperv_frac
            #             Basin_imperv_stor = Basin_imperv_stor + DBLE(Imperv_stor(i)*Hruarea_imperv )
            #           ENDIF
            #           Hru_sroffi(i) = Sri*Imperv_frac
            #           Basin_sroffi = Basin_sroffi + DBLE( Sri*Hruarea_imperv )
            #         ENDIF
            #         ENDIF
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

            # cdl -- saving sroff here
            #         IF ( dprst_chk==ACTIVE ) Dprst_stor_hru(i) = (Dprst_vol_open(i)+Dprst_vol_clos(i))/Hruarea_dble
            if dprst_chk == ACTIVE:
                self.dprst_stor_hru[i] = (
                    self.dprst_vol_open[i] + self.dprst_vol_clos[i]
                ) / hruarea
            #
            #         Sroff(i) = srunoff
            self.sroff[i] = srunoff
            #         Hortonian_flow(i) = srunoff
            #         Basin_hortonian = Basin_hortonian + DBLE( srunoff*Hruarea )
            #         Basin_sroff = Basin_sroff + DBLE( srunoff*Hruarea )
            #       ENDDO
            # cdl -- this is the end of the hru do loop

            # cdl -- nothing to do here
            # !******Compute basin weighted averages (to units of inches/dt)
            #       !rsr, should be land_area???
            #       Basin_sroff = Basin_sroff*Basin_area_inv
            #       Basin_imperv_evap = Basin_imperv_evap*Basin_area_inv
            #       Basin_imperv_stor = Basin_imperv_stor*Basin_area_inv
            #       Basin_infil = Basin_infil*Basin_area_inv
            #       ! doesn't include CFGI runoff
            #       Basin_sroffp = Basin_sroffp*Basin_area_inv
            #       Basin_sroffi = Basin_sroffi*Basin_area_inv
            #       Basin_hortonian = Basin_hortonian*Basin_area_inv
            #       Basin_contrib_fraction = Basin_contrib_fraction*Basin_area_inv
            #       Basin_apply_sroff = Basin_apply_sroff*Basin_area_inv
            #       IF ( Cascade_flag>CASCADE_OFF ) THEN
            #         Basin_hortonian_lakes = Basin_hortonian_lakes*Basin_area_inv
            #         Basin_sroff_down = Basin_sroff_down*Basin_area_inv
            #         Basin_sroff_upslope = Basin_sroff_upslope*Basin_area_inv
            #       ENDIF
            #
            #       IF ( Dprst_flag==ACTIVE ) THEN
            #         Basin_dprst_volop = Basin_dprst_volop*Basin_area_inv
            #         Basin_dprst_volcl = Basin_dprst_volcl*Basin_area_inv
            #         Basin_dprst_evap = Basin_dprst_evap*Basin_area_inv
            #         Basin_dprst_seep = Basin_dprst_seep*Basin_area_inv
            #         Basin_dprst_sroff = Basin_dprst_sroff*Basin_area_inv
            #       ENDIF

        return

    def basin_init(self):
        """
        This is trying to replicate the prms basin_init function that calculates some of the
        variables needed here by runoff.  This should probably go somewhere else at some point
        as I suspect other components may need similar information.

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
                        dprst_clos_flag = ACTIVE  # cdl -- todo: this variable should be stored to self?
                    if self.dprst_area_open_max[i] > 0.0:
                        dprst_open_flag = ACTIVE  # cdl -- todo: this variable should be stored to self?
                    perv_area = perv_area - self.dprst_area_max[i]

            self.hru_perv[i] = perv_area
            self.hru_frac_perv[i] = perv_area / harea
        return

    #       SUBROUTINE imperv_et(Imperv_stor, Potet, Imperv_evap, Sca, Avail_et)
    #       USE PRMS_SRUNOFF, ONLY: Imperv_frac
    #       IMPLICIT NONE
    # ! Arguments
    #       REAL, INTENT(IN) :: Potet, Sca, Avail_et
    #       REAL, INTENT(INOUT) :: Imperv_stor, Imperv_evap
    # !***********************************************************************
    #       IF ( Sca<1.0 ) THEN
    #         IF ( Potet<Imperv_stor ) THEN
    #           Imperv_evap = Potet*(1.0-Sca)
    #         ELSE
    #           Imperv_evap = Imperv_stor*(1.0-Sca)
    #         ENDIF
    #         IF ( Imperv_evap*Imperv_frac>Avail_et ) Imperv_evap = Avail_et/Imperv_frac
    #         Imperv_stor = Imperv_stor - Imperv_evap
    #       ENDIF
    #       END SUBROUTINE imperv_et
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

    #       SUBROUTINE compute_infil(Net_rain, Net_ppt, Imperv_stor, Imperv_stor_max, Snowmelt, &
    #      &                         Snowinfil_max, Net_snow, Pkwater_equiv, Infil, Hru_type, Intcp_changeover)
    #       USE PRMS_CONSTANTS, ONLY: NEARZERO, DNEARZERO, LAND, ACTIVE, CASCADE_OFF
    #       USE PRMS_MODULE, ONLY: Cascade_flag
    #       USE PRMS_SRUNOFF, ONLY: Sri, Hruarea_imperv, Upslope_hortonian, Ihru, Srp, Isglacier
    #       USE PRMS_SNOW, ONLY: Pptmix_nopack
    #       IMPLICIT NONE
    # ! Arguments
    #       INTEGER, INTENT(IN) :: Hru_type
    #       REAL, INTENT(IN) :: Net_rain, Net_ppt, Imperv_stor_max
    #       REAL, INTENT(IN) :: Snowmelt, Snowinfil_max, Net_snow, Intcp_changeover
    #       DOUBLE PRECISION, INTENT(IN) :: Pkwater_equiv
    #       REAL, INTENT(INOUT) :: Imperv_stor, Infil
    def compute_infil(
        self,
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
        ihru,
    ):

        # cdl -- did not implement cascade
        #       hru_flag = 0
        #       IF ( Hru_type==LAND .OR. Isglacier==ACTIVE ) hru_flag = 1 ! land or glacier
        # ! compute runoff from cascading Hortonian flow
        #       IF ( Cascade_flag>CASCADE_OFF ) THEN
        #         avail_water = SNGL( Upslope_hortonian(Ihru) )
        #         IF ( avail_water>0.0 ) THEN
        #           Infil = Infil + avail_water
        #           IF ( hru_flag==1 ) CALL perv_comp(avail_water, avail_water, Infil, Srp)
        #         ENDIF
        #       ELSE
        #         avail_water = 0.0
        #       ENDIF
        isglacier = False  # todo -- hardwired
        cascade_active = False
        hru_flag = 0
        if hru_type == LAND or isglacier:
            hru_flag = 1
        if cascade_active:
            raise Exception("bad bad bad")
        else:
            avail_water = 0.0

        # ! compute runoff from canopy changeover water
        #       IF ( Intcp_changeover>0.0 ) THEN
        #         avail_water = avail_water + Intcp_changeover
        #         Infil = Infil + Intcp_changeover
        #         IF ( hru_flag==1 ) CALL perv_comp(Intcp_changeover, Intcp_changeover, Infil, Srp)
        #       ENDIF
        if intcp_changeover > 0.0:
            avail_water = avail_water + intcp_changeover
            infil = infil + intcp_changeover
            if hru_flag == 1:
                infil, srp = self.perv_comp(
                    intcp_changeover, intcp_changeover, infil, srp, ihru
                )

        # !******if rain/snow event with no antecedent snowpack,
        # !******compute the runoff from the rain first and then proceed with the
        # !******snowmelt computations
        #
        #       IF ( Pptmix_nopack(Ihru)==ACTIVE ) THEN
        #         avail_water = avail_water + Net_rain
        #         Infil = Infil + Net_rain
        #         IF ( hru_flag==1 ) CALL perv_comp(Net_rain, Net_rain, Infil, Srp)
        #       ENDIF
        if self.pptmix_nopack[ihru] == ACTIVE:
            avail_water = avail_water + net_rain
            infil = infil + net_rain
            if hru_flag == 1:
                infil, srp = self.perv_comp(
                    net_rain, net_rain, infil, srp, ihru
                )

        # !******If precipitation on snowpack, all water available to the surface is
        # !******considered to be snowmelt, and the snowmelt infiltration
        # !******procedure is used.  If there is no snowpack and no precip,
        # !******then check for melt from last of snowpack.  If rain/snow mix
        # !******with no antecedent snowpack, compute snowmelt portion of runoff.
        #       IF ( Snowmelt>0.0 ) THEN ! includes glacier melt, if any
        #         avail_water = avail_water + Snowmelt
        #         Infil = Infil + Snowmelt
        #         IF ( hru_flag==1 ) THEN
        #           IF ( Pkwater_equiv>0.0D0 .OR. Net_ppt-Net_snow<NEARZERO ) THEN
        # !******Pervious area computations
        #             CALL check_capacity(Snowinfil_max, Infil)
        # !******Snowmelt occurred and depleted the snowpack
        #           ELSE
        #             CALL perv_comp(Snowmelt, Net_ppt, Infil, Srp)
        #           ENDIF
        #         ENDIF
        if snowmelt > 0.0:
            avail_water = avail_water + snowmelt
            infil = infil + snowmelt
            if hru_flag == 1:
                if pkwater_equiv > 0.0 or net_ppt - net_snow < NEARZERO:
                    infil, srp = self.check_capacity(
                        snowinfil_max, infil, srp, ihru
                    )
                else:
                    infil, srp = self.perv_comp(
                        snowmelt, net_ppt, infil, srp, ihru
                    )

        # !******There was no snowmelt but a snowpack may exist.  If there is
        # !******no snowpack then check for rain on a snowfree HRU.
        #
        #       ELSEIF ( Pkwater_equiv<DNEARZERO ) THEN
        #
        # !       If no snowmelt and no snowpack but there was net snow then
        # !       snowpack was small and was lost to sublimation.
        #
        #         IF ( Net_snow<NEARZERO .AND. Net_rain>0.0 ) THEN
        # ! no snow, some rain
        #           avail_water = avail_water + Net_rain
        #           Infil = Infil + Net_rain
        #           IF ( hru_flag==1 ) CALL perv_comp(Net_rain, Net_rain, Infil, Srp)
        #         ENDIF
        elif pkwater_equiv < DNEARZERO:
            if net_snow < NEARZERO and net_rain > 0.0:
                avail_water = avail_water + net_rain
                infil = infil + net_rain
                if hru_flag == 1:
                    infil, srp = self.perv_comp(
                        net_rain, net_rain, infil, srp, ihru
                    )
        #
        # !***** Snowpack exists, check to see if infil exceeds maximum daily
        # !***** snowmelt infiltration rate. Infil results from rain snow mix
        # !***** on a snowfree surface.
        #
        #       ELSEIF ( Infil>0.0 ) THEN
        #         IF ( hru_flag==1 ) CALL check_capacity(Snowinfil_max, Infil)
        #       ENDIF
        elif infil > 0.0:
            if hru_flag == 1:
                infil, srp = self.check_capacity(
                    snowinfil_max, infil, srp, ihru
                )

        # !******Impervious area computations
        #       IF ( Hruarea_imperv>0.0 ) THEN
        #         Imperv_stor = Imperv_stor + avail_water
        #         IF ( hru_flag==1 ) THEN
        #           IF ( Imperv_stor>Imperv_stor_max ) THEN
        #             Sri = Imperv_stor - Imperv_stor_max
        #             Imperv_stor = Imperv_stor_max
        #           ENDIF
        #         ENDIF
        #       ENDIF
        if hruarea_imperv > 0.0:
            imperv_stor = imperv_stor + avail_water
            if hru_flag == 1:
                if imperv_stor > imperv_stor_max:
                    sri = imperv_stor - imperv_stor_max
                    imperv_stor = imperv_stor_max

        return sri, srp, imperv_stor, infil

    #       SUBROUTINE perv_comp(Pptp, Ptc, Infil, Srp)
    def perv_comp(self, pptp, ptc, infil, srp, ihru):
        #       USE PRMS_CONSTANTS, ONLY: smidx_module !, CLOSEZERO
        #       USE PRMS_MODULE, ONLY: Sroff_flag
        #       USE PRMS_SRUNOFF, ONLY: Ihru, Smidx_coef, Smidx_exp, &
        #      &    Carea_max, Carea_min, Carea_dif, Contrib_fraction
        #       USE PRMS_FLOWVARS, ONLY: Soil_moist, Soil_rechr, Soil_rechr_max
        #       IMPLICIT NONE
        # ! Arguments
        #       REAL, INTENT(IN) :: Pptp, Ptc
        #       REAL, INTENT(INOUT) :: Infil, Srp
        # ! Local Variables
        #       REAL :: smidx, srpp, ca_fraction
        # !***********************************************************************
        # !******Pervious area computations
        #       IF ( Sroff_flag==smidx_module ) THEN
        #         ! antecedent soil_moist
        #         smidx = Soil_moist(Ihru) + (0.5*Ptc)
        #         IF ( smidx>25.0) THEN
        #           ca_fraction = Carea_max(Ihru)
        #         ELSE
        #           ca_fraction = Smidx_coef(Ihru)*10.0**(Smidx_exp(Ihru)*smidx)
        #         ENDIF
        #       ELSE
        #         ! antecedent soil_rechr
        #         ca_fraction = Carea_min(Ihru) + Carea_dif(Ihru)*(Soil_rechr(Ihru)/Soil_rechr_max(Ihru))
        #       ENDIF
        smidx_module = True
        if smidx_module:
            smidx = self.soil_moist_prev[ihru] + 0.5 * ptc
            if smidx > 25.0:
                ca_fraction = self.carea_max[ihru]
            else:
                ca_fraction = self.smidx_coef[ihru] * 10.0 ** (
                    self.smidx_exp[ihru] * smidx
                )
        else:
            raise Exception("you did a bad thing...")

        #       IF ( ca_fraction>Carea_max(Ihru) ) ca_fraction = Carea_max(Ihru)
        #       srpp = ca_fraction*Pptp
        #       Contrib_fraction(Ihru) = ca_fraction
        # !      IF ( srpp<0.0 ) THEN
        # !        PRINT *, 'negative srp', srpp
        # !        srpp = 0.0
        # !      ENDIF
        #       Infil = Infil - srpp
        #       Srp = Srp + srpp
        #       !IF ( Srp<CLOSEZERO ) Srp = 0.0
        #
        #       END SUBROUTINE perv_comp
        if ca_fraction > self.carea_max[ihru]:
            ca_fraction = self.carea_max[ihru]
        srpp = ca_fraction * pptp
        self.contrib_fraction[ihru] = ca_fraction
        infil = infil - srpp
        srp = srp + srpp

        return infil, srp

    # !***********************************************************************
    # ! fill soil to soil_moist_max, if more than capacity restrict
    # ! infiltration by snowinfil_max, with excess added to runoff
    # !***********************************************************************
    #       SUBROUTINE check_capacity(Snowinfil_max, Infil)
    #       USE PRMS_FLOWVARS, ONLY: Soil_moist_max, Soil_moist
    #       USE PRMS_SRUNOFF, ONLY: Ihru, Srp
    #       IMPLICIT NONE
    # ! Arguments
    #       REAL, INTENT(IN) :: Snowinfil_max
    #       REAL, INTENT(INOUT) :: Infil
    # ! Local Variables
    #       REAL :: capacity, excess
    # !***********************************************************************
    #       capacity = Soil_moist_max(Ihru) - Soil_moist(Ihru)
    #       excess = Infil - capacity
    #       IF ( excess>Snowinfil_max ) THEN
    #         Srp = Srp + excess - Snowinfil_max
    #         Infil = Snowinfil_max + capacity
    #       ENDIF
    #
    #       END SUBROUTINE check_capacity
    def check_capacity(self, snowinfil_max, infil, srp, ihru):
        capacity = self.soil_moist_max[ihru] - self.soil_moist_prev[ihru]
        excess = infil - capacity
        if excess > snowinfil_max:
            srp = srp + excess - snowinfil_max
            infil = snowinfil_max + capacity
        return infil, srp

    #       SUBROUTINE dprst_init()
    #       USE PRMS_SRUNOFF
    #       USE PRMS_CONSTANTS, ONLY: ACTIVE, NEARZERO
    #       USE PRMS_MODULE, ONLY: Init_vars_from_file, Nhru, PRMS4_flag, Inputerror_flag
    #       USE PRMS_BASIN, ONLY: Dprst_clos_flag, Dprst_frac, &
    #      &    Dprst_area_clos_max, Dprst_area_open_max, Basin_area_inv, &
    #      &    Hru_area_dble, Active_hrus, Hru_route_order, Dprst_open_flag
    #       USE PRMS_FLOWVARS, ONLY: Dprst_vol_open, Dprst_vol_clos
    #       IMPLICIT NONE
    # ! Functions
    #       INTRINSIC :: EXP, LOG, DBLE, SNGL
    #       INTEGER, EXTERNAL :: getparam
    # ! Local Variables
    #       INTEGER :: i, j
    #       REAL :: frac_op_ar, frac_cl_ar, open_vol_r, clos_vol_r
    # !***********************************************************************
    def dprst_init(self):

        #       Dprst_evap_hru = 0.0
        #       Dprst_seep_hru = 0.0D0
        #       Dprst_sroff_hru = 0.0D0
        #       Dprst_insroff_hru = 0.0
        #       IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==7 ) THEN
        #         IF ( getparam(MODNAME, 'dprst_frac_init', Nhru, 'real', Dprst_frac_init)/=0 ) CALL read_error(2, 'dprst_frac_init')
        #       ENDIF
        #       IF ( getparam(MODNAME, 'dprst_flow_coef', Nhru, 'real', Dprst_flow_coef)/=0 ) CALL read_error(2, 'dprst_flow_coef')
        #       IF ( Dprst_open_flag==ACTIVE ) THEN
        #         IF ( getparam(MODNAME, 'dprst_seep_rate_open', Nhru, 'real', Dprst_seep_rate_open)/=0 )  &
        #      &       CALL read_error(2, 'dprst_seep_rate_open')
        #         IF ( getparam(MODNAME, 'va_open_exp', Nhru, 'real', Va_open_exp)/=0 ) CALL read_error(2, 'va_open_exp')
        #         IF ( getparam(MODNAME, 'op_flow_thres', Nhru, 'real', Op_flow_thres)/=0 ) CALL read_error(2, 'op_flow_thres')
        #       ELSE
        #         Dprst_seep_rate_open = 0.0
        #         Va_open_exp = 0.0
        #         Op_flow_thres = 0.0
        #       ENDIF
        #       IF ( PRMS4_flag==ACTIVE ) THEN
        #         IF ( getparam(MODNAME, 'sro_to_dprst', Nhru, 'real', Sro_to_dprst_perv)/=0 ) CALL read_error(2, 'sro_to_dprst')
        #       ELSE
        #         IF ( getparam(MODNAME, 'sro_to_dprst_perv', Nhru, 'real', Sro_to_dprst_perv)/=0 ) CALL read_error(2, 'sro_to_dprst_perv')
        #       ENDIF
        #       IF ( getparam(MODNAME, 'sro_to_dprst_imperv', Nhru, 'real', Sro_to_dprst_imperv)/=0 ) &
        #      &     CALL read_error(2, 'sro_to_dprst_imperv')
        #       IF ( getparam(MODNAME, 'dprst_depth_avg', Nhru, 'real', Dprst_depth_avg)/=0 ) CALL read_error(2, 'dprst_depth_avg')
        #       IF ( getparam(MODNAME, 'dprst_et_coef', Nhru, 'real', Dprst_et_coef)/=0 ) CALL read_error(2, 'dprst_et_coef')
        #       IF ( Dprst_clos_flag==ACTIVE ) THEN
        #         IF ( getparam(MODNAME, 'dprst_seep_rate_clos', Nhru, 'real', Dprst_seep_rate_clos)/=0 ) &
        #      &       CALL read_error(2, 'dprst_seep_rate_clos')
        #         IF ( getparam(MODNAME, 'va_clos_exp', Nhru, 'real', Va_clos_exp)/=0 ) CALL read_error(2, 'va_clos_exp')
        #       ELSE
        #         Dprst_seep_rate_clos = 0.0
        #         Va_clos_exp = 0.0
        #       ENDIF
        #       Dprst_in = 0.0D0
        #       Dprst_area_open = 0.0
        #       Dprst_area_clos = 0.0
        #       Dprst_stor_hru = 0.0D0
        #       Dprst_vol_thres_open = 0.0D0
        #       Dprst_vol_open_max = 0.0D0
        #       Dprst_vol_clos_max = 0.0D0
        #       Dprst_vol_frac = 0.0
        #       Dprst_vol_open_frac = 0.0
        #       Dprst_vol_clos_frac = 0.0
        #       Basin_dprst_volop = 0.0D0
        #       Basin_dprst_volcl = 0.0D0
        #       DO j = 1, Active_hrus
        #         i = Hru_route_order(j)
        for j in range(self.nhru):
            i = j
            #
            #         IF ( Dprst_frac(i)>0.0 ) THEN
            if self.dprst_frac[i] > 0.0:

                #           IF ( Dprst_depth_avg(i)==0.0 ) THEN
                #             PRINT *, 'ERROR, dprst_frac>0 and dprst_depth_avg==0 for HRU:', i, '; dprst_frac:', Dprst_frac(i)
                #             Inputerror_flag = 1
                #             CYCLE
                #           ENDIF
                if self.dprst_depth_avg[i] == 0.0:
                    raise Exception(
                        f"dprst_fac > and dprst_depth_avg == 0 for HRU {i}"
                    )

                # !         calculate open and closed volumes (acre-inches) of depression storage by HRU
                # !         Dprst_area_open_max is the maximum open depression area (acres) that can generate surface runoff:
                #           IF ( Dprst_clos_flag==ACTIVE ) Dprst_vol_clos_max(i) = DBLE( Dprst_area_clos_max(i)*Dprst_depth_avg(i) )
                #           IF ( Dprst_open_flag==ACTIVE ) Dprst_vol_open_max(i) = DBLE( Dprst_area_open_max(i)*Dprst_depth_avg(i) )
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

                #
                # !         calculate the initial open and closed depression storage volume:
                #           IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==7 ) THEN
                #             IF ( Dprst_open_flag==ACTIVE ) Dprst_vol_open(i) = DBLE(Dprst_frac_init(i))*Dprst_vol_open_max(i)
                #             IF ( Dprst_clos_flag==ACTIVE ) Dprst_vol_clos(i) = DBLE(Dprst_frac_init(i))*Dprst_vol_clos_max(i)
                #           ENDIF
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

                #
                # !         threshold volume is calculated as the % of maximum open
                # !         depression storage above which flow occurs *  total open depression storage volume
                #           Dprst_vol_thres_open(i) = DBLE(Op_flow_thres(i))*Dprst_vol_open_max(i)
                self.dprst_vol_thres_open[i] = (
                    self.op_flow_thres[i] * self.dprst_vol_open_max[i]
                )
                #
                # !         initial open and closed storage volume as fraction of total open and closed storage volume
                #
                # !         Open depression surface area for each HRU:
                #           IF ( Dprst_vol_open(i)>0.0D0 ) THEN
                #             open_vol_r = SNGL( Dprst_vol_open(i)/Dprst_vol_open_max(i) )
                #             IF ( open_vol_r<NEARZERO ) THEN
                #               frac_op_ar = 0.0
                #             ELSEIF ( open_vol_r>1.0 ) THEN
                #               frac_op_ar = 1.0
                #             ELSE
                #               frac_op_ar = EXP(Va_open_exp(i)*LOG(open_vol_r))
                #             ENDIF
                #             Dprst_area_open(i) = Dprst_area_open_max(i)*frac_op_ar
                #             IF ( Dprst_area_open(i)>Dprst_area_open_max(i) ) Dprst_area_open(i) = Dprst_area_open_max(i)
                # !            IF ( Dprst_area_open(i)<NEARZERO ) Dprst_area_open(i) = 0.0
                #           ENDIF
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

                #
                # !         Closed depression surface area for each HRU:
                #           IF ( Dprst_vol_clos(i)>0.0D0 ) THEN
                #             clos_vol_r = SNGL( Dprst_vol_clos(i)/Dprst_vol_clos_max(i) )
                #             IF ( clos_vol_r<NEARZERO ) THEN
                #               frac_cl_ar = 0.0
                #             ELSEIF ( clos_vol_r>1.0 ) THEN
                #               frac_cl_ar = 1.0
                #             ELSE
                #               frac_cl_ar = EXP(Va_clos_exp(i)*LOG(clos_vol_r))
                #             ENDIF
                #             Dprst_area_clos(i) = Dprst_area_clos_max(i)*frac_cl_ar
                #             IF ( Dprst_area_clos(i)>Dprst_area_clos_max(i) ) Dprst_area_clos(i) = Dprst_area_clos_max(i)
                # !            IF ( Dprst_area_clos(i)<NEARZERO ) Dprst_area_clos(i) = 0.0
                #           ENDIF
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
                # !         calculate the basin open and closed depression storage volumes
                #           Basin_dprst_volop = Basin_dprst_volop + Dprst_vol_open(i)
                #           Basin_dprst_volcl = Basin_dprst_volcl + Dprst_vol_clos(i)
                #           Dprst_stor_hru(i) = (Dprst_vol_open(i)+Dprst_vol_clos(i))/Hru_area_dble(i)
                self.dprst_stor_hru[i] = (
                    self.dprst_vol_open[i] + self.dprst_vol_clos[i]
                ) / self.hru_area[i]
                self.dprst_stor_hru_old[i] = self.dprst_stor_hru[i]

                #           IF ( Dprst_vol_open_max(i)>0.0 ) Dprst_vol_open_frac(i) = SNGL( Dprst_vol_open(i)/Dprst_vol_open_max(i) )
                if self.dprst_vol_open_max[i] > 0.0:
                    self.dprst_vol_open_frac[i] = (
                        self.dprst_vol_open[i] / self.dprst_vol_open_max[i]
                    )

                #           IF ( Dprst_vol_clos_max(i)>0.0 ) Dprst_vol_clos_frac(i) = SNGL( Dprst_vol_clos(i)/Dprst_vol_clos_max(i) )
                if self.dprst_vol_clos_max[i] > 0.0:
                    self.dprst_vol_clos_frac[i] = (
                        self.dprst_vol_clos[i] / self.dprst_vol_clos_max[i]
                    )

                #           if (Dprst_vol_open_max(i)+Dprst_vol_clos_max(i) > 0.0 ) then  !RGN added to avoid divide by zero 1/20/2022
                #             Dprst_vol_frac(i) = SNGL( (Dprst_vol_open(i)+Dprst_vol_clos(i))/(Dprst_vol_open_max(i)+Dprst_vol_clos_max(i)) )
                #           end if
                if (
                    self.dprst_vol_open_max[i] + self.dprst_vol_clos_max[i]
                    > 0.0
                ):
                    self.dprst_vol_frac[i] = (
                        self.dprst_vol_open[i] + self.dprst_vol_clos[i]
                    ) / (
                        self.dprst_vol_open_max[i] + self.dprst_vol_clos_max[i]
                    )

        #         ENDIF
        #       ENDDO
        #       Basin_dprst_volop = Basin_dprst_volop*Basin_area_inv
        #       Basin_dprst_volcl = Basin_dprst_volcl*Basin_area_inv
        #
        #       IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==7 ) DEALLOCATE ( Dprst_frac_init )
        #
        #       END SUBROUTINE dprst_init
        return

    #       SUBROUTINE dprst_comp(Dprst_vol_clos, Dprst_area_clos_max, Dprst_area_clos, &
    #      &           Dprst_vol_open_max, Dprst_vol_open, Dprst_area_open_max, Dprst_area_open, &
    #      &           Dprst_sroff_hru, Dprst_seep_hru, Sro_to_dprst_perv, Sro_to_dprst_imperv, Dprst_evap_hru, &
    #      &           Avail_et, Net_rain, Dprst_in)
    #       USE PRMS_CONSTANTS, ONLY: ERROR_water_use, NEARZERO, DNEARZERO, OFF, ACTIVE ! , DEBUG_less
    #       USE PRMS_MODULE, ONLY: Cascade_flag, Dprst_add_water_use, Dprst_transfer_water_use, &
    #      &    Nowyear, Nowmonth, Nowday !, Print_debug
    #       USE PRMS_SRUNOFF, ONLY: Srp, Sri, Ihru, Perv_frac, Imperv_frac, Hruarea, Dprst_et_coef, &
    #      &    Dprst_seep_rate_open, Dprst_seep_rate_clos, Va_clos_exp, Va_open_exp, Dprst_flow_coef, &
    #      &    Dprst_vol_thres_open, Dprst_vol_clos_max, Dprst_insroff_hru, Upslope_hortonian, &
    #      &    Basin_dprst_volop, Basin_dprst_volcl, Basin_dprst_evap, Basin_dprst_seep, Basin_dprst_sroff, &
    #      &    Dprst_vol_open_frac, Dprst_vol_clos_frac, Dprst_vol_frac, Dprst_stor_hru, Hruarea_dble
    #       USE PRMS_BASIN, ONLY: Dprst_frac_open, Dprst_frac_clos
    #       USE PRMS_WATER_USE, ONLY: Dprst_transfer, Dprst_gain
    #       USE PRMS_SET_TIME, ONLY: Cfs_conv
    #       USE PRMS_INTCP, ONLY: Net_snow
    #       USE PRMS_CLIMATEVARS, ONLY: Potet
    #       USE PRMS_FLOWVARS, ONLY: Pkwater_equiv
    #       USE PRMS_SNOW, ONLY: Snowmelt, Pptmix_nopack, Snowcov_area
    #       IMPLICIT NONE
    # ! Functions
    #       INTRINSIC :: EXP, LOG, MAX, DBLE, SNGL
    # ! Arguments
    #       REAL, INTENT(IN) :: Dprst_area_open_max, Dprst_area_clos_max, Net_rain
    #       REAL, INTENT(IN) :: Sro_to_dprst_perv, Sro_to_dprst_imperv
    #       DOUBLE PRECISION, INTENT(IN) :: Dprst_vol_open_max
    #       DOUBLE PRECISION, INTENT(INOUT) :: Dprst_vol_open, Dprst_vol_clos, Dprst_in
    #       REAL, INTENT(INOUT) :: Avail_et
    #       REAL, INTENT(OUT) :: Dprst_area_open, Dprst_area_clos, Dprst_evap_hru
    #       DOUBLE PRECISION, INTENT(OUT) :: Dprst_sroff_hru, Dprst_seep_hru
    # ! Local Variables
    #       REAL :: inflow, dprst_avail_et
    #       REAL :: dprst_srp, dprst_sri
    #       REAL :: dprst_srp_open, dprst_srp_clos, dprst_sri_open, dprst_sri_clos
    #       REAL :: frac_op_ar, frac_cl_ar, open_vol_r, clos_vol_r, unsatisfied_et
    #       REAL :: tmp, dprst_evap_open, dprst_evap_clos
    #       DOUBLE PRECISION :: seep_open, seep_clos, tmp1
    # !***********************************************************************
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

        # !     add the hortonian flow to the depression storage volumes:
        #       IF ( Cascade_flag>OFF ) THEN
        #         inflow = SNGL( Upslope_hortonian(Ihru) )
        #       ELSE
        #         inflow = 0.0
        #       ENDIF
        cascade_flag = OFF  # cdl -- todo: hardwired
        if cascade_flag > OFF:
            raise Exception("i am brokin")
        else:
            inflow = 0.0

        #
        #       IF ( Pptmix_nopack(Ihru)==ACTIVE ) inflow = inflow + Net_rain
        if self.pptmix_nopack[ihru]:
            inflow = inflow + net_rain

        #
        # !******If precipitation on snowpack all water available to the surface is considered to be snowmelt
        # !******If there is no snowpack and no precip,then check for melt from last of snowpack.
        # !******If rain/snow mix with no antecedent snowpack, compute snowmelt portion of runoff.
        #
        #       IF ( Snowmelt(Ihru)>0.0 ) THEN
        #         inflow = inflow + Snowmelt(Ihru)
        if self.snowmelt[ihru]:
            inflow = inflow + self.snowmelt[ihru]

        #
        # !******There was no snowmelt but a snowpack may exist.  If there is
        # !******no snowpack then check for rain on a snowfree HRU.
        #       ELSEIF ( Pkwater_equiv(Ihru)<DNEARZERO ) THEN
        elif self.pkwater_equiv[ihru] < DNEARZERO:
            #
            # !      If no snowmelt and no snowpack but there was net snow then
            # !      snowpack was small and was lost to sublimation.
            #         IF ( Net_snow(Ihru)<NEARZERO .AND. Net_rain>0.0 ) THEN
            #           inflow = inflow + Net_rain
            #         ENDIF
            #       ENDIF
            if self.net_snow[ihru] < NEARZERO and net_rain > 0.0:
                inflow = inflow + net_rain

        # cdl -- not implemented
        #
        #       IF ( Dprst_add_water_use==ACTIVE ) THEN
        #         IF ( Dprst_gain(Ihru)>0.0 ) inflow = inflow + Dprst_gain(Ihru) / SNGL( Cfs_conv )
        #       ENDIF
        dprst_add_water_use = OFF
        if dprst_add_water_use == ACTIVE:
            raise Exception("not supported")

        #
        #       Dprst_in = 0.0D0
        #       IF ( Dprst_area_open_max>0.0 ) THEN
        #         Dprst_in = DBLE( inflow*Dprst_area_open_max ) ! inch-acres
        #         Dprst_vol_open = Dprst_vol_open + Dprst_in
        #       ENDIF
        dprst_in = 0.0
        if dprst_area_open_max > 0.0:
            dprst_in = inflow * dprst_area_open_max
            dprst_vol_open = dprst_vol_open + dprst_in

        #
        #       IF ( Dprst_area_clos_max>0.0 ) THEN
        #         tmp1 = DBLE( inflow*Dprst_area_clos_max ) ! inch-acres
        #         Dprst_vol_clos = Dprst_vol_clos + tmp1
        #         Dprst_in = Dprst_in + tmp1
        #       ENDIF
        #       Dprst_in = Dprst_in/Hruarea_dble ! inches over HRU
        if dprst_area_clos_max > 0.0:
            tmp1 = inflow * dprst_area_clos_max
            dprst_vol_clos = dprst_vol_clos + tmp1
            dprst_in = dprst_in + tmp1
        dprst_in = dprst_in / self.hru_area[ihru]

        #
        #       ! add any pervious surface runoff fraction to depressions
        #       dprst_srp = 0.0
        #       dprst_sri = 0.0
        #       IF ( Srp>0.0 ) THEN
        #         tmp = Srp*Perv_frac*Sro_to_dprst_perv*Hruarea
        #         IF ( Dprst_area_open_max>0.0 ) THEN
        #           dprst_srp_open = tmp*Dprst_frac_open(Ihru) ! acre-inches
        #           dprst_srp = dprst_srp_open/Hruarea
        #           Dprst_vol_open = Dprst_vol_open + DBLE( dprst_srp_open )
        #         ENDIF
        #         IF ( Dprst_area_clos_max>0.0 ) THEN
        #           dprst_srp_clos = tmp*Dprst_frac_clos(Ihru)
        #           dprst_srp = dprst_srp + dprst_srp_clos/Hruarea
        #           Dprst_vol_clos = Dprst_vol_clos + DBLE( dprst_srp_clos )
        #         ENDIF
        #         Srp = Srp - dprst_srp/Perv_frac
        #         IF ( Srp<0.0 ) THEN
        #           IF ( Srp<-NEARZERO ) PRINT *, 'dprst srp<0.0', Srp, dprst_srp
        #           ! may need to adjust dprst_srp and volumes
        #           Srp = 0.0
        #         ENDIF
        #       ENDIF
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

        #
        #       IF ( Sri>0.0 ) THEN
        #         tmp = Sri*Imperv_frac*Sro_to_dprst_imperv*Hruarea
        #         IF ( Dprst_area_open_max>0.0 ) THEN
        #           dprst_sri_open = tmp*Dprst_frac_open(Ihru)
        #           dprst_sri = dprst_sri_open/Hruarea
        #           Dprst_vol_open = Dprst_vol_open + DBLE( dprst_sri_open )
        #         ENDIF
        #         IF ( Dprst_area_clos_max>0.0 ) THEN
        #           dprst_sri_clos = tmp*Dprst_frac_clos(Ihru)
        #           dprst_sri = dprst_sri + dprst_sri_clos/Hruarea
        #           Dprst_vol_clos = Dprst_vol_clos + DBLE( dprst_sri_clos )
        #         ENDIF
        #         Sri = Sri - dprst_sri/Imperv_frac
        #         IF ( Sri<0.0 ) THEN
        #           IF ( Sri<-NEARZERO ) PRINT *, 'dprst sri<0.0', Sri, dprst_sri
        #           ! may need to adjust dprst_sri and volumes
        #           Sri = 0.0
        #         ENDIF
        #       ENDIF
        #       Dprst_insroff_hru(Ihru) = dprst_srp + dprst_sri
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

        # cdl -- todo: did not implement
        #       IF ( Dprst_transfer_water_use==ACTIVE ) THEN
        #         IF ( Dprst_area_open_max>0.0 ) THEN
        #           IF ( Dprst_transfer(Ihru)>0.0 ) THEN
        #             IF ( SNGL(Dprst_vol_open*Cfs_conv)<Dprst_transfer(Ihru) ) THEN
        #               PRINT *, 'ERROR, not enough storage for transfer from open surface-depression storage:', &
        #      &                  Ihru, ' Date:', Nowyear, Nowmonth, Nowday
        #               PRINT *, '       storage: ', Dprst_vol_open, '; transfer: ', Dprst_transfer(Ihru)/Cfs_conv
        #               ERROR STOP ERROR_water_use
        #             ENDIF
        #             Dprst_vol_open = Dprst_vol_open - DBLE( Dprst_transfer(Ihru)*Dprst_area_open_max ) / Cfs_conv
        #           ENDIF
        #         ENDIF
        #         IF ( Dprst_area_clos_max>0.0 ) THEN
        #           IF ( Dprst_transfer(Ihru)>0.0 ) THEN
        #             IF ( SNGL(Dprst_area_clos_max*Cfs_conv)<Dprst_transfer(Ihru) ) THEN
        #               PRINT *, 'ERROR, not enough storage for transfer from closed surface-depression storage:', &
        #      &                  Ihru, ' Date:', Nowyear, Nowmonth, Nowday
        #               PRINT *, '       storage: ', Dprst_vol_clos, '; transfer: ', Dprst_transfer(Ihru)/Cfs_conv
        #               ERROR STOP ERROR_water_use
        #             ENDIF
        #             Dprst_vol_clos = Dprst_vol_clos - DBLE( Dprst_transfer(Ihru)*Dprst_area_clos_max ) / Cfs_conv
        #           ENDIF
        #         ENDIF
        #       ENDIF
        #
        # !     Open depression surface area for each HRU:
        #       Dprst_area_open = 0.0
        #       IF ( Dprst_vol_open>0.0D0 ) THEN
        #         open_vol_r = SNGL( Dprst_vol_open/Dprst_vol_open_max )
        #         IF ( open_vol_r<NEARZERO ) THEN
        #           frac_op_ar = 0.0
        #         ELSEIF ( open_vol_r>1.0 ) THEN
        #           frac_op_ar = 1.0
        #         ELSE
        #           frac_op_ar = EXP(Va_open_exp(Ihru)*LOG(open_vol_r))
        #         ENDIF
        #         Dprst_area_open = Dprst_area_open_max*frac_op_ar
        #         IF ( Dprst_area_open>Dprst_area_open_max ) Dprst_area_open = Dprst_area_open_max
        # !        IF ( Dprst_area_open<NEARZERO ) Dprst_area_open = 0.0
        #       ENDIF
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

        #
        # !     Closed depression surface area for each HRU:
        #       IF ( Dprst_area_clos_max>0.0 ) THEN
        #         Dprst_area_clos = 0.0
        #         IF ( Dprst_vol_clos>0.0D0 ) THEN
        #           clos_vol_r = SNGL( Dprst_vol_clos/Dprst_vol_clos_max(Ihru) )
        #           IF ( clos_vol_r<NEARZERO ) THEN
        #             frac_cl_ar = 0.0
        #           ELSEIF ( clos_vol_r>1.0 ) THEN
        #             frac_cl_ar = 1.0
        #           ELSE
        #             frac_cl_ar = EXP(Va_clos_exp(Ihru)*LOG(clos_vol_r))
        #           ENDIF
        #           Dprst_area_clos = Dprst_area_clos_max*frac_cl_ar
        #           IF ( Dprst_area_clos>Dprst_area_clos_max ) Dprst_area_clos = Dprst_area_clos_max
        # !          IF ( Dprst_area_clos<NEARZERO ) Dprst_area_clos = 0.0
        #         ENDIF
        #       ENDIF
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
        #       ! evaporate water from depressions based on snowcov_area
        #       ! dprst_evap_open & dprst_evap_clos = inches-acres on the HRU
        #       unsatisfied_et = Avail_et
        #       dprst_avail_et = (Potet(Ihru)*(1.0-Snowcov_area(Ihru)))*Dprst_et_coef(Ihru)
        #       Dprst_evap_hru = 0.0
        #       IF ( dprst_avail_et>0.0 ) THEN
        #         dprst_evap_open = 0.0
        #         dprst_evap_clos = 0.0
        #         IF ( Dprst_area_open>0.0 ) THEN
        #           dprst_evap_open = MIN(Dprst_area_open*dprst_avail_et, SNGL(Dprst_vol_open))
        #           IF ( dprst_evap_open/Hruarea>unsatisfied_et ) THEN
        #             dprst_evap_open = unsatisfied_et*Hruarea
        #           ENDIF
        #           IF ( dprst_evap_open>SNGL(Dprst_vol_open) ) dprst_evap_open = SNGL( Dprst_vol_open )
        #           unsatisfied_et = unsatisfied_et - dprst_evap_open/Hruarea
        #           Dprst_vol_open = Dprst_vol_open - DBLE( dprst_evap_open )
        #         ENDIF
        #         IF ( Dprst_area_clos>0.0 ) THEN
        #           dprst_evap_clos = MIN(Dprst_area_clos*dprst_avail_et, SNGL(Dprst_vol_clos))
        #           IF ( dprst_evap_clos/Hruarea>unsatisfied_et ) THEN
        #             dprst_evap_clos = unsatisfied_et*Hruarea
        #           ENDIF
        #           IF ( dprst_evap_clos>SNGL(Dprst_vol_clos) ) dprst_evap_clos = SNGL( Dprst_vol_clos )
        #           Dprst_vol_clos = Dprst_vol_clos - DBLE( dprst_evap_clos )
        #         ENDIF
        #         Dprst_evap_hru = (dprst_evap_open + dprst_evap_clos)/Hruarea
        #       ENDIF
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

        #
        #       ! compute seepage
        #       Dprst_seep_hru = 0.0D0
        #       IF ( Dprst_vol_open>0.0D0 ) THEN
        #         seep_open = Dprst_vol_open*DBLE( Dprst_seep_rate_open(Ihru) )
        #         Dprst_vol_open = Dprst_vol_open - seep_open
        #         IF ( Dprst_vol_open<0.0D0 ) THEN
        # !          IF ( Dprst_vol_open<-DNEARZERO ) PRINT *, 'negative dprst_vol_open:', Dprst_vol_open, ' HRU:', Ihru
        #           seep_open = seep_open + Dprst_vol_open
        #           Dprst_vol_open = 0.0D0
        #         ENDIF
        #         Dprst_seep_hru = seep_open/Hruarea_dble
        #       ENDIF
        dprst_seep_hru = 0.0
        if dprst_vol_open > 0.0:
            seep_open = dprst_vol_open * self.dprst_seep_rate_open[ihru]
            dprst_vol_open = dprst_vol_open - seep_open
            if dprst_vol_open < 0.0:
                seep_open = seep_open + dprst_vol_open
                dprst_vol_open = 0.0
            dprst_seep_hru = seep_open / self.hru_area[ihru]

        #
        #       ! compute open surface runoff
        #       Dprst_sroff_hru = 0.0D0
        #       IF ( Dprst_vol_open>0.0D0 ) THEN
        #         Dprst_sroff_hru = MAX( 0.0D0, Dprst_vol_open-Dprst_vol_open_max )
        #         Dprst_sroff_hru = Dprst_sroff_hru + &
        #      &                    MAX( 0.0D0, (Dprst_vol_open-Dprst_sroff_hru-Dprst_vol_thres_open(Ihru))*DBLE(Dprst_flow_coef(Ihru)) )
        #         Dprst_vol_open = Dprst_vol_open - Dprst_sroff_hru
        #         Dprst_sroff_hru = Dprst_sroff_hru/Hruarea_dble
        #         ! sanity checks
        #         IF ( Dprst_vol_open<0.0D0 ) THEN
        # !          IF ( Dprst_vol_open<-DNEARZERO ) PRINT *, 'issue, dprst_vol_open<0.0', Dprst_vol_open
        #           Dprst_vol_open = 0.0D0
        #         ENDIF
        #       ENDIF
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

        #
        #       IF ( Dprst_area_clos_max>0.0 ) THEN
        #         IF ( Dprst_area_clos>NEARZERO ) THEN
        #           seep_clos = Dprst_vol_clos*DBLE( Dprst_seep_rate_clos(Ihru) )
        #           Dprst_vol_clos = Dprst_vol_clos - seep_clos
        #           IF ( Dprst_vol_clos<0.0D0 ) THEN
        # !            IF ( Dprst_vol_clos<-DNEARZERO ) PRINT *, 'issue, dprst_vol_clos<0.0', Dprst_vol_clos
        #             seep_clos = seep_clos + Dprst_vol_clos
        #             Dprst_vol_clos = 0.0D0
        #           ENDIF
        #           Dprst_seep_hru = Dprst_seep_hru + seep_clos/Hruarea_dble
        #         ENDIF
        #       ENDIF
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

        #
        #       Basin_dprst_volop = Basin_dprst_volop + Dprst_vol_open
        #       Basin_dprst_volcl = Basin_dprst_volcl + Dprst_vol_clos
        #       Basin_dprst_evap = Basin_dprst_evap + DBLE( Dprst_evap_hru*Hruarea )
        #       Basin_dprst_seep = Basin_dprst_seep + Dprst_seep_hru*Hruarea_dble
        #       Basin_dprst_sroff = Basin_dprst_sroff + Dprst_sroff_hru*Hruarea_dble
        #       Avail_et = Avail_et - Dprst_evap_hru
        #       IF ( Dprst_vol_open_max>0.0 ) Dprst_vol_open_frac(Ihru) = SNGL( Dprst_vol_open/Dprst_vol_open_max )
        #       IF ( Dprst_vol_clos_max(Ihru)>0.0 ) Dprst_vol_clos_frac(Ihru) = SNGL( Dprst_vol_clos/Dprst_vol_clos_max(Ihru) )
        #       if (Dprst_vol_open_max+Dprst_vol_clos_max(Ihru) > 0.0 ) then  !RGN added to avoid divide by zero 1/20/2022
        #         Dprst_vol_frac(Ihru) = SNGL( (Dprst_vol_open+Dprst_vol_clos)/(Dprst_vol_open_max+Dprst_vol_clos_max(Ihru)) )
        #       end if
        #       Dprst_stor_hru(Ihru) = (Dprst_vol_open+Dprst_vol_clos)/Hruarea_dble
        #
        #       END SUBROUTINE dprst_comp
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
        )
