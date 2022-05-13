from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control
from ..constants import HruType, zero

adaptable = Union[str, np.ndarray, Adapter]

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
        params: PrmsParameters,
        soil_moist: adaptable,
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
        verbose: bool = False,
    ):
        self.name = "PRMSRunoff"
        verbose = True
        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
            subclass_name=self.name,
        )

        # store dependencies
        self._input_variables_dict = {}
        self._input_variables_dict["soil_moist"] = adapter_factory(
            soil_moist, "soil_moist"
        )
        self._input_variables_dict["net_rain"] = adapter_factory(
            net_rain, "net_rain"
        )
        self._input_variables_dict["net_ppt"] = adapter_factory(
            net_ppt, "net_ppt"
        )
        self._input_variables_dict["net_snow"] = adapter_factory(
            net_snow, "net_snow"
        )
        self._input_variables_dict["potet"] = adapter_factory(potet, "potet")
        self._input_variables_dict["snowmelt"] = adapter_factory(
            snowmelt, "snowmelt"
        )
        self._input_variables_dict["snow_evap"] = adapter_factory(
            snow_evap, "snow_evap"
        )
        self._input_variables_dict["pkwater_equiv"] = adapter_factory(
            pkwater_equiv, "pkwater_equiv"
        )
        self._input_variables_dict["pptmix_nopack"] = adapter_factory(
            pptmix_nopack, "pptmix_nopack"
        )
        self._input_variables_dict["snowcov_area"] = adapter_factory(
            snowcov_area, "snowcov_area"
        )
        self._input_variables_dict["hru_intcpevap"] = adapter_factory(
            hru_intcpevap, "hru_intcpevap"
        )
        self._input_variables_dict["intcp_changeover"] = adapter_factory(
            intcp_changeover, "intcp_changeover"
        )

        # call the basin_init hack to calculate basin
        # variables
        self.basin_init()

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
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get input variables

        Returns:
            variables: input variables

        """
        return (
            "soil_moist",
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
        }

    def advance(self):
        """Advance box
        Returns:
            None

        """
        # self.var1_stor_old = self.var1_stor
        self._itime_step += 1

        for key, value in self._input_variables_dict.items():
            value.advance()
            v = getattr(self, key)
            v[:] = value.current

        return

    def calculate(self, time_length, vectorized=False):
        """Calculate terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """
        if vectorized:
            self.calculate_vectorized(time_length)
        else:
            # self.calculate_procedural(time_length)
            self.calculate_prms_style()
        return

    def calculate_procedural(self, time_length):

        # Make calculations
        for i in range(self.nhru):

            self.infil[i] = 0.0
            srp = 0.0

            # If rain/snow event with no antecedent snowpack, compute the runoff from the
            # rain first and then proceed with the snowmelt computations.
            self.infil[i] = self.infil[i] + self.net_rain[i]

            precip = self.net_rain[i]
            self.contrib_fraction[i] = self.calculate_contrib_fraction(
                self.carea_max[i],
                self.smidx_coef[i],
                self.smidx_exp[i],
                self.soil_moist[i],
                precip,
            )

            self.infil[i], srp = self.compute_infil_srp(
                self.contrib_fraction[i], self.net_rain[i], self.infil[i], srp
            )

            runoff = runoff + srp * perv_area + sri * Hruarea_imperv
            srunoff = runoff / Hruarea_dble

        return

    def calculate_vectorized(self, time_length):

        # Make calculations
        precip = self.net_rain
        self.contrib_fraction = self.calculate_contrib_fraction(
            self.carea_max,
            self.smidx_coef,
            self.smidx_exp,
            self.soil_moist,
            precip,
        )

        return

    @staticmethod
    def calculate_contrib_fraction(
        carea_max, smidx_coef, smidx_exp, soil_moist, precip
    ):
        smidx = soil_moist + (0.5 * precip)
        contrib_frac = smidx_coef * 10.0 ** (smidx_exp * smidx)
        contrib_frac = min(carea_max, contrib_frac)
        return contrib_frac

    @staticmethod
    def compute_infil_srp(contrib_frac, precip, infil, perv_runoff):
        adjusted_precip = contrib_frac * precip
        infil = infil - adjusted_precip
        perv_runoff = perv_runoff + adjusted_precip
        return infil, perv_runoff

    def calculate_prms_style(self):

        dprst_chk = 0
        Infil = 0.0
        # DO k = 1, Active_hrus
        for k in range(self.nhru):

            # cdl i = Hru_route_order(k)
            i = k
            ihru = i

            hruarea = self.hru_area[i]
            Hruarea_dble = hruarea
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
            if hruarea_imperv > 0.0:
                imperv_frac = self.hru_percent_imperv[i]
                self.hru_sroffi[i] = 0.0
                self.imperv_evap[i] = 0.0
                # cdl -- not sure what this is:  Hru_impervevap(i) = 0.0

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

            # cdl -- TODO depression storage needs to be implemented for NHM
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
            smidx = self.soil_moist[ihru] + 0.5 * ptc
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
        capacity = self.soil_moist_max[ihru] - self.soil_moist[ihru]
        excess = infil - capacity
        if excess > snowinfil_max:
            srp = srp + excess - snowinfil_max
            infil = snowinfil_max + capacity
        return infil, srp
