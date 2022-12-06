import numpy as np

from pynhm.base.storageUnit import StorageUnit

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import HruType, epsilon32, epsilon64, inch2cm, nan, one, zero

# These are constants used like variables (on self) in PRMS6
# They dont appear on any LHS, so it seems they are constants
# These should probably be in the parameters?
acum_init = [
    0.80,
    0.77,
    0.75,
    0.72,
    0.70,
    0.69,
    0.68,
    0.67,
    0.66,
    0.65,
    0.64,
    0.63,
    0.62,
    0.61,
    0.60,
]

amlt_init = [
    0.72,
    0.65,
    0.60,
    0.58,
    0.56,
    0.54,
    0.52,
    0.50,
    0.48,
    0.46,
    0.44,
    0.43,
    0.42,
    0.41,
    0.40,
]

maxalb = 15
ONETHIRD = 1.0 / 3.0

tcind = 0

dbgind = 434


class PRMSSnow(StorageUnit):
    """PRMS snow pack."""

    def __init__(
        self,
        control: Control,
        orad_hru: adaptable,
        soltab_horad_potsw: adaptable,
        swrad: adaptable,
        hru_intcpevap: adaptable,
        hru_ppt: adaptable,
        potet: adaptable,
        pptmix: adaptable,
        prmx: adaptable,
        tavgc: adaptable,
        tmaxc: adaptable,
        tminc: adaptable,
        net_ppt: adaptable,
        net_rain: adaptable,
        net_snow: adaptable,
        transp_on: adaptable,
        budget_type: str = None,
        verbose: bool = False,
        load_n_time_batches: int = 1,
    ) -> "PRMSSnow":

        super().__init__(
            control=control,
            verbose=verbose,
            load_n_time_batches=load_n_time_batches,
        )
        self.name = "PRMSSnow"
        self._set_inputs(locals())
        self._set_budget(budget_type)

        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get snow pack parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nhru",
            "nmonth",
            "ndeplval",
            "cov_type",
            "covden_win",
            "covden_sum",
            "hru_type",
            "albset_rna",
            "albset_rnm",
            "albset_sna",
            "albset_snm",
            "den_init",
            "den_max",
            "settle_const",
            "emis_noppt",
            "freeh2o_cap",
            "hru_deplcrv",
            "melt_force",
            "melt_look",
            "potet_sublim",
            "rad_trncf",
            "snarea_curve",
            "snarea_thresh",
            "snowpack_init",
            "tmax_allsnow",
            "cecn_coef",
            "tstorm_mo",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get snow pack input variables

        Returns:
            variables: input variables

        """
        return (
            "hru_ppt",
            "hru_intcpevap",
            "net_ppt",
            "net_rain",
            "net_snow",
            "orad_hru",
            "potet",
            "pptmix",
            "prmx",
            "soltab_horad_potsw",
            "swrad",
            "tavgc",
            "tmaxc",
            "tminc",
            "transp_on",
        )

    @staticmethod
    def get_init_values() -> dict:
        """Get snow pack initial values

        Returns:
            dict: initial values for named variables
        """

        return {
            "ai": zero,
            "albedo": zero,
            "frac_swe": zero,
            "freeh2o": zero,
            "freeh2o_change": nan,
            "freeh2o_prev": nan,
            "iasw": False,
            "int_alb": one,
            "iso": one,
            "lso": zero,
            "lst": False,
            "mso": one,
            "pk_def": zero,
            "pk_den": zero,
            "pk_depth": zero,
            "pk_ice": zero,
            "pk_ice_change": nan,
            "pk_ice_prev": nan,
            "pk_precip": zero,
            "pk_temp": zero,
            "pksv": zero,
            "pkwater_ante": nan,
            "pkwater_equiv": nan,
            "pkwater_equiv_change": nan,
            "pptmix_nopack": False,
            "pss": nan,
            "pst": nan,
            "salb": zero,
            "scrv": zero,
            "slst": zero,
            "snow_evap": zero,
            "snowcov_area": zero,
            "snowcov_areasv": zero,
            "snowmelt": zero,
            "snsv": zero,
            "tcal": zero,
            "through_rain": nan,
        }

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": ["net_rain", "net_snow"],
            "outputs": [
                "snow_evap",
                "snowmelt",
                "through_rain",
            ],
            "storage_changes": [
                "freeh2o_change",
                "pk_ice_change"
                # eventaully pakwater_equiv_change should be removed
                # entirely from the code and the metadata.
                # But there are significant gaps between the above and
                # pkwater_equiv
                # "pkwater_equiv_change",
            ],
        }

    @staticmethod
    def get_restart_variables() -> tuple:
        """Get snow pack restart variables

        Returns:
            variables: restart variables
        """

        return (
            "albedo",
            "freeh2o",
            "iasw",
            "int_alb",
            "iso",
            "lso",
            "lst",
            "mso",
            "pk_def",
            "pk_depth",
            "pk_den",
            "pk_ice",
            "pk_temp",
            "pksv",
            "pss",
            "pst",
            "salb",
            "scrv",
            "slst",
            "snowcov_area",
            "snowcov_areasv",
            "snsv",
        )

    def _set_initial_conditions(self):
        """Initialize PRMSSnow snowpack variables."""

        # Derived parameters
        self.tmax_allsnow_c = (self.tmax_allsnow - 32.0) / 1.8
        del self.tmax_allsnow

        # Deninv and denmaxinv not in variables nor in metadata but we can set
        # them on self for convenience
        if (self.den_init.shape == ()) or (
            self.den_init.shape[0] != self.nhru
        ):
            den_init = np.ones(self.nhru, dtype=np.float64) * self.den_init
        else:
            den_init = self.den_init

        self.deninv = one / den_init
        self.denmaxinv = one / self.den_max.copy()

        self.pkwater_equiv[:] = self.snowpack_init.copy()

        sd = int(self.ndeplval / 11)
        self.snarea_curve_2d = np.reshape(self.snarea_curve, (sd, 11))

        if self.control.config["init_vars_from_file"] in [0, 2, 3]:

            # The super().__init__ already set_initial_conditions using its
            # set_initial_conditions
            # Below Im just following PRMS6, will reconcile later with the
            # super (may be redundant).
            vars_init = [
                "albedo",
                "iasw",
                "int_alb",
                "iso",
                "lso",
                "lst",
                "mso",
                "pk_def",
                "pk_temp",
                "pksv",
                "salb",
                "scrv",
                "slst",
                "snowcov_areasv",
                "snsv",
            ]
            for vv in vars_init:
                self._initialize_var(vv)

            pkweq_gt_zero = self.pkwater_equiv > zero
            wh_pkweq_gt_zero = np.where(pkweq_gt_zero)
            self.pk_depth[wh_pkweq_gt_zero] = (
                self.pkwater_equiv[wh_pkweq_gt_zero]
                * self.deninv[wh_pkweq_gt_zero]
            )
            self.pk_den[wh_pkweq_gt_zero] = (
                self.pkwater_equiv[wh_pkweq_gt_zero]
                / self.pk_depth[wh_pkweq_gt_zero]
            )
            self.pk_ice[wh_pkweq_gt_zero] = self.pkwater_equiv[
                wh_pkweq_gt_zero
            ]
            self.freeh2o[wh_pkweq_gt_zero] = (
                self.pk_ice[wh_pkweq_gt_zero]
                * self.freeh2o_cap[wh_pkweq_gt_zero]
            )
            self.ai[wh_pkweq_gt_zero] = self.pkwater_equiv[
                wh_pkweq_gt_zero
            ]  # inches

            ai_gt_snarea_thresh = self.ai > self.snarea_thresh
            wh_pkweq_gt_zero_and_ai_gt_snth = np.where(
                pkweq_gt_zero & ai_gt_snarea_thresh
            )
            self.ai[wh_pkweq_gt_zero_and_ai_gt_snth] = self.snarea_thresh[
                wh_pkweq_gt_zero_and_ai_gt_snth
            ]

            for ww in range(len(wh_pkweq_gt_zero[0])):
                self.snowcov_area[wh_pkweq_gt_zero[ww]] = self.sca_deplcrv(
                    self.snarea_curve_2d[
                        1:11, self.hru_deplcrv[wh_pkweq_gt_zero[ww] - 1]
                    ],
                    self.frac_swe[wh_pkweq_gt_zero[ww]],
                )

            self.pss[:] = self.pkwater_equiv.copy()
            self.pst[:] = self.pkwater_equiv.copy()

        else:

            raise RuntimeError("Snow restart capability not implemented")
            # JLM: a list of restart variables dosent shed light on what states actually have memory.
            # JLM: could there be a diagnostic part of the advance?

        return

    def _advance_variables(self) -> None:
        self.pkwater_ante[:] = self.pkwater_equiv
        self.freeh2o_prev[:] = self.freeh2o
        self.pk_ice_prev[:] = self.pk_ice
        return

    def _calculate(self, simulation_time):
        """Calculate snow pack terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None
        """

        self._simulation_time = simulation_time

        self.canopy_covden = self.covden_win
        wh_transp_on = np.where(self.transp_on)
        self.canopy_covden[wh_transp_on] = self.covden_sum[wh_transp_on]

        cals = zero  # JLM this is unnecessary.

        # newsnow
        self.newsnow = np.full([self.nhru], False, dtype=bool)
        net_snow_gt_zero = self.net_snow > zero
        wh_net_snow_gt_zero = np.where(net_snow_gt_zero)
        self.newsnow[wh_net_snow_gt_zero] = True

        # default assumption
        self.pptmix_nopack[:] = False

        for jj in range(self.nhru):

            if self.hru_type[jj] == HruType.LAKE:
                continue

            # <
            trd = self.orad_hru[jj] / self.soltab_horad_potsw[jj]

            # If it's the first julian day of the water year, several
            # variables need to be reset:
            # - reset the previous snow water eqivalent plus new snow to 0
            # - reset flags to indicate it is not melt season or potetential melt season
            # - reset the counter for the number of days a snowpack is at 0 deg Celsius
            # TODO: rsr, do we want to reset all HRUs, what about Southern Hemisphere
            if self.control.current_dowy == 1:
                self.pss[jj] = zero  # [inches]
                self.iso[jj] = 1  # [flag]
                self.mso[jj] = 1  # [flag]
                self.lso[jj] = 0  # [counter]

            # <
            # HRU SET-UP - SET DEFAULT VALUES AND/OR BASE CONDITIONS FOR THIS TIME PERIOD
            # **************************************************************

            # By default, the precipitation added to snowpack, snowmelt,
            # and snow evaporation are 0.
            # JLM: this could happen outside the loop
            self.pk_precip[jj] = zero  # [inches]
            self.snowmelt[jj] = zero  # [inches]
            self.snow_evap[jj] = zero  # [inches]
            self.frac_swe[jj] = zero
            self.ai[jj] = zero
            self.tcal[jj] = zero

            # If the day of the water year is beyond the forced melt day indicated by
            # the parameter, then set the flag indicating melt season
            # TODO: rsr, need to rethink this at some point
            if self.control.current_doy == self.melt_force[jj]:
                self.iso[jj] = 2  # [flag]  # JLM: why not use enumerator here?

            # <
            # If the day of the water year is beyond the first day to look for melt
            # season indicated by the parameter, then set the flag indicating to watch
            # for melt season.
            # TODO: rsr, need to rethink this at some point
            if self.control.current_doy == self.melt_look[jj]:
                self.mso[jj] = 2  # [flag]  # could emume this BEFORE/AFTER

            # <
            # if jj == dbgind:
            #     print(f"self.pkwater_equiv 0 : {self.pkwater_equiv[dbgind]}")
            #     print(f"self.pk_ice 0 : {self.pk_ice[dbgind]}")
            #     print(f"self.tcal 0 : {self.tcal[dbgind]}")

            if self.pkwater_equiv[jj] < epsilon64:
                # No existing snowpack
                if not self.newsnow[jj]:
                    # Skip the HRU if there is no snowpack and no new snow
                    # Reset to be sure it is zero if snowpack melted on last timestep.
                    self.snowcov_area[jj] = zero
                    continue
                else:
                    # We ahave new snow; the initial snow-covered area is complete (1)
                    # JLM: why set this here? just for the case of no existing snow?
                    # This might be removable.
                    self.snowcov_area[jj] = one  # [fraction of area]

            # <<
            # HRU STEP 1 - DEAL WITH PRECIPITATION AND ITS EFFECT ON THE WATER
            #              CONTENT AND HEAT CONTENT OF SNOW PACK
            # ***********************************************************************
            self.ppt_to_pack(jj)

            # if jj == dbgind:
            #     print(f"self.pkwater_equiv 1 : {self.pkwater_equiv[dbgind]}")
            #     print(f"self.net_rain: {self.net_rain[dbgind]}", flush=True)
            #     print(f"self.tcal 1 : {self.tcal[dbgind]}")

            # <
            if self.pkwater_equiv[jj] > zero:
                # If there is still snow after precipitation
                # HRU STEP 2 - CALCULATE THE NEW SNOW COVERED AREA from depletion curve
                # **********************************************************
                self.snowcov(jj)
                # if jj == dbgind:
                #     print(
                #         f"self.pkwater_equiv 2 : {self.pkwater_equiv[dbgind]}"
                #     )
                #     print(f"self.snowcov_area: {self.snowcov_area[dbgind]}")
                #     print(f"self.tcal 2 : {self.tcal[dbgind]}")

                # if self.control._itime_step == 29 and jj == dbgind:
                #    pdb.set_trace()

                # HRU STEP 3 - COMPUTE THE NEW ALBEDO
                # **********************************************************
                self.snalbedo(jj)
                # if jj == dbgind:
                #     print(
                #         f"self.pkwater_equiv 3 : {self.pkwater_equiv[dbgind]}"
                #     )
                #     print(f"self.tcal 3 : {self.tcal[dbgind]}")

                # HRU STEP 4 - DETERMINE RADIATION FLUXES AND SNOWPACK
                #              STATES NECESSARY FOR ENERGY BALANCE
                # **********************************************************
                self.step_4(jj, trd)
                # if jj == dbgind:
                #     print(
                #         f"self.pkwater_equiv 4 : {self.pkwater_equiv[dbgind]}"
                #     )
                #     print(f"self.tcal 4 : {self.tcal[dbgind]}")

                #  HRU STEP 5 - CALCULATE SNOWPACK LOSS TO EVAPORATION
                # ********************************************************
                if self.pkwater_equiv[jj] > zero:
                    # Snow can evaporate when transpiration is not occuring or when
                    # transpiration is occuring with cover types of bare soil or grass.
                    # if (transp_on[jj] == 0 .or. (transp_on[jj] == 1 .and. cov_type[jj] < 2):
                    if (not self.transp_on[jj]) or (
                        self.transp_on[jj] and self.cov_type[jj] < 2
                    ):
                        (
                            self.freeh2o[jj],
                            self.pk_def[jj],
                            self.pk_ice[jj],
                            self.pk_temp[jj],
                            self.pkwater_equiv[jj],
                            self.snow_evap[jj],
                        ) = self._calc_snowevap(
                            self.freeh2o[jj],
                            self.hru_intcpevap[jj],
                            self.pk_def[jj],
                            self.pk_ice[jj],
                            self.pk_temp[jj],
                            self.pkwater_equiv[jj],
                            self.potet[jj],
                            self.potet_sublim[jj],
                            self.snow_evap[jj],
                            self.snowcov_area[jj],
                            self.verbose,
                        )

                # <<
                elif self.pkwater_equiv[jj] < zero:
                    if self.verbose:
                        if self.pkwater_equiv[jj] < (-1 * epsilon64):
                            print(
                                f"snowpack issue 3, negative pkwater_equiv, "
                                f"HRU: {jj}, value: {self.pkwater_equiv[jj]}"
                            )

                    # <<
                    # just to be sure negative values are ignored
                    self.pkwater_equiv[jj] = zero

                # if jj == dbgind:
                #     print(
                #         f"self.pkwater_equiv 5 : {self.pkwater_equiv[dbgind]}"
                #     )
                #     print(f"self.tcal 5 : {self.tcal[dbgind]}")
                # <
                # HRU CLEAN-UP - ADJUST FINAL HRU SNOWPACK STATES AND
                #                INCREMENT THE BASIN TOTALS
                # *********************************************************

                # Final state of the snowpack depends on whether it still exists after
                # all the processing above.
                # 2 options below (if-then, else)

                if self.pkwater_equiv[jj] > zero:
                    # (1) Snow pack still exists
                    # Snowpack still exists

                    if self.pk_den[jj] > zero:
                        self.pk_depth[jj] = (
                            self.pkwater_equiv[jj] / self.pk_den[jj]
                        )
                    else:
                        self.pk_den[jj] = self.den_max
                        self.pk_depth[jj] = (
                            self.pkwater_equiv[jj] * self.denmaxinv
                        )

                    # <
                    self.pss[jj] = self.pkwater_equiv[jj]

                    # If it is during the melt period and snowfall was insufficient to
                    # reset albedo, then reduce the cumulative new snow by the amount
                    # melted during the period (but don't let it be negative).
                    # if (self.lst[jj] > 0:
                    if self.lst[jj]:
                        self.snsv[jj] = self.snsv[jj] - self.snowmelt[jj]

                        if self.snsv[jj] < zero:
                            self.snsv[jj] = zero

                # <<<<
                # LAST check to clear out all arrays if packwater is gone
                if self.pkwater_equiv[jj] <= zero:
                    if self.verbose:
                        if self.pkwater_equiv[jj] < -epsilon64:
                            print(
                                f"Snowpack problem, pkwater_equiv negative, HRU: {jj}, value: {self.pkwater_equiv[jj]}"
                            )

                    # <<
                    self.pkwater_equiv[
                        jj
                    ] = zero  # just to be sure negative values are ignored

                    # Snowpack has been completely depleted, reset all states to no-snowpack values
                    self.pk_depth[jj] = zero
                    self.pss[jj] = zero
                    self.snsv[jj] = zero
                    self.lst[jj] = False
                    self.pst[jj] = zero
                    self.iasw[jj] = False
                    self.albedo[jj] = zero
                    self.pk_den[jj] = zero
                    self.snowcov_area[jj] = zero
                    self.pk_def[jj] = zero
                    self.pk_temp[jj] = zero
                    self.pk_ice[jj] = zero
                    self.freeh2o[jj] = zero
                    self.snowcov_areasv[jj] = zero  # rsr, not in original code
                    self.ai[jj] = zero
                    self.frac_swe[jj] = zero

            # <<

        self.pkwater_equiv_change[:] = self.pkwater_equiv - self.pkwater_ante
        self.freeh2o_change[:] = self.freeh2o - self.freeh2o_prev
        self.pk_ice_change[:] = self.pk_ice - self.pk_ice_prev

        wh_through = (
            ((self.pk_ice_prev + self.freeh2o_prev) <= epsilon64)
            & ~self.newsnow
        ) | (self.pptmix_nopack == 1)
        self.through_rain[:] = np.where(wh_through, self.net_rain, zero)

        return

    @staticmethod
    def sca_deplcrv(snarea_curve: np.ndarray, frac_swe: float) -> float:
        """Interpolate along snow covered area depletion curve"""
        if frac_swe > one:
            res = snarea_curve[-1]
        else:
            # Get the indices (as integers) of the depletion curve that bracket the
            # given frac_swe (next highest and next lowest).
            idx = int(10.0 * (frac_swe + 0.2))  # [index]
            jdx = idx - 1  # [index]
            if idx > (11):
                idx = 11
            # Calculate the fraction of the distance (from the next lowest) the given
            # frac_swe is between the next highest and lowest curve values.
            dify = (frac_swe * 10.0) - float(jdx - 1)  # [fraction]
            # Calculate the difference in snow covered area represented by next
            # highest and lowest curve values.
            difx = snarea_curve[idx - 1] - snarea_curve[jdx - 1]
            # Linearly interpolate a snow covered area between those represented by
            # the next highest and lowest curve values.
            res = snarea_curve[jdx - 1] + dify * difx

        # <
        return res

    def ppt_to_pack(self, jj):
        """Add rain and/or snow to snowpack."""

        # WARNING: pan - wouldn't this be pkwater_equiv > DNEARZERO?
        ppt_through = ~(
            (self.pkwater_equiv[jj] > zero and self.net_ppt[jj] > zero)
            or self.net_snow[jj] > zero
        )
        if ppt_through:
            return

        caln = zero
        calpr = zero
        calps = zero
        pndz = zero

        month_ind = self.control.current_month - 1

        tsnow = self.tavgc[jj]  # [degrees C]

        if self.pptmix[jj] == 1:
            # (1) If precipitation is mixed...
            # If there is any rain, the rain temperature is halfway between the maximum
            # temperature and the allsnow temperature.
            train = (
                self.tmaxc[jj] + self.tmax_allsnow_c[month_ind, jj]
            ) * 0.5  # [degrees C]

            # Temperatures will differ, depending on the presence of existing snowpack.
            # should this be epsilon32?
            if self.pkwater_equiv[jj] > zero:
                # If there is a snowpack, snow temperature is halfway between the minimum
                # daily temperature and maximum temperature for which all precipitation is snow.
                tsnow = (
                    self.tminc[jj] + self.tmax_allsnow_c[month_ind, jj]
                ) * 0.5  # [degrees C]

            elif self.pkwater_equiv[jj] < zero:
                # If no existing snowpack, snow temperature is the average temperature for the day.
                self.pkwater_equiv[
                    jj
                ] = zero  # To be sure negative snowpack is ignored

        # <<
        else:
            # (2) If precipitation is all snow or all rain...
            # If there is any rain, the rain temperature is the average temperature.
            train = self.tavgc[jj]  # [degrees C]
            if train < epsilon32:
                # If average temperature is close to freezing, the rain temperature is
                # halfway between the maximum daily temperature and maximum temperature
                # for which all precipitation is snow.
                train = (
                    self.tmaxc[jj] + self.tmax_allsnow_c[month_ind, jj]
                ) * 0.5  # [degrees C]

        # <<
        if train < zero:
            train = zero  # [degrees C] # train can't be < 0
        if tsnow > zero:
            tsnow = zero  # [degrees C] # tsnow can't be > 0

        # Leavesley comments...
        # If snowpack already exists, add rain first, then add snow.  If no
        # antecedent snowpack, rain is already taken care of, so start snowpack with
        # snow.  This subroutine assumes that in a mixed event, the rain will be
        # first and turn to snow as the temperature drops.

        # Rain can only add to the snowpack if a previous snowpack exists, so rain
        # or a mixed event is processed differently when a snowpack exists.
        # 2 options below (if-then, elseif)

        if self.pkwater_equiv[jj] > zero:
            # (1) There is net rain on an existing snowpack...

            if self.net_rain[jj] > zero:
                # Add rain water to pack (rain on snow) and increment the precipitation
                # on the snowpack by the rain water.
                self.pkwater_equiv[jj] = (
                    self.pkwater_equiv[jj] + self.net_rain[jj]
                )  # [inches]
                self.pk_precip[jj] = (
                    self.pk_precip[jj] + self.net_rain[jj]
                )  # [inches]

                # Incoming rain water carries heat that must be added to the snowpack.
                # This heat could both warm the snowpack and melt snow. Handling of this
                # heat depends on the current thermal condition of the snowpack.
                # 2 options below (if-then, else)

                # (1.1) If the snowpack is colder than freezing it has a heat deficit
                # (requires heat to be brought to isothermal at 0 degC)...
                if self.pk_def[jj] > zero:
                    # Calculate the number of calories given up per inch of rain when
                    # cooling it from the current rain temperature to 0 deg C and then
                    # freezing it (liquid to solid state latent heat).
                    # This calculation assumes a volume of an inch of rain over a square cm of area
                    # 80 cal come from freezing 1 cm3 at 0 C (latent heat of fusion is 80 cal/cm^3),
                    # 1 cal from cooling 1cm3 for every degree C (specific heat of water is 1 cal/(cm^3 degC)),
                    # convert from 1 cm depth over 1 square cm to 1 inch depth over 1 square cm (INCH2CM = 2.54 cm/in)
                    caln = (80.000 + train) * inch2cm  # [cal / (in cm^2)]

                    # Calculate the amount of rain in inches (at the current rain temperature)
                    # needed to bring the snowpack to isothermal at 0.
                    pndz = self.pk_def[jj] / caln  # [inches]

                    # The effect of rain on the snowpack depends on if there is not enough,
                    # enough, or more than enough heat in the rain to bring the snowpack
                    # to isothermal at 0 degC or not 3 options below (if-then, elseif, else).

                    if abs(self.net_rain[jj] - pndz) < epsilon32:
                        # (1.1.1) Exactly enough rain to bring pack to isothermal...
                        # Heat deficit and temperature of the snowpack go to 0.
                        self.pk_def[jj] = zero  # [cal/cm^2]
                        self.pk_temp[jj] = zero  # [degrees C]

                        # In the process of giving up its heat, all the net rain freezes and
                        # becomes pack ice.
                        self.pk_ice[jj] = (
                            self.pk_ice[jj] + self.net_rain[jj]
                        )  # [inches]

                    elif self.net_rain[jj] < pndz:
                        # (1.1.2) Rain not sufficient to bring pack to isothermal...
                        # The snowpack heat deficit decreases by the heat provided by rain
                        # and a new snowpack temperature is calculated.
                        # 1.27 is the specific heat of ice (0.5 cal/(cm^3 degC))
                        # times the conversion of cm to inches (2.54 cm/in)
                        self.pk_def[jj] = self.pk_def[jj] - (
                            caln * self.net_rain[jj]
                        )  # [cal/(in cm^3)]
                        self.pk_temp[jj] = (
                            -1
                            * self.pk_def[jj]
                            / (self.pkwater_equiv[jj] * 1.27)
                        )

                        # All the net rain freezes and becomes pack ice
                        self.pk_ice[jj] = self.pk_ice[jj] + self.net_rain[jj]

                    else:
                        # (1.1.3) Rain in excess of amount required to bring pack to isothermal...
                        # Heat deficit and temperature of the snowpack go to 0.
                        self.pk_def[jj] = zero
                        self.pk_temp[jj] = zero
                        # The portion of net rain that brings the snowpack to isothermal freezes.
                        self.pk_ice[jj] = self.pk_ice[jj] + pndz

                        # The rest of the net rain becomes free water in the snowpack.
                        # Note that there cannot be previous freeh2o because the snowpack
                        # had a heat deficit (all water was ice) before this condition was reached.
                        self.freeh2o[jj] = self.net_rain[jj] - pndz

                        # Calculate the excess heat per area added by the portion of rain
                        # that does not bring the snowpack to isothermal (using specific heat of water)
                        calpr = (
                            train * (self.net_rain[jj] - pndz) * inch2cm
                        )  # [cal/cm^2]

                        # Add the new heat to the snow pack (the heat in this excess rain
                        # will melt some of the pack ice when the water cools to 0 degC).

                        # if (chru == 5):
                        #     print *, ' ', pkwater_equiv[jj], month, train, calpr, pndz, net_rain[jj], net_snow[jj]

                        # endif

                        (
                            self.freeh2o[jj],
                            self.iasw[jj],
                            self.pk_def[jj],
                            self.pk_den[jj],
                            self.pk_ice[jj],
                            self.pk_depth[jj],
                            self.pk_temp[jj],
                            self.pss[jj],
                            self.pst[jj],
                            self.snowmelt[jj],
                            self.pkwater_equiv[jj],
                        ) = self._calc_calin(
                            calpr,
                            self.den_max,
                            self.denmaxinv,
                            self.freeh2o[jj],
                            self.freeh2o_cap[jj],
                            self.iasw[jj],
                            self.pk_def[jj],
                            self.pk_den[jj],
                            self.pk_depth[jj],
                            self.pk_ice[jj],
                            self.pk_temp[jj],
                            self.pkwater_equiv[jj],
                            self.pss[jj],
                            self.pst[jj],
                            self.snowcov_area[jj],
                            self.snowmelt[jj],
                        )
                        # if (chru == 5):
                        #     print *, '*', pkwater_equiv[jj]
                        # endif

                # <<
                else:
                    # (1.2) Rain on snowpack that is isothermal at 0 degC (no heat deficit)...
                    # All net rain is added to free water in the snowpack.
                    self.freeh2o[jj] = self.freeh2o[jj] + self.net_rain[jj]

                    # Calculate the heat per area added by the rain (using specific heat of water).
                    calpr = train * self.net_rain[jj] * inch2cm  # [cal/cm^2]

                    # Add the new heat to the snow pack (the heat in rain will melt some
                    # of the pack ice when the water cools to 0 degC).
                    (
                        self.freeh2o[jj],
                        self.iasw[jj],
                        self.pk_def[jj],
                        self.pk_den[jj],
                        self.pk_ice[jj],
                        self.pk_depth[jj],
                        self.pk_temp[jj],
                        self.pss[jj],
                        self.pst[jj],
                        self.snowmelt[jj],
                        self.pkwater_equiv[jj],
                    ) = self._calc_calin(
                        calpr,
                        self.den_max,
                        self.denmaxinv,
                        self.freeh2o[jj],
                        self.freeh2o_cap[jj],
                        self.iasw[jj],
                        self.pk_def[jj],
                        self.pk_den[jj],
                        self.pk_depth[jj],
                        self.pk_ice[jj],
                        self.pk_temp[jj],
                        self.pkwater_equiv[jj],
                        self.pss[jj],
                        self.pst[jj],
                        self.snowcov_area[jj],
                        self.snowmelt[jj],
                    )

        # <<<
        elif self.net_rain[jj] > zero:
            # (2) If there is net rain but no snowpack, set flag for a mix on no snowpack.

            # Be careful with the code here.
            # If this subroutine is called when there is an all-rain day on no
            # existing snowpack (currently, it will not), then the flag here will be
            # set inappropriately.
            self.pptmix_nopack[jj] = True  # [flag]

        # <
        # At this point, the subroutine has handled all conditions where there is
        # net rain, so if there is net snow (doesn't matter if there is a pack or not)...
        if self.net_snow[jj] > zero:
            # Add the new snow to the pack water equivalent, precip, and ice
            self.pkwater_equiv[jj] = self.pkwater_equiv[jj] + self.net_snow[jj]
            self.pk_precip[jj] = self.pk_precip[jj] + self.net_snow[jj]
            self.pk_ice[jj] = self.pk_ice[jj] + self.net_snow[jj]

            # The temperature of the new snow will determine its effect on snowpack heat deficit
            # 2 options below (if-then, else)
            if tsnow >= zero:
                # (1) if the new snow is at least 0 degC...

                # Incoming snow does not change the overall heat content of the snowpack.
                # However, the temperature will change, because the total heat content
                # of the snowpack will be "spread out" among more snow.  Calculate the
                # snow pack temperature from the heat deficit, specific heat of snow,
                # and the new total snowpack water content.
                self.pk_temp[jj] = (
                    -1 * self.pk_def[jj] / (self.pkwater_equiv[jj] * 1.27)
                )  # [degrees C]
                # JLM: i dont see where pk_def is set if there was not existing snowpk

            # <
            else:
                # (2) If the new snow is colder than 0 degC...
                # Calculate the amount of heat the new snow will absorb if warming it
                # to 0C (negative number). This is the negative of the heat deficit of
                # the new snow.
                calps = tsnow * self.net_snow[jj] * 1.27  # [cal/cm^2]

                # The heat to warm the new snow can come from different sources
                # depending on the state of the snowpack.
                # 2 options below (if-then, else)
                if self.freeh2o[jj] > zero:
                    # (2.1) If there is free water in the pack (at least some of it is going to freeze)...
                    (
                        self.freeh2o[jj],
                        self.pk_def[jj],
                        self.pk_ice[jj],
                        self.pk_temp[jj],
                        self.pkwater_equiv[jj],
                    ) = self._calc_caloss(
                        calps,
                        self.freeh2o[jj],
                        self.pk_def[jj],
                        self.pk_ice[jj],
                        self.pk_temp[jj],
                        self.pkwater_equiv[jj],
                    )

                else:
                    # (2.2) If there is no free water (snow pack has a heat deficit greater
                    #       than or equal to 0)...
                    # Heat deficit increases because snow is colder than pack (minus a
                    # negative number = plus) and calculate the new pack temperature.
                    self.pk_def[jj] = self.pk_def[jj] - calps  # [cal/cm^2]
                    self.pk_temp[jj] = (
                        -1 * self.pk_def[jj] / (self.pkwater_equiv[jj] * 1.27)
                    )  # [degrees C]

        # <<<
        return

    @staticmethod
    def _calc_calin(
        cal,
        den_max,
        denmaxinv,
        freeh2o,
        freeh2o_cap,
        iasw,
        pk_def,
        pk_den,
        pk_depth,
        pk_ice,
        pk_temp,
        pkwater_equiv,
        pss,
        pst,
        snowcov_area,
        snowmelt,
    ):
        """Compute changes in snowpack when a net gain in heat energy has occurred."""

        # Local Variables: all doubles
        # apk_ice: Pack-ice per area [inches]
        # apmlt: Actual potential snowmelt [inches]
        # dif: Difference between incoming calories and calories needed to
        #    bring the pack to isothermal at 0 [cal/cm^2]
        # dif_water: Amount of free water in excess of the capacity to hold
        #    free water.
        # pmlt: Potential amount of snowmelt from excess heat in rain [inches]
        # pwcap: Capacity of the snowpack to hold free water [inches]

        # Calculate the difference between the incoming calories and the
        # calories needed to bring the pack to isothermal at 0 (heat deficit).
        dif = cal - pk_def  # [cal/cm^2]

        # The way incoming heat is handled depends on whether there is not
        # enough, just enough, or more than enough heat to overcome the heat
        # deficit of the snowpack.
        # 3 choices below (if-then, elseif, else)

        if dif < zero:
            # (1) Not enough heat to overcome heat deficit...
            # Reduce the heat deficit by the amount of incoming calories and
            # adjust to
            # the new temperature based on new heat deficit.
            pk_def = pk_def - cal  # [cal/cm^2]
            pk_temp = -1 * pk_def / (pkwater_equiv * 1.27)  # [degrees C]

        elif abs(dif) < epsilon32:
            # JLM: I moved this from an else at the bottom to this
            #      point in the conditional chain.
            # JLM: The test had been equality with zero, changed to
            #      less than epsilon32.
            # (2) Just enough heat to overcome heat deficit
            # Set temperature and heat deficit to zero. the pack is "ripe"
            pk_temp = zero  # [degrees C]
            pk_def = zero  # [cal/cm^2]

        elif dif > zero:
            # (3) More than enough heat to overcome heat deficit (melt ice)...
            # Calculate the potential amount of snowmelt from excess heat in
            # rain it takes 203.2 calories / (in cm^2) to melt snow (latent
            # heat of fusion).
            # Convert from 1 cm depth over 1 square cm to
            # 1 inch depth over 1 square cm 80.0*(INCH2CM = 2.54 cm/in) = 203.2
            pmlt = dif / 203.2  # [inches]

            # Actual snowmelt can only come from snow covered area, so to
            # calculate the actual potential snowmelt, the potential snowmelt
            # from snowcovered area must be re-normalized to HRU area (rather
            # than snowcover area). In effect, the potential snowmelt per area
            # is reduced by the fraction of the watershed that is actually
            # covered by snow.
            apmlt = pmlt * snowcov_area  # [inches]

            # Set the heat deficit and temperature of the remaining snowpack
            # to 0.
            pk_def = zero  # [cal/cm^2]
            pk_temp = zero  # [degrees C]

            # The only pack ice that is melted is in the snow covered area, so
            # the pack ice needs to be re-normalized to the snowcovered area
            # (rather than HRU area). In effect, the pack ice per area is
            # increased by the fraction of the watershed that is actually
            # covered by snow.
            if snowcov_area > zero:
                apk_ice = pk_ice / snowcov_area  # [inches]
            else:
                # print *, 'snowcov_area really small, melt all ice',
                # snowcov_area, ' pmlt:', pmlt, ' dif:', dif, ' pk_ice:', pk_ice
                apk_ice = zero

            # <
            # If snow is melting, the heat is handled based on whether all or
            # only part of the pack ice melts.
            # 2 options below (if-then, else)

            if pmlt > apk_ice:
                # (3.1) Heat applied to snow covered area is sufficient to melt all the
                #       ice in that snow pack.
                # All pack water equivalent becomes meltwater.
                snowmelt = snowmelt + pkwater_equiv  # [inches]
                pkwater_equiv = zero  # [inches]
                # iasw = 0  # [flag]
                iasw = False  # [flag]

                # Set all snowpack states to 0.
                # snowcov_area = zero  # [fraction of area] # shouldn't be
                # changed with melt
                pk_def = zero  # [cal / cm^2]
                pk_temp = zero  # [degreees C]
                pk_ice = zero  # [inches]
                freeh2o = zero  # [inches]
                pk_depth = zero  # [inches]
                pss = zero  # [inches]
                pst = zero  # [inches]
                pk_den = zero  # [fraction of depth]

            else:
                # (3.2) Heat only melts part of the ice in the snow pack.
                # Remove actual melt from frozen water and add melt to free
                # water.
                pk_ice = pk_ice - apmlt  # [inches]
                freeh2o = freeh2o + apmlt  # [inches]

                # Calculate the capacity of the snowpack to hold free water
                # according to
                # its current level of frozen water.
                pwcap = freeh2o_cap * pk_ice  # [inches]

                # Calculate the amount of free water in excess of the capacity
                # to hold free water.
                dif_water = freeh2o - pwcap  # [inches]

                # If there is more free water than the snowpack can hold, then
                # there is going to be melt...

                if dif_water > zero:

                    if dif_water > pkwater_equiv:
                        dif_water = pkwater_equiv

                    # <
                    # total packwater decreases by the excess and a new depth
                    # is calculated based on density.
                    pkwater_equiv = pkwater_equiv - dif_water  # [inches]

                    # Free water is at the current capacity.
                    freeh2o = pwcap  # [inches]
                    if pk_den > zero:
                        pk_depth = pkwater_equiv / pk_den  # [inches]
                        # RAPCOMMENT - Added the conditional statement to make
                        # sure there is no division by zero (this can happen if
                        # there is a mixed event on no existing snowpack
                        # because a pack density has not been calculated, yet).

                    else:
                        # rsr, this should not happen, remove later
                        # if print_debug > -1:
                        #    print *, 'snow density problem', pk_depth,
                        # pk_den, pss, pkwater_equiv
                        # call print_date(1)

                        pk_den = den_max
                        pk_depth = pkwater_equiv * denmaxinv  # [inches]

                    # <
                    # Snowmelt increases by the excess free water.
                    snowmelt = snowmelt + dif_water  # [inches]

                    # Reset the previous-snowpack-plus-new-snow to the current
                    # pack water equivalent.
                    pss = pkwater_equiv  # [inches]

        # <<<
        return (
            freeh2o,
            iasw,
            pk_def,
            pk_den,
            pk_ice,
            pk_depth,
            pk_temp,
            pss,
            pst,
            snowmelt,
            pkwater_equiv,
        )

    @staticmethod
    def _calc_caloss(
        cal,
        freeh2o,
        pk_def,
        pk_ice,
        pk_temp,
        pkwater_equiv,
    ):
        """Compute change in snowpack when a net loss in heat energy has
        occurred."""

        # Loss of heat is handled differently if there is liquid water in the
        # snowpack or not.
        # 2 options below (if-then, else)
        if freeh2o < epsilon32:
            # (1) No free water exists in pack
            # Heat deficit increases because snow is colder than pack
            # (minus a negative number = plus).
            pk_def = pk_def - cal  # [cal/cm^2]

        else:
            # (2) Free water exists in pack
            # Calculate the total amount of heat per area that can be released
            # by free water freezing.
            calnd = freeh2o * 203.2  # [cal/cm^2]

            # Determine the difference between heat in free water and the heat
            # that can be absorbed by new snow (without melting). Remember
            # that cal is a negative number.
            dif = cal + calnd  # [cal/cm^2]

            # The effect of freezing water depends on whether all or only part
            # of the liquid water freezes.
            # 2 options below (if-then, else)
            if dif > zero:
                # (2) Only part of free water freezes
                # The calories absorbed by the new snow freezes some of the free water
                # (increase in ice, decrease in free water).
                pk_ice = pk_ice + (-cal / 203.2)  # [inches]
                freeh2o = freeh2o - (-cal / 203.2)  # [inches]
                return (
                    freeh2o,
                    pk_def,
                    pk_ice,
                    pk_temp,
                    pkwater_equiv,
                )

            # <
            else:  # if ( dif<=zero :
                # (1) All free water freezes
                # If all the water freezes, then the remaining heat that can be
                # absorbed by new snow (that which is not provided by freezing
                # free water) becomes the new pack heat deficit.
                if dif < zero:
                    pk_def = -dif  # [cal/cm^2]
                # Free pack water becomes ice.
                pk_ice = pk_ice + freeh2o  # [inches]
                freeh2o = zero  # [inches]

        # <<
        if pkwater_equiv > zero:
            # There is still a snowpack, calculate the new temperature.
            pk_temp = -1 * pk_def / (pkwater_equiv * 1.27)  # [degrees C]

        elif pkwater_equiv < zero:
            # if ( pkwater_equiv<-DNEARZERO ) &
            #   print *, 'snowpack issue 4, negative pkwater_equiv',
            # pkwater_equiv
            pkwater_equiv = zero

        # <
        return (
            freeh2o,
            pk_def,
            pk_ice,
            pk_temp,
            pkwater_equiv,
        )

    def snowcov(self, jj):
        """Compute snow-covered area"""

        # Local Variables: all doubles
        # snowcov_area_ante: Antecedent snow-covered area [fraction]
        # difx: Difference between the maximum snow-covered area and the
        #       snow-covered area before the last new snow [inches]
        # dify: Difference between the water equivalent before the last new
        #       snow and the previous water equivalent [inches]
        # fracy: Ratio of the unmelted amount of previous new snow in the snow
        #        pack to the value of 3/4 of previous new snow [fraction]

        # self variables
        # ai(RW), frac_swe(RW), iasw(RW), pksv(RW), scrv(RW), snowcov_area(RW),
        # snowcov_areasv(RW),

        # JLM: why is the portion of new snow (3/4) used in the depletion curve
        #      not a parameter? Or it at least seems like it would be related
        #      to self.hru_delpcrv[jj] (lower values for higher curves)

        snowcov_area_ante = self.snowcov_area[jj]

        # Reset snowcover area to the maximum
        # JLM: dont do this, it's not clear in multiple places. just save the
        # max as a local variable and use that.
        self.snowcov_area[jj] = self.snarea_curve_2d[
            self.hru_deplcrv[jj] - 1, 11 - 1
        ]  # [fraction of area]

        # Track the maximum pack water equivalent for the current snow pack.
        if self.pkwater_equiv[jj] > self.pst[jj]:
            self.pst[jj] = self.pkwater_equiv[jj]  # [inches]

        # <
        # Set ai to the maximum packwater equivalent, but no higher than the
        # threshold for complete snow cover.
        self.ai[jj] = self.pst[jj]  # [inches]
        if self.ai[jj] > self.snarea_thresh[jj]:
            self.ai[jj] = self.snarea_thresh[jj]  # [inches]

        # <
        # Calculate the ratio of the current packwater equivalent to the maximum
        # packwater equivalent for the given snowpack.
        if self.ai[jj] == zero:
            self.frac_swe[jj] = zero
        else:
            self.frac_swe[jj] = (
                self.pkwater_equiv[jj] / self.ai[jj]
            )  # [fraction]

        # <
        # There are 3 potential conditions for the snow area curve:
        # A. snow is accumulating and the pack is currently at its maximum level.
        # B. snow is depleting and the area is determined by the snow area curve.
        # C. new snow has occured on a depleting pack, temporarily resetting to 100% cover.
        # For case (C), the snow covered area is linearly interpolated between 100%
        # and the snow covered area before the new snow.
        # In general, 1/4 of the new snow has to melt before the snow covered area
        # goes below 100%, and then the remaining 3/4 has to melt to return to the
        # previous snow covered area.

        # First, the code decides whether snow is accumulating (A) or not (B/C).
        # 2 options below (if-then, else)
        if self.pkwater_equiv[jj] >= self.ai[jj]:
            # (1) The pack water equivalent is at the maximum
            # Stay on the snow area curve (it will be at the maximum because the pack
            # water equivalent is equal to ai and it can't be higher).
            # iasw = 0
            self.iasw[jj] = False

        else:
            # (2) The pack water equivalent is less than the maximum
            # If the snowpack isn't accumulating to a new maximum, it is either on the
            # curve (condition B above) or being interpolated between the previous
            # place on the curve and 100% (condition C above).
            # 2 options below (if-then, elseif)

            if self.newsnow[jj]:
                # (2.1) There was new snow...
                # New snow will always reset the snow cover to 100%. However, different
                # states change depending  on whether the previous snow area condition
                # was on the curve or being interpolated between the curve and 100%.
                # 2 options below (if-then, else)
                # if (iasw > 0) then

                if self.iasw[jj]:
                    # (2.1.1) The snow area is being interpolated between 100%
                    #         and a previous location on the curve...
                    # The location on the interpolated line is based on how much of the
                    # new snow has melted.  Because the first 1/4 of the new snow doesn't
                    # matter, it has to keep track of the current snow pack plus 3/4 of
                    # the new snow.
                    self.scrv[jj] = self.scrv[jj] + (
                        0.75 * self.net_snow[jj]
                    )  # [inches]
                    # scrv = pkwater_equiv - (0.25D0*dble(net_snow))) # [inches]
                    # RAPCOMMENT - CHANGED TO INCREMENT THE SCRV VALUE if ALREADY
                    #             INTERPOLATING BETWEEN CURVE AND 100%
                    # JLM: why NOT use pkwater_equiv on the RHS? it makes scrv
                    # appear prognostic and the logic is complicated enough
                    # to be uncertain if it is/should be drifting from
                    # pkwater_equiv

                else:
                    # (2.1.2) The current snow area is on the curve...
                    # If switching from the snow area curve to interpolation between the
                    # curve and 100%, the current state of the snow pack has to be saved
                    # so that the interpolation can continue until back to the original
                    # conditions.
                    # First, set the flag to indicate interpolation between 100% and the
                    # previous area should be done.
                    # iasw = 1  # [flag]
                    self.iasw[jj] = True  # [flag]

                    # Save the current snow covered area (before the new net snow).
                    self.snowcov_areasv[
                        jj
                    ] = snowcov_area_ante  # [inches] PAN: this is [fraction]

                    # Save the current pack water equivalent (before the new net snow).
                    self.pksv[jj] = (
                        self.pkwater_equiv[jj] - self.net_snow[jj]
                    )  # [inches]

                    # The location on the interpolated line is based on how much of the
                    # new snow has melted.  Because the first 1/4 of the new snow doesn't
                    # matter, it has to keep track of the current snow pack plus 3/4 of
                    # the new snow.
                    self.scrv[jj] = self.pkwater_equiv[jj] - (
                        0.25 * self.net_snow[jj]
                    )  # [inches]

                # <
                # The subroutine terminates here because the snow covered area always
                # starts at 100% if there is any new snow (no need to reset it from the
                # maximum value set at the beginning of the subroutine).
                return

            # <
            elif self.iasw[jj]:
                # (2.2) There was no new snow, but the snow covered area is currently
                #       being interpolated between 100% from a previous new snow and the
                #       snow covered area before that previous new snow...
                # JLM: not handling the case of no newsnow but being ON the
                #      curve is not great, apparently it's just what happens
                #      after all conditions

                # If the first 1/4 of the previous new snow has not melted yet, then the
                # snow covered area is still 100% and the subroutine can terminate.
                if self.pkwater_equiv[jj] > self.scrv[jj]:
                    return

                # <
                # At this point, the program is almost sure it is interpolating between
                # the previous snow covered area and 100%, but it is possible that
                # enough snow has melted to return to the snow covered area curve instead.
                # 2 options below (if-then, else)
                if self.pkwater_equiv[jj] >= self.pksv[jj]:

                    # (2.2.1) The snow pack still has a larger water equivalent than before
                    #         the previous new snow.  I.e., new snow has not melted back to
                    #         original area...

                    # Do the interpolation between 100% and the snow covered area before
                    # the previous new snow.

                    # Calculate the difference between the maximum snow covered area
                    # (remember that snowcov_area is always set to the maximum value at
                    # this point) and the snow covered area before the last new snow.
                    # JLM: use snowcov_frac_max instead of relying on the max being
                    #      temporarily set on the variable for clarity
                    # JLM: difx and dify are misleading, use better names
                    #      i'd assume x is swe/pkwater_equiv and y is
                    #      the associate sca, but apparently they are just
                    #      dummy names
                    difx = self.snowcov_area[jj] - self.snowcov_areasv[jj]

                    # Calculate the difference between the water equivalent before the
                    # last new snow and the previous water equivalent plus 3/4 of the last
                    # new snow. In effect, get the value of 3/4 of the previous new snow.
                    dify = self.scrv[jj] - self.pksv[jj]  # [inches]   #gl1098

                    # If 3/4 of the previous new snow is significantly different from
                    # zero, then calculate the ratio of the unmelted amount of previous
                    # new snow in the snow pack to the value of 3/4 of previous new snow.
                    # In effect, this is the fraction of the previous new snow that
                    # determines the current interpolation of snow covered area.
                    fracy = zero  # [fraction]   #gl1098
                    if dify > zero:
                        fracy = (
                            self.pkwater_equiv[jj] - self.pksv[jj]
                        ) / dify  # [fraction]

                    # <
                    # Linearly interpolate the new snow covered area.
                    self.snowcov_area[jj] = (
                        self.snowcov_areasv[jj] + fracy * difx
                    )  # [fraction of area]

                    # Terminate the subroutine.
                    return

                else:
                    # (2.2.2) The snow pack has returned to the snow water equivalent before
                    #         the previous new snow. I.e. back to original area before new snow.

                    # Reset the flag to use the snow area curve
                    # iasw = 0  # [flag]
                    self.iasw[jj] = False  # [flag]

            # <<
            # If this subroutine is still running at this point, then the program
            # knows that the snow covered area needs to be adjusted according to the
            # snow covered area curve.  So at this point it must interpolate between
            # points on the snow covered area curve (not the same as interpolating
            # between 100% and the previous spot on the snow area depletion curve).
            # JLM: better to just call this explicitly above, with each regime
            #      and make a case for no new snow and not interpolating?
            #      could also make this a function...
            self.snowcov_area[jj] = self.sca_deplcrv(
                self.snarea_curve_2d[self.hru_deplcrv[jj] - 1, :],
                self.frac_swe[jj],
            )

        # <
        return

    def snalbedo(self, jj):
        """Compute snowpack albedo"""

        # Local Variables
        # l: Number of days (or effective days) since last snowfall [days]

        # self variables
        # albedo, int_alb, iso, lst, slst, snsv,

        # The albedo is always reset to a new initial (high) value when there is
        # new snow above a threshold (parameter). Albedo is then a function of the
        # number of days since the last new snow. Intermediate conditions apply
        # when there is new snow below the threshold to reset the albedo to its
        # highest value.
        # The curve for albedo change (decreasing) is different for the snow
        # accumulation season and the snow melt season.
        # The albedo first depends on if there is no new snow during the current
        # time step, if there is new snow during accumulation season, or if there
        # is new snow during melt season.

        # 3 options below (if-then, elseif, else)
        if not self.newsnow[jj]:
            # (1) There is no new snow

            # If no new snow, check if there was previous new snow that
            # was not sufficient to reset the albedo (lst=1)
            # lst can only be greater than 0 during melt season (see below)
            # if (lst > 0) then
            if self.lst[jj]:
                # slst is the number of days (float) since the last new snowfall.
                # Set the albedo curve back three days from the number of days since
                # the previous snowfall (see salb assignment below).
                # (note that "shallow new snow" indicates new snow that is insufficient
                # to completely reset the albedo curve)
                # In effect, a shallow new snow sets the albedo curve back a few days,
                # rather than resetting it entirely.
                self.slst[jj] = self.salb[jj] - 3.0  # [days]

                # Make sure the number of days since last new snow isn't less than 1.
                if self.slst[jj] < one:
                    self.slst[jj] = one  # [days]

                # <
                if self.iso[jj] != 2:
                    # If not in melt season.

                    # NOTE: This code is unreachable in its current state. This code
                    # is only run during melt season due to the fact that lst can only
                    # be set to 1 in the melt season. Therefore, iso is always going to
                    # be equal to 2.
                    # Make sure the maximum point on the albedo curve is 5.
                    # In effect, if there is any new snow, the albedo can only get so
                    # low in accumulation season, even if the new snow is insufficient
                    # to reset albedo entirely.
                    if self.slst[jj] > 5.0:
                        self.slst[jj] = 5.0  # [days]

                # <
                # Reset the shallow new snow flag and cumulative shallow snow variable (see below).
                # lst = 0  # [flag]
                self.lst[jj] = False  # [flag]
                self.snsv[jj] = zero  # [inches]

        # <<
        elif self.iso[jj] == 2:
            # (2) New snow during the melt season
            # RAPCOMMENT - CHANGED TO ISO FROM MSO

            # If there is too much rain in a precipitation mix, albedo will not be reset.
            # New snow changes albedo only if the fraction rain is less than the
            # threshold above which albedo is not reset.
            if self.prmx[jj] < self.albset_rnm:
                # If the fraction rain doesn't prevent the albedo from being reset,
                # then how the albedo changes depends on whether the snow amount is
                # above or below the threshold for resetting albedo.
                # 2 options below (if-then, else)

                if self.net_snow[jj] > self.albset_snm:
                    # (2.1) If there is enough new snow to reset the albedo
                    # Reset number of days since last new snow to 0.
                    self.slst[jj] = zero  # [days]
                    # lst = 0  # [flag]
                    self.lst[jj] = False  # [flag]
                    # Reset the saved new snow to 0.
                    self.snsv[jj] = zero  # [inches]

                else:
                    # (2.2) If there is not enough new snow this time period to reset the
                    #       albedo on its own.

                    # snsv tracks the amount of snow that has fallen as long as the
                    # total new snow is not enough to reset the albedo.
                    self.snsv[jj] = (
                        self.snsv[jj] + self.net_snow[jj]
                    )  # [inches]

                    # Even if the new snow during this time period is insufficient to
                    # reset the albedo, it may still reset the albedo if it adds enough
                    # to previous shallow snow accumulation.  The change in albedo
                    # depends on if the total amount of accumulated shallow snow has
                    # become enough to reset the albedo or not.
                    # 2 options below (if-then, else)

                    if self.snsv[jj] > self.albset_snm:
                        # (2.2.1) If accumulated shallow snow is enough to reset the albedo.
                        # Reset the albedo states.
                        self.slst[jj] = zero  # [days]
                        # lst = 0  # [flag]
                        self.lst[jj] = False  # [flag]
                        self.snsv[jj] = zero  # [inches]

                    else:
                        # (2.2.2) If the accumulated shallow snow is not enough to
                        #         reset the albedo curve.

                        # salb records the number of days since the last new snow
                        # that reset albedo.
                        # if (lst == 0) then
                        if not self.lst[jj]:
                            self.salb[jj] = self.slst[jj]  # [days]

                        # <
                        # Reset the number of days since new snow
                        self.slst[jj] = zero  # [days]

                        # Set the flag indicating that there is shallow new snow
                        # (i.e. not enough new snow to reset albedo).
                        # lst = 1  # [flag]
                        self.lst[jj] = True  # [flag]

        # <<<<
        else:
            # (3) New snow during the accumulation season.
            # The change in albedo depends on if the precipitation is a mix, if the
            # rain is above a threshold,  or if the snow is above a threshold.
            # 4 options below (if-then, elseif, elseif, else)

            if self.pptmix[jj] < one:
                # (3.1) This is not a mixed event...
                # During the accumulation season, the threshold for resetting the
                # albedo does not apply if there is a snow-only event. Therefore, no
                # matter how little snow there is, it will always reset the albedo
                # curve the the maximum, if it occurs during the accumulation season.
                # Reset the time since last snow to 0.
                self.slst[jj] = zero  # [days]
                # There is no new shallow snow
                # lst = 0  # [flag]
                self.lst[jj] = False  # [flag]

            elif self.prmx[jj] >= self.albset_rna:
                # (3.2) This is a mixed event and the fraction rain is above
                #       the threshold above which albedo is not reset...
                # There is no new shallow snow.
                # lst = 0  # [flag]
                self.lst[jj] = False  # [flag]
                # Albedo continues to decrease on the curve

            elif self.net_snow[jj] >= self.albset_sna:
                # (3.3) If it is a mixed event and there is enough new snow to reset albedo...
                # Reset the albedo.
                self.slst[jj] = zero  # [days]
                # There is no new shallow snow.
                # lst = 0  # [flag]
                self.lst[jj] = False  # [flag]

            else:
                # (3.4) This is a mixed event and the new snow was not enough to reset the albedo...
                # Set the albedo curve back 3 days (increasing the albedo).
                self.slst[jj] = self.slst[jj] - 3.0  # [days]

                # Make sure the number of days since last new snow is not less than 0.
                if self.slst[jj] < zero:
                    self.slst[jj] = zero  # [days]

                # <
                # Make sure the number of days since last new snow is not greater than 5.
                # In effect, if there is any new snow, the albedo can only get so low
                # in accumulation season, even if the new snow is insufficient to
                # reset albedo entirely
                if self.slst[jj] > 5.0:
                    self.slst[jj] = 5.0  # [days]

                # <
                # lst = 0  # [flag]
                self.lst[jj] = False  # [flag]

            # <
            self.snsv[jj] = zero  # [inches]

        # <
        # At this point, the subroutine knows where on the curve the albedo should
        # be based on current conditions and the new snow (determined by value
        # of slst variable).

        # Get the integer value for days (or effective days) since last snowfall.
        ll = int(self.slst[jj] + 0.5)  # [days]

        # Increment the state variable for days since the last snowfall.
        self.slst[jj] = self.slst[jj] + 1.0  # [days]

        # ******Compute albedo
        # Albedo will only be different from the max (default value) if it has
        # been more than 0 days since the last new snow capable of resetting the
        # albedo.  If albedo is at the maximum, the maximum is different for
        # accumulation and melt season.
        # 3 options below (if-then, elseif, else)

        if ll > 0:
            # (1) It has been more than 0 days since the last new snow.
            # Albedo depends on whether it is currently on the accumulation season
            # curve or on the melt season curve.
            # 3 options below (if-then, elseif, else)

            if self.int_alb[jj] == 2:
                # (1.1) Currently using the melt season curve (Old snow - Spring melt period)...
                # Don't go past the last possible albedo value.
                if ll > maxalb:
                    ll = maxalb  # [days]

                # <
                # Get the albedo number from the melt season curve.
                self.albedo[jj] = amlt_init[ll - 1]  # [fraction of radiation]

            elif ll <= maxalb:
                # (1.2) Currently using the accumulation season curve (Old snow - Winter accumulation period)...
                #       and not past the maximum curve index.
                # Get the albedo number from the accumulation season curve.
                self.albedo[jj] = acum_init[ll - 1]  # [fraction of radiation]

            else:
                # (1.3) Currently using the accumulation season curve and past the maximum curve index...
                # Start using the the MELT season curve at 12 days previous to the
                # current number of days since the last new snow.
                ll = ll - 12  # [days]
                # Keep using the melt season curve until its minimum value
                # (maximum index) is reached or until there is new snow.
                if ll > maxalb:
                    ll = maxalb  # [days]

                # <
                # Get the albedo value from the melt season curve.
                self.albedo[jj] = amlt_init[ll - 1]  # [fraction of radiation]

        # <<
        elif self.iso[jj] == 2:
            # (2) New snow has reset the albedo and it is melt season.
            # Set albedo to initial value during melt season.
            # NOTE: RAPCOMMENT - CHANGED TO ISO FROM MSO
            # albedo = 0.81  # [fraction of radiation] original value
            # [fraction of radiation] value Rob suggested
            self.albedo[jj] = 0.72
            # int_alb is a flag to indicate use of the melt season curve (2)
            # or accumulation season curve (1).
            # Set flag to indicate melt season curve.
            self.int_alb[jj] = 2  # [flag]

        else:
            # (3) New snow has reset the albedo and it is accumulation season.
            # Set albedo to initial value during accumulation season.
            self.albedo[jj] = 0.91  # [fraction of radiation]
            # Set flag to indicate accumulation season curve.
            self.int_alb[jj] = 1  # [flag]

        # <
        return

    def step_4(self, jj, trd):
        """SNOWPACK RADIATION FLUXES ENERGY BALANCE"""
        # JLM:  This can/should probably be refactored to be radiation and energy balances
        # separately but for sakes of debuggin i'm keeing it as it is.

        # Set the emissivity of the air to the emissivity when there is no precipitation
        emis = self.emis_noppt[jj]  # [fraction of radiation]

        # If there is any precipitation in the HRU, reset the emissivity to 1
        if self.hru_ppt[jj] > zero:
            emis = one  # [fraction of radiation]

        # <
        # Save the current value of emissivity
        esv = emis  # [fraction of radiation]

        # The incoming shortwave radiation is the HRU radiation adjusted by the
        # albedo (some is reflected back into the atmoshphere) and the
        # transmission coefficient (some is intercepted by the winter vegetative canopy)
        swn = (
            self.swrad[jj] * (one - self.albedo[jj]) * self.rad_trncf[jj]
        )  # [cal/cm^2] or [Langleys]

        # Set the convection-condensation for a half-day interval
        cec = (
            self.cecn_coef[self.control.current_month - 1, jj] * 0.5
        )  # [cal/(cm^2 degC)] or [Langleys/degC]

        # If the land cover is trees, reduce the convection-condensation
        # parameter by half
        if self.cov_type[jj] > 2:
            # RSR: cov_type==4 is valid for trees (coniferous)
            cec = cec * 0.5  # [cal/(cm^2 degC)] or [Langleys/degC]

        # <
        # Calculate the new snow depth (Riley et al. 1973)
        # RSR: the following 3 lines of code were developed by Rob Payn, 7/10/2013
        # The snow depth depends on the previous snow pack water equivalent plus
        # the current net snow.
        self.pss[jj] = self.pss[jj] + self.net_snow[jj]  # [inches]
        dpt_before_settle = (
            self.pk_depth[jj] + self.net_snow[jj] * self.deninv[jj]
        )
        dpt1 = dpt_before_settle + self.settle_const * (
            (self.pss[jj] * self.denmaxinv) - dpt_before_settle
        )

        # RAPCOMMENT - CHANGED TO THE APPROPRIATE FINITE DIFFERENCE APPROXIMATION OF SNOW DEPTH
        # JLM: pk_depth is prognostic here
        self.pk_depth[jj] = dpt1  # [inches]

        # Calculate the snowpack density
        if dpt1 > zero:
            self.pk_den[jj] = self.pkwater_equiv[jj] / dpt1
        else:
            self.pk_den[jj] = zero  # [inch water equiv / inch depth]

        # <
        # The effective thermal conductivity is approximated (empirically)
        # as zero077 times (snowpack density)^2 [cal / (sec g degC)] Therefore,
        # the effective conductivity term (inside the square root) in the
        # equation for conductive heat exchange can be calculated as follows:
        #   (zero077 * pk_den^2) / (pk_den * 0.5)
        # where 0.5 is the specific heat of ice [cal / (g degC)]
        # this simplifies to the following
        effk = 0.0154 * self.pk_den[jj]  # [unitless]

        # 13751 is the number of seconds in 12 hours over pi
        # So for a half day, to calculate the conductive heat exchange per cm
        # of snow per cm^2 area per degree temperature difference is the following
        # In effect, multiplying cst times the temperature gradient gives the
        # heatexchange by heat conducted (calories) per square cm of snowpack
        cst = self.pk_den[jj] * (
            np.sqrt(effk * 13751.0)
        )  # [cal/(cm^2 degC)] or [Langleys / degC]

        # Check whether to force spring melt
        # Spring melt is forced if time is before the melt-force day and after
        # the melt-look day (parameters). If between these dates, the spring melt
        # applies if the snowpack temperature is above or equal to 0 for more
        # than 4 cycles of the snorun function.
        if self.iso[jj] == 1:
            # Before the first melt-force day

            if self.mso[jj] == 2:
                # After the first melt-look day
                # Melt season is determined by the number of days the snowpack is
                # above 0 degrees C. The first time that the snowpack is isothermal
                # at 0 degrees C for more than 4 days is the beginning of snowmelt season.
                # 2 options below (if-then, else)

                if self.pk_temp[jj] >= zero:
                    # (1) The snowpack temperature is 0 degrees
                    # Increment the number of days that the snowpack has been
                    # isothermal at 0 degrees C
                    self.lso[jj] = self.lso[jj] + 1  # [days]

                    # If the snowpack temperature has been 0 or greater for more than 4 cycles
                    if self.lso[jj] > 4:
                        # Set the melt-force flag and reset counter
                        self.iso[jj] = 2  # [flag]
                        self.lso[jj] = 0  # [days]

                # <<
                else:
                    # (2) The snowpack temperature is less than 0 degrees
                    # Reset the counter for days snowpack temperature is above 0
                    self.lso[jj] = 0  # [days]

        # <<<
        # Compute energy balance for night period
        # niteda is a flag indicating nighttime (1) or daytime (2)
        # Set the flag indicating night time
        niteda = 1  # [flag]

        # No shortwave (solar) radiation at night.
        sw = zero  # [cal / cm^2] or [Langleys]

        # Temperature is halfway between the minimum and average temperature
        # for the day.
        temp = (self.tminc[jj] + self.tavgc[jj]) * 0.5

        # Track total heat flux from both night and day periods

        # Calculate the night time energy balance
        (
            self.tcal[jj],
            self.freeh2o[jj],
            self.pk_def[jj],
            self.pk_ice[jj],
            self.pk_temp[jj],
            self.pkwater_equiv[jj],
        ) = self._calc_snowbal(
            niteda=niteda,
            cec=cec,
            cst=cst,
            esv=esv,
            sw=sw,
            temp=temp,
            trd=trd,
            calin=self._calc_calin,  # fn
            caloss=self._calc_caloss,  # fn
            canopy_covden=self.canopy_covden[jj],
            den_max=self.den_max,
            denmaxinv=self.denmaxinv,
            emis_noppt=self.emis_noppt[jj],
            freeh2o=self.freeh2o[jj],
            freeh2o_cap=self.freeh2o_cap[jj],
            hru_ppt=self.hru_ppt[jj],
            iasw=self.iasw[jj],
            pk_def=self.pk_def[jj],
            pk_den=self.pk_den[jj],
            pk_depth=self.pk_depth[jj],
            pk_ice=self.pk_ice[jj],
            pk_temp=self.pk_temp[jj],
            pkwater_equiv=self.pkwater_equiv[jj],
            pss=self.pss[jj],
            pst=self.pst[jj],
            snowcov_area=self.snowcov_area[jj],
            snowmelt=self.snowmelt[jj],
            tcal=self.tcal[jj],
            tstorm_mo=self.tstorm_mo[self.control.current_month - 1, jj],
        )

        # [cal/cm^2] or [Langleys]

        # Compute energy balance for day period (if the snowpack still exists)
        # THIS SHOULD HAPPEN IN SNOBAL
        if self.pkwater_equiv[jj] > zero:
            # Set the flag indicating daytime
            niteda = 2  # [flag]

            # Set shortwave radiation as calculated earlier
            sw = swn  # [cal/cm^2] or [Langleys]

            # Temperature is halfway between the maximum and average temperature
            # for the day.
            temp = (self.tmaxc[jj] + self.tavgc[jj]) * 0.5  # [degrees C]
            (
                cals,
                self.freeh2o[jj],
                self.pk_def[jj],
                self.pk_ice[jj],
                self.pk_temp[jj],
                self.pkwater_equiv[jj],
            ) = self._calc_snowbal(
                niteda=niteda,
                cec=cec,
                cst=cst,
                esv=esv,
                sw=sw,
                temp=temp,
                trd=trd,
                calin=self._calc_calin,  # fn
                caloss=self._calc_caloss,  # fn
                canopy_covden=self.canopy_covden[jj],
                den_max=self.den_max,
                denmaxinv=self.denmaxinv,
                emis_noppt=self.emis_noppt[jj],
                freeh2o=self.freeh2o[jj],
                freeh2o_cap=self.freeh2o_cap[jj],
                hru_ppt=self.hru_ppt[jj],
                iasw=self.iasw[jj],
                pk_def=self.pk_def[jj],
                pk_den=self.pk_den[jj],
                pk_depth=self.pk_depth[jj],
                pk_ice=self.pk_ice[jj],
                pk_temp=self.pk_temp[jj],
                pkwater_equiv=self.pkwater_equiv[jj],
                pss=self.pss[jj],
                pst=self.pst[jj],
                snowcov_area=self.snowcov_area[jj],
                snowmelt=self.snowmelt[jj],
                tcal=self.tcal[jj],
                tstorm_mo=self.tstorm_mo[self.control.current_month - 1, jj],
            )

            # Track total heat flux from both night and day periods
            self.tcal[jj] = self.tcal[jj] + cals  # [cal/cm^2] or [Langleys]

        # <
        return

    @staticmethod
    def _calc_snowbal(
        niteda,
        cec,
        cst,
        esv,
        sw,
        temp,
        trd,
        #
        calin,  # fn
        caloss,  # fn
        canopy_covden,
        den_max,
        denmaxinv,
        emis_noppt,
        freeh2o,
        freeh2o_cap,
        hru_ppt,
        iasw,
        pk_def,
        pk_den,
        pk_depth,
        pk_ice,
        pk_temp,
        pkwater_equiv,
        pss,
        pst,
        snowcov_area,
        snowmelt,
        tcal,
        tstorm_mo,
    ) -> np.ndarray:
        """Snowpack energy balance: 1st call is for night period, 2nd call
        for day period."""

        # *******************************************************************
        # Calculate the potential long wave energy from air based on
        # temperature (assuming perfect black-body emission).
        # Stefan Boltzmann/2 = (11.71E-8)/2 = 0.585E-7
        # because add for day and night
        air = 0.585e-7 * ((temp + 273.16) ** 4.0)  # [cal/cm^2] or [Langleys]

        # Set emissivity, which is the fraction of perfect black-body emission
        # that is actually applied.
        emis = esv  # [fraction of radiation]
        # JLM this is bonkers, it's just the inverse of above (setp 4). not
        # clear we need esv could just have emis_local if we dont want to
        # track that. this should be combined with setting emis below as and
        # else condition?

        # The snowpack surface temperature and long-wave radiation FROM the
        # snowpack depend on the air temperature (effectively, snowpack
        # temperature cannot be larger than 0 degC).

        # 2 options below (if-then, else)
        if temp < zero:
            # (1) If the temperature is below freezing, surface snow
            # temperature and long wave energy are determined by temperature...
            ts = temp  # [degrees C]
            sno = air  # [cal/cm^2] or [Langleys]

        else:
            # (2) If the temperature is at or above freezing, snow temperature
            # and long wave energy are set to values corresponding to a
            # temperature of 0 degC...
            ts = zero  # [degrees C]
            sno = 325.7  # [cal/cm^2] or [Langleys]

        # <
        if hru_ppt > zero:
            # If precipitation over the time period was due to convective
            # thunderstorms, then the emissivity should be reset.

            if tstorm_mo == 1:
                # The emissivity of air depends on if it is day or night and
                # the fraction of measured short wave radiation to potential
                # short wave radiation is used as a surrogate to the duration
                # of the convective storms.
                # 2 options below (if-then, else)

                if niteda == 1:
                    # (1) Night
                    # Set the default emissivity
                    emis = 0.85  # [fraction of radiation]

                    # If measured radiation is greater than 1/3 potential
                    # radiation through the time period, then the emissivity
                    # is set to the "no precipitation" value.
                    if trd > ONETHIRD:
                        emis = emis_noppt  # [fraction of radiation]

                # <<
                else:
                    # (2) Day
                    # If measured radiation is greater than 1/3 potential
                    # radiation but less than 1/2, then the emissivity is
                    # interpolated between 1.0 and 0.85.
                    if trd > ONETHIRD:
                        emis = 1.29 - (0.882 * trd)  # [fraction of radiation]

                    # If measured radiation is greater than 1/2 potential
                    # radiation, then the emissivity is interpolated between
                    # 0.85 and 0.75.
                    if trd >= 0.5:
                        emis = 0.95 - (0.2 * trd)  # [fraction of radiation]

        # <<<<
        # Calculate the net incoming long wave radiation coming from the sky or
        # canopy in the uncovered or covered portions of the snowpack,
        # respectively. Note that the canopy is assumed to be a perfect
        # blackbody (emissivity=1) and the air has emissivity as determined
        # from previous calculations.
        sky = (one - canopy_covden) * (
            (emis * air) - sno
        )  # [cal/cm^2] or [Langleys]
        can = canopy_covden * (air - sno)  # [cal/cm^2] or [Langleys]

        # RAPCOMMENT - CHECK THE INTERECEPT MODULE FOR CHANGE. What if the land
        #              cover is grass? Is this automatically covered by
        #              canopy_covden being zero if the cover type is grass?

        # If air temperature is above 0 degC then set the energy from
        # condensation and convection, otherwise there is no energy from
        # convection or condensation.
        cecsub = zero  # [cal/cm^2] or [Langleys]

        if (temp > zero) and (hru_ppt > zero):
            cecsub = cec * temp  # [cal/cm^2] or [Langleys]

        # <
        # Total energy potentially available from atmosphere: longwave,
        # shortwave, and condensation/convection.
        cal = sky + can + cecsub + sw  # [cal/cm^2] or [Langleys]

        # If the surface temperature of the snow is 0 degC, and there is net
        # incoming energy, then energy conduction has to be from the surface
        # into the snowpack. Therefore, the energy from the atmosphere is
        # applied to the snowpack and subroutine terminates.
        if (ts >= zero) and (cal > zero):
            (
                freeh2o,
                iasw,
                pk_def,
                pk_den,
                pk_ice,
                pk_depth,
                pk_temp,
                pss,
                pst,
                snowmelt,
                pkwater_equiv,
            ) = calin(
                cal,
                den_max,
                denmaxinv,
                freeh2o,
                freeh2o_cap,
                iasw,
                pk_def,
                pk_den,
                pk_depth,
                pk_ice,
                pk_temp,
                pkwater_equiv,
                pss,
                pst,
                snowcov_area,
                snowmelt,
            )
            return (
                cal,
                freeh2o,
                pk_def,
                pk_ice,
                pk_temp,
                pkwater_equiv,
            )

        # <
        # If the program gets to this point, then either the surface
        # temperature is less than 0 degC, or the total energy from the
        # atmosphere is not providing energy to the snowpack.

        # Because the temperature of the surface of the snowpack is assumed to
        # be controlled by air temperature, there is a potential heat flux due
        # to conduction between the deeper snowpack and its surface.

        # Calculate conductive heat flux as a function of the temperature
        # gradient then set new snowpack conditions depending on the direction
        # of heat flow.
        qcond = cst * (ts - pk_temp)  # [cal/cm^2] or [Langleys]

        # RAPCOMMENT - The original equation in the paper implies that the this
        #              equation should be relative to the temperature gradient
        #              in degF, not degC (Anderson 1968).  Which is correct?

        # The energy flow depends on the direction of conduction and the
        # temperature of the surface of the snowpack. The total energy from the
        # atmosphere can only penetrate into the snow pack if the temperature
        # gradient allows conduction from the surface into the snowpack.
        # 4 options below (if-then, elseif, elseif, else)

        if qcond < zero:
            # (1) Heat is conducted from the snowpack to the surface
            #     (atmospheric energy is NOT applied to snowpack)...
            # If the temperature of the snowpack is below 0 degC,
            # add to the heat deficit.  Otherwise, remove heat
            # from the 0 degC isothermal snow pack.

            if pk_temp < zero:
                # Increase the heat deficit (minus a negative) and adjust
                # temperature.
                pk_def = pk_def - qcond  # [cal/cm^2] or [Langleys]
                pk_temp = -1 * pk_def / (pkwater_equiv * 1.27)  # [degrees C]

            else:
                # Remove heat from the snowpack.
                (freeh2o, pk_def, pk_ice, pk_temp, pkwater_equiv,) = caloss(
                    qcond,
                    freeh2o,
                    pk_def,
                    pk_ice,
                    pk_temp,
                    pkwater_equiv,
                )
                # NOTE: Even though cal is not applied to the snowpack under
                # this condition, it maintains its value and the referencing
                # code uses it to calculate the total energy balance of the
                # snowpack. Right now cal isn't used for anything outside this
                # subroutine, but care should be taken if it is.
        # <<
        elif qcond < epsilon32:
            # (2)  There is no heat conduction, qcond = zero
            # If the pack temperature is isothermal at 0 degC, then apply any
            # incoming radiation, condensation (latent heat), and convection
            # heat to the snowpack.

            if pk_temp >= zero:
                # It does not appear that the interior of the following if
                # statement is reachable in its current form, because if these
                # conditions are true, then the code for surface temperature=0
                # and cal=positive number would have run and the subroutine
                # will have terminated.
                if cal > zero:
                    (
                        freeh2o,
                        iasw,
                        pk_def,
                        pk_den,
                        pk_ice,
                        pk_depth,
                        pk_temp,
                        pss,
                        pst,
                        snowmelt,
                        pkwater_equiv,
                    ) = calin(
                        cal,
                        den_max,
                        denmaxinv,
                        freeh2o,
                        freeh2o_cap,
                        iasw,
                        pk_def,
                        pk_den,
                        pk_depth,
                        pk_ice,
                        pk_temp,
                        pkwater_equiv,
                        pss,
                        pst,
                        snowcov_area,
                        snowmelt,
                    )

        # <<<
        elif ts >= zero:
            # (3) conduction is from the surface to the snowpack and the
            #     surface temperature is 0 degrees C...
            # Note that cal must be <= 0 for this condition to apply.
            # Otherwise, the program wouldn't have gotten to this point.
            # Determine if the conductive heat is enough to overcome the
            # current heat deficit.
            pk_defsub = pk_def - qcond

            if pk_defsub < zero:
                # Deficit is overcome and snowpack becomes isothermal at 0degC.
                pk_def = zero  # [cal/cm^2] or [Langleys]
                pk_temp = zero  # [degrees C]

            else:
                # Deficit is decreased by conducted heat and temperature is
                # recalculated.
                pk_def = pk_defsub  # [cal/cm^2] or [Langleys]
                pk_temp = -pk_defsub / (pkwater_equiv * 1.27)  # [degrees C]

        # <<
        else:
            # (4) conduction is from the surface to the snowpack and the
            #     surface temperature is less than 0 degrees C...
            # Calculate the pack deficit if the snowpack was all at the surface
            # temperature, then calculate how many calories to shift the pack
            # to that deficit (pks will be a positive number because the
            # conduction direction is from the surface into the snowpack).
            pkt = -ts * (pkwater_equiv * 1.27)  # [cal/cm^2] or [Langleys]
            pks = pk_def - pkt  # [cal/cm^2] or [Langleys]

            # Determine if the conducted heat is enough to shift the pack to the
            # deficit relative to the surface temperature.
            pk_defsub = pks - qcond  # [cal/cm^2] or [Langleys]

            # The effect of incoming conducted heat depends on whether it is
            # enough to bring the snowpack to the same temperature as the
            # surface or not.
            # 2 options below (if-then, else)
            if pk_defsub < zero:
                # (4.1) There is enough conducted heat to bring the deep
                #       snowpack to the surface temperature...
                # There is enough conduction to change to the new pack deficit.
                pk_def = pkt  # [cal/cm^2] or [Langleys]
                pk_temp = ts  # [degrees C]

            else:
                # (4.2) There is not enough conducted heat to bring the deep
                #       snowpack to the surface temperature...
                # The pack deficit doesn't make it all the way to the surface
                # deficit,
                # but is decreased relative to the conducted heat. Note that
                # the next
                # statement is equivalent to pk_def = pk_def - qcond
                pk_def = pk_defsub + pkt  # [cal/cm^2] or [Langleys]
                pk_temp = -1 * pk_def / (pkwater_equiv * 1.27)  # [degrees C]

        # <
        return (
            cal,
            freeh2o,
            pk_def,
            pk_ice,
            pk_temp,
            pkwater_equiv,
        )

    @staticmethod
    def _calc_snowevap(
        freeh2o,
        hru_intcpevap,
        pk_def,
        pk_ice,
        pk_temp,
        pkwater_equiv,
        potet,
        potet_sublim,
        snow_evap,
        snowcov_area,
        verbose,
    ):
        # Local Variables
        # avail_et:
        # cal:  Amount of heat deficit that is removed by sublimating ice
        # [cal/cm^2]
        # ez: Amount of evaporation affecting the snowpack [inches]

        # The amount of evaporation affecting the snowpack is the total
        # evaporation potential minus the evaporation from the interception
        # storage.
        ez = potet_sublim * potet * snowcov_area - hru_intcpevap  # [inches]

        # The effects of evaporation depend on whether there is any potential
        # for evaporation, and if the potential evapotation is enough to
        # completely deplete the snow pack or not.
        # 3 options below (if-then, elseif, else)
        if ez < epsilon32:
            # (1) There is no potential for evaporation...
            snow_evap = 0.0  # [inches]

        elif ez >= pkwater_equiv:
            # (2) Enough potential evaporation to entirely deplete the
            # snowpack...
            # Set the evaporation to the pack water equivalent and set all
            # snowpack
            # variables to no-snowpack values.
            snow_evap = pkwater_equiv  # [inches]

            snow_evap = pkwater_equiv  # [inches]
            pkwater_equiv = zero  # [inches]
            pk_ice = zero  # [inches]
            pk_def = zero  # [cal/cm^2]
            freeh2o = zero  # [inches]
            pk_temp = zero  # [degrees C]

        else:
            # (3) Potential evaporation only partially depletes snowpack...
            # Evaporation depletes the amount of ice in the snowpack
            # (sublimation).
            pk_ice = pk_ice - ez
            # Change the pack conditions according to whether there is any
            # ice left in the snowpack.
            if pk_ice < zero:
                # RAPCOMMENT - CHANGED TO CHECK FOR NEGATIVE PACK ICE
                # If all pack ice is removed, then there cannot be a heat
                # deficit.
                pk_ice = zero
                pk_def = zero
                pk_temp = zero
            else:
                # Calculate the amount of heat deficit that is removed by the
                # sublimating ice. Note that this only changes the heat
                # deficit if the pack temperature is less than 0 degC.
                cal = pk_temp * ez * 1.27
                pk_def = pk_def + cal

            # <
            # Remove the evaporated water from the pack water equivalent.
            pkwater_equiv = pkwater_equiv - ez
            snow_evap = ez

        # <
        if snow_evap < zero:
            pkwater_equiv = pkwater_equiv - snow_evap

            if pkwater_equiv < zero:
                if verbose:
                    if pkwater_equiv < -epsilon64:
                        print(
                            "snowpack issue, negative pkwater_equiv in "
                            f"snowevap: {pkwater_equiv}"
                        )

                    #  <
                    pkwater_equiv = zero

            # <<
            snow_evap = zero

        # <
        avail_et = potet - hru_intcpevap - snow_evap
        if avail_et < zero:
            snow_evap = snow_evap + avail_et
            pkwater_equiv = pkwater_equiv - avail_et

            if snow_evap < zero:
                pkwater_equiv = pkwater_equiv - snow_evap

                if pkwater_equiv < zero:
                    if verbose:
                        if pkwater_equiv < -epsilon64:
                            print(
                                "snowpack issue 2, negative pkwater_equiv in "
                                f"snowevap: {pkwater_equiv}"
                            )

                    # <<
                    # To be sure negative snowpack is ignored
                    pkwater_equiv = zero

                # <
                snow_evap = zero

        # <<
        return (
            freeh2o,
            pk_def,
            pk_ice,
            pk_temp,
            pkwater_equiv,
            snow_evap,
        )

    @staticmethod
    def set_snow_zero():
        pkwater_equiv = zero
        pk_depth = zero
        pss = zero
        snsv = zero
        lst = False
        pst = zero
        iasw = False
        albedo = zero
        pk_den = zero
        snowcov_area = zero
        pk_def = zero
        pk_temp = zero
        pk_ice = zero
        freeh2o = zero
        snowcov_areasv = zero
        ai = zero
        frac_swe = zero

        return (
            pkwater_equiv,
            pk_depth,
            pss,
            snsv,
            lst,
            pst,
            iasw,
            albedo,
            pk_den,
            snowcov_area,
            pk_def,
            pk_temp,
            pk_ice,
            freeh2o,
            snowcov_areasv,
            ai,
            frac_swe,
        )
