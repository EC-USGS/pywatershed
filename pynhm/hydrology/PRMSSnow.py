from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control
from ..constants import HruType, one, zero, epsilon, inch2cm

adaptable = Union[str, np.ndarray, Adapter]

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


class PRMSSnow(StorageUnit):
    """PRMS snow pack

    Args:
        params: parameter object
        atm: atmosphere object

    """

    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        orad_hru: adaptable,
        soltab_horad_potsw: adaptable,
        swrad: adaptable,
        hru_ppt: adaptable,
        prmx: adaptable,
        tavgc: adaptable,
        tmaxc: adaptable,
        tminc: adaptable,
        net_ppt: adaptable,
        net_rain: adaptable,
        net_snow: adaptable,
        transp_on: adaptable,
        verbose: bool = False,
    ) -> "PRMSSnow":

        self.name = "PRMSSnow"
        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
            subclass_name=self.name,
        )

        # Adapt every input
        self._input_variables_dict = {}
        for ii in self.inputs:
            self._input_variables_dict[ii] = adapter_factory(locals()[ii], ii)

        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get snow pack parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nhru",
            "ndeplval",
            "nmonths",
            # "cov_type",
            "hru_type",
            # "hru_route_order",  # really necessary? does it matter for column calcs?
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
            "net_ppt",
            "net_rain",
            "net_snow",
            "orad_hru",
            "prmx",
            "soltab_horad_potsw",
            "swrad",
            "tavgc",
            "tmaxc",
            "tminc",
            "transp_on",
        )

    @staticmethod
    def get_variables() -> tuple:
        """Get snow pack variables

        Returns:
            variables: variables

        """
        return (
            "ai",
            "albedo",
            "frac_swe",
            "freeh2o",
            "iasw",
            "int_alb",
            "iso",
            "lso",
            "lst",
            "mso",
            "pk_def",
            "pk_den",
            "pk_depth",
            "pk_ice",
            "pk_precip",
            "pk_temp",
            "pksv",
            "pkwater_ante",
            "pkwater_equiv",
            "pptmix_nopack",
            "pss",
            "pst",
            "salb",
            "scrv",
            "slst",
            "snow_evap",
            "snowcov_area",
            "snowcov_areasv",
            "snowmelt",
            "snsv",
            "tcal",
            # "newsnow",  # precipitation module these are just flags
            # "pptmix",  # precipitation module
            # acum, amlt, deninv, denmaxinv are not worth tracking on self as in PRMS6
        )

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

    @staticmethod
    def get_init_values() -> dict:
        """Get snow pack inital values

        Returns:
            dict: inital values for named variables
        """

        return {
            "ai": zero,
            "albedo": zero,
            "frac_swe": zero,
            "freeh2o": zero,
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
            "pk_precip": zero,
            "pk_temp": zero,
            "pksv": zero,
            "pptmix_nopack": False,
            "salb": zero,
            "scrv": zero,
            "slst": zero,
            "snow_evap": zero,
            "snowcov_area": zero,
            "snowcov_areasv": zero,
            "snowmelt": zero,
            "snsv": zero,
            "tcal": zero,
        }

    def set_initial_conditions(self):
        """Initialize PRMSSnow snowpack variables."""

        # Deninv and denmaxinv not in variables nor in metadata but we can set them on self
        # for convenience
        self.deninv = one / self.den_init.copy()
        self.denmaxinv = one / self.den_max.copy()

        self.pkwater_equiv = self.snowpack_init.copy()
        self.pkwater_equiv_ante = None

        sd = int(self.ndeplval / 11)
        self.snarea_curve_2d = np.reshape(self.snarea_curve, (sd, 11))

        if self.control.config["init_vars_from_file"] in [0, 2, 3]:

            # The super().__init__ already set_initial_conditions using its set_inital_conditions
            # Below Im just following PRMS6, will reconcile later with the super (may be redundant).
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
                self.initialize_var(vv)

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
                        1:11, self.hru_deplcrv[wh_pkweq_gt_zero[ww]]
                    ],
                    self.frac_swe[wh_pkweq_gt_zero[ww]],
                )

            self.pss = self.pkwater_equiv
            self.pst = self.pkwater_equiv

        else:

            raise RuntimeError("Snow restart capability not implemented")
            # JLM: a list of restart variables dosent shed light on what states actually have memory.
            # JLM: could there be a diagnostic part of the advance?

        return

    def advance(self) -> None:
        """Advance the snow pack
        Returns:
            None

        """
        self.pkwater_ante = self.pkwater_equiv
        self._itime_step += 1

        for key, value in self._input_variables_dict.items():
            # JLM: This is only because adapter advances dont all take current time.
            if key == "soltab_horad_potsw":
                value.advance(self.control.current_time)
            else:
                value.advance()
            v = getattr(self, key)
            v[:] = value.current

        return

    def calculate(self, simulation_time):
        """Calculate snow pack terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None
        """

        self._simulation_time = simulation_time

        # curr_month => model_time%Nowmonth, &
        # day_of_year => model_time%day_of_year, &
        # day_of_water_year => model_time%day_of_water_year, &

        # nhru => model_basin%nhru, &
        # nmonths => model_basin%nmonths, &
        # active_hrus => model_basin%active_hrus, &
        # active_mask => model_basin%active_mask, &
        # cov_type => model_basin%cov_type, &
        # hru_type => model_basin%hru_type, &
        # hru_route_order => model_basin%hru_route_order, &

        cals = zero

        # these dont really need to be set on self
        newsnow = np.full([self.nhru], False, dtype=bool)
        pptmix = np.zeros([self.nhru])

        net_snow_gt_zero = self.net_snow > zero
        wh_net_snow_gt_zero = np.where(net_snow_gt_zero)
        newsnow[wh_net_snow_gt_zero] = True

        net_rain_gt_zero = self.net_rain > zero
        wh_net_snow_gt_zero_and_net_rain_gt_zero = np.where(
            net_rain_gt_zero & net_rain_gt_zero
        )
        pptmix[wh_net_snow_gt_zero_and_net_rain_gt_zero] = 1

        pptmix_nopack = False

        for jj in range(self.nhru):

            if self.hru_type[jj] == HruType.LAKE:
                continue

            trd = self.orad_hru[jj] / self.soltab_horad_potsw[jj]

            # If it's the first julian day of the water year, several
            # variables need to be reset:
            # - reset the previous snow water eqivalent plus new snow to 0
            # - reset flags to indicate it is not melt season or potetential melt season
            # - reset the counter for the number of days a snowpack is at 0 deg Celsius
            # TODO: rsr, do we want to reset all HRUs, what about Southern Hemisphere
            if self.control.current_dowy == 1:
                self.pss[jj] = zero  # [inches]
                self.iso[jj] = True  # [flag]
                self.mso[jj] = True  # [flag]
                self.lso[jj] = 0  # [counter]

            # HRU SET-UP - SET DEFAULT VALUES AND/OR BASE CONDITIONS FOR THIS TIME PERIOD
            # **************************************************************

            self.pkwater_ante[jj] = self.pkwater_equiv[jj]

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

            # If the day of the water year is beyond the first day to look for melt
            # season indicated by the parameter, then set the flag indicating to watch
            # for melt season.
            # TODO: rsr, need to rethink this at some point
            if self.control.current_doy == self.melt_look[jj]:
                self.mso[jj] = 2  # [flag]  # could emume this BEFORE/AFTER

            if self.pkwater_equiv[jj] < epsilon:
                # No existing snowpack
                if newsnow[jj]:
                    # Skip the HRU if there is no snowpack and no new snow
                    # Reset to be sure it is zero if snowpack melted on last timestep.
                    self.snowcov_area[jj] = zero
                    continue
                else:
                    # We have new snow; the initial snow-covered area is complete (1)
                    self.snowcov_area[jj] = one  # [fraction of area]

            # HRU STEP 1 - DEAL WITH PRECIPITATION AND ITS EFFECT ON THE WATER
            #              CONTENT AND HEAT CONTENT OF SNOW PACK
            # ***********************************************************************
            # WARNING: pan - wouldn't this be pkwater_equiv > DNEARZERO?
            if (
                self.pkwater_equiv[jj] > zero and self.net_ppt[jj] > zero
            ) or self.net_snow[jj] > zero:
                self.ppt_to_pack(jj)

            # HRU STEP 2 - CALCULATE THE NEW SNOW COVERED AREA
            # **********************************************************

            # HRU STEP 3 - COMPUTE THE NEW ALBEDO
            # **********************************************************

            # HRU STEP 4 - DETERMINE RADIATION FLUXES AND SNOWPACK
            #              STATES NECESSARY FOR ENERGY BALANCE
            # **********************************************************

            #  HRU STEP 5 - CALCULATE SNOWPACK LOSS TO EVAPORATION
            # ********************************************************

            #  HRU CLEAN-UP - ADJUST FINAL HRU SNOWPACK STATES AND
            #                 INCREMENT THE BASIN TOTALS
            # *********************************************************

        return

    @staticmethod
    def sca_deplcrv(snarea_curve: np.ndarray, frac_swe: float) -> float:
        """Interpolate along snow covered area depletion curve"""
        if frac_swe > one:
            res = snarea_curve[-1]
        else:
            # Get the indices (as integers) of the depletion curve that bracket the
            # given frac_swe (next highest and next lowest).
            idx = int(10.0 * (frac_swe + 0.2)) - 1  # [index]
            jdx = idx - 1 - 1  # [index]
            if idx > (11 - 1):
                idx = 11 - 1
            # Calculate the fraction of the distance (from the next lowest) the given
            # frac_swe is between the next highest and lowest curve values.
            dify = (frac_swe * 10.0) - float(jdx - 1)  # [fraction]
            # Calculate the difference in snow covered area represented by next
            # highest and lowest curve values.
            difx = snarea_curve[idx] - snarea_curve[jdx]
            # Linearly interpolate a snow covered area between those represented by
            # the next highest and lowest curve values.
            res = snarea_curve[jdx] + dify * difx

        return res

    def ppt_to_pack(jj):
        """Add rain and/or snow to snowpack."""
        # map from prms6 to here
        # model_climate holds self.pkwater_equiv
        # model_precip holds tmax_allsnow_c  ---   ???? where is this from

        caln = zero
        calpr = zero
        calps = zero
        pndz = zero

        tsnow = tavg[jj]  # [degrees C]
        month = self.control.current_month
        month_ind = month - 1

        if self.pptmix[jj] == 1:
            # (1) If precipitation is mixed...
            # If there is any rain, the rain temperature is halfway between the maximum
            # temperature and the allsnow temperature.
            train = (
                self.tmaxc[jj] + tmax_allsnow[jj, month_ind]
            ) * 0.5  # [degrees C]

            # Temperatures will differ, depending on the presence of existing snowpack.
            # should this be epsilon?
            if self.pkwater_equiv[jj] > zero:
                # If there is a snowpack, snow temperature is halfway between the minimum
                # daily temperature and maximum temperature for which all precipitation is snow.
                tsnow = (
                    self.tminc[jj] + self.tmax_allsnow[jj, month_ind]
                ) * 0.5  # [degrees C]
            elif self.pkwater_equiv[jj] < zero:
                # If no existing snowpack, snow temperature is the average temperature for the day.
                self.pkwater_equiv[
                    jj
                ] = zero  # To be sure negative snowpack is ignored
        else:
            # (2) If precipitation is all snow or all rain...
            # If there is any rain, the rain temperature is the average temperature.
            train = tavgc[jj]  # [degrees C]
            if train < epsilon:
                # If average temperature is close to freezing, the rain temperature is
                # halfway between the maximum daily temperature and maximum temperature
                # for which all precipitation is snow.
                train = (
                    tmax[jj] + self.tmax_allsnow[jj, month_ind]
                ) * 0.5  # [degrees C]

        # <-  <-
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

        # -----------------------------------------------------------------------------
        # -----------------------------------------------------------------------------
        if self.pkwater_equiv[jj] > zero:
            # (1) There is net rain on an existing snowpack...
            if net_rain[jj] > zero:
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

                    if abs(self.net_rain[jj] - pndz) < epsilon:
                        # (1.1.1) Exactly enough rain to bring pack to isothermal...
                        # Heat deficit and temperature of the snowpack go to 0.
                        self.pk_def[jj] = zero  # [cal/cm^2]
                        self.pk_temp[jj] = zero  # [degrees C]

                        # In the process of giving up its heat, all the net rain freezes and
                        # becomes pack ice.
                        self.pk_ice[jj] = (
                            self.pk_ice[jj] + net_rain[jj]
                        )  # [inches]

                    elif net_rain[jj] < pndz:
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

                        self.calin(calpr, jj)

                        # if (chru == 5):
                        #     print *, '*', pkwater_equiv[jj]
                        # endif

                else:
                    # (1.2) Rain on snowpack that is isothermal at 0 degC (no heat deficit)...
                    # All net rain is added to free water in the snowpack.
                    self.freeh2o[jj] = self.freeh2o[jj] + self.net_rain[jj]

                    # Calculate the heat per area added by the rain (using specific heat of water).
                    calpr = train * net_rain[jj] * inch2cm  # [cal/cm^2]

                    # Add the new heat to the snow pack (the heat in rain will melt some
                    # of the pack ice when the water cools to 0 degC).

                    self.calin(calpr, jj)

            elif self.net_rain[jj] > zero:
                # (2) If there is net rain but no snowpack, set flag for a mix on no snowpack.

                # Be careful with the code here.
                # If this subroutine is called when there is an all-rain day on no
                # existing snowpack (currently, it will not), then the flag here will be
                # set inappropriately.
                self.pptmix_nopack[jj] = True  # [flag]

        return

    def calin(self, cal, jj):
        """Compute changes in snowpack when a net gain in heat energy has occurred."""

        # Local Variables: all doubles
        # apk_ice: Pack-ice per area [inches]
        # apmlt: Actual potential snowmelt [inches]
        # dif: Difference between incoming calories and calories needed to bring the pack to isothermal at 0 [cal/cm^2]
        # dif_water: Amount of free water in excess of the capacity to hold free water.
        # pmlt: Potential amount of snowmelt from excess heat in rain [inches]
        # pwcap: Capacity of the snowpack to hold free water [inches]

        # Calculate the difference between the incoming calories and the calories
        # needed to bring the pack to isothermal at 0 (heat deficit).
        dif = cal - self.pk_def  # [cal/cm^2]

        # The way incoming heat is handled depends on whether there is not enough,
        # just enough, or more than enough heat to overcome the heat deficit of the snowpack.
        # 3 choices below (if-then, elseif, else)

        if dif < zero:
            # (1) Not enough heat to overcome heat deficit...
            # Reduce the heat deficit by the amount of incoming calories and adjust to
            # the new temperature based on new heat deficit.
            self.pk_def = selfpk_def - cal  # [cal/cm^2]
            self.pk_temp = (
                -1 * self.pk_def / (self.pkwater_equiv * 1.27)
            )  # [degrees C]

        elif abs(dif) < epsilon:
            # JLM: The test had been equality with zero, changed to less than epsilon.
            # (2) Just enough heat to overcome heat deficit
            # rsr 1/27/2016 why not set all snow states to 0 ???
            # Set temperature and heat deficit to zero.
            self.pk_temp = zero  # [degrees C]
            self.pk_def = zero  # [cal/cm^2]

        elif dif > zero:
            # (3) More than enough heat to overcome heat deficit (melt ice)...
            # Calculate the potential amount of snowmelt from excess heat in rain it
            # takes 203.2 calories / (in cm^2) to melt snow (latent heat of fusion).
            # Convert from 1 cm depth over 1 square cm to
            # 1 inch depth over 1 square cm 80.0*(INCH2CM = 2.54 cm/in) = 203.2
            pmlt = dif / 203.2  # [inches]

            # Actual snowmelt can only come from snow covered area, so to calculate
            # the actual potential snowmelt, the potential snowmelt from snowcovered
            # area must be re-normalized to HRU area (rather than snowcover area).
            # In effect, the potential snowmelt per area is reduced by the fraction of
            # the watershed that is actually covered by snow.
            apmlt = pmlt * self.snowcov_area  # [inches]

            # Set the heat deficit and temperature of the remaining snowpack to 0.
            self.pk_def = zero  # [cal/cm^2]
            self.pk_temp = zero  # [degrees C]

            # The only pack ice that is melted is in the snow covered area, so the
            # pack ice needs to be re-normalized to the snowcovered area (rather than
            # HRU area). In effect, the pack ice per area is increased by the fraction
            # of the watershed that is actually covered by snow.
            if self.snowcov_area > zero:
                apk_ice = self.pk_ice / self.snowcov_area  # [inches]
            else:
                # print *, 'snowcov_area really small, melt all ice', snowcov_area, ' pmlt:', pmlt, ' dif:', dif, ' pk_ice:', pk_ice
                apk_ice = zero

            # If snow is melting, the heat is handled based on whether all or only
            # part of the pack ice melts.
            # 2 options below (if-then, else)

            if pmlt > apk_ice:
                # (3.1) Heat applied to snow covered area is sufficient to melt all the
                #       ice in that snow pack.
                # All pack water equivalent becomes meltwater.
                self.snowmelt = self.snowmelt + self.pkwater_equiv  # [inches]
                self.pkwater_equiv = zero  # [inches]
                # iasw = 0  # [flag]
                self.iasw = False  # [flag]

                # Set all snowpack states to 0.
                # snowcov_area = zero  # [fraction of area] # shouldn't be changed with melt
                self.pk_def = zero  # [cal / cm^2]
                self.pk_temp = zero  # [degreees C]
                self.pk_ice = zero  # [inches]
                self.freeh2o = zero  # [inches]
                self.pk_depth = zero  # [inches]
                self.pss = zero  # [inches]
                self.pst = zero  # [inches]
                self.pk_den = zero  # [fraction of depth]

            else:
                # (3.2) Heat only melts part of the ice in the snow pack.
                # Remove actual melt from frozen water and add melt to free water.
                self.pk_ice = self.pk_ice - apmlt  # [inches]
                self.freeh2o = self.freeh2o + apmlt  # [inches]

                # Calculate the capacity of the snowpack to hold free water according to
                # its current level of frozen water.
                pwcap = self.freeh2o_cap * self.pk_ice  # [inches]

                # Calculate the amount of free water in excess of the capacity to hold
                # free water.
                dif_water = freeh2o - pwcap  # [inches]

                # If there is more free water than the snowpack can hold, then there is
                # going to be melt...
                if dif_water > zero:
                    if dif_water > self.pkwater_equiv:
                        dif_water = self.pkwater_equiv

                    # total packwater decreases by the excess and a new depth is calculated
                    # based on density.
                    self.pkwater_equiv = (
                        self.pkwater_equiv - dif_water
                    )  # [inches]

                    # Free water is at the current capacity.
                    self.freeh2o = pwcap  # [inches]
                    if self.pk_den > zero:
                        self.pk_depth = (
                            self.pkwater_equiv / self.pk_den
                        )  # [inches]
                        # RAPCOMMENT - Added the conditional statement to make sure there is
                        #              no division by zero (this can happen if there is a
                        #              mixed event on no existing snowpack because a pack
                        #              density has not been calculated, yet).
                    else:
                        # rsr, this should not happen, remove later
                        # if print_debug > -1:
                        #    print *, 'snow density problem', pk_depth, pk_den, pss, pkwater_equiv
                        # call print_date(1)

                        self.pk_den = self.den_max
                        self.pk_depth = (
                            self.pkwater_equiv * self.denmaxinv
                        )  # [inches]

                    # Snowmelt increases by the excess free water.
                self.snowmelt = self.snowmelt + dif_water  # [inches]

                # Reset the previous-snowpack-plus-new-snow to the current pack water equivalent.
                self.pss = self.pkwater_equiv  # [inches]

        return
