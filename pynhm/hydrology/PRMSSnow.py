from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control
from ..constants import HruType, one, zero

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
        """Get snow pack initial values

        Returns:
            dict: initial values for named variables
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
        self.deninv = one / self.den_init.copy()
        self.denmaxinv = one / self.den_max.copy()

        self.pkwater_equiv = self.snowpack_init.copy()
        self.pkwater_equiv_ante = None

        sd = int(self.ndeplval / 11)
        self.snarea_curve_2d = np.reshape(self.snarea_curve, (sd, 11))

        if self.control.config["init_vars_from_file"] in [0, 2, 3]:

            # The super().__init__ already set_initial_conditions using its set_initial_conditions
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

    def _advance_variables(self) -> None:
        self.pkwater_ante = self.pkwater_equiv
        return

    def _advance_inputs(self) -> None:
        """Advance the snow pack inputs
        Returns:
            None
        """
        # JLM: This method is only because adapter advances dont all take current time.
        for key, value in self._input_variables_dict.items():
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
        newsnow = np.zeros([self.nhru], dtype=int)
        pptmix = np.zeros([self.nhru])

        net_snow_gt_zero = self.net_snow > zero
        wh_net_snow_gt_zero = np.where(net_snow_gt_zero)
        newsnow[wh_net_snow_gt_zero] = 1

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

        # HRU STEP 1 - DEAL WITH PRECIPITATION AND ITS EFFECT ON THE WATER
        #              CONTENT AND HEAT CONTENT OF SNOW PACK
        # ***********************************************************************

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
