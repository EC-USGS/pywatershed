from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control
from ..constants import one, zero

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
        # soltab_horad_potsw: adaptable,
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
    ) -> "PRMSGroundwater":

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
        """Get groundwater reservoir parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nhru",
            "ndeplval",
            "nmonths",
            # "cov_type",
            # "hru_type",
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
        """Get groundwater reservoir input variables

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
            # "soltab_horad_potsw",
            "swrad",
            "tavgc",
            "tmaxc",
            "tminc",
            "transp_on",
        )

    @staticmethod
    def get_variables() -> tuple:
        """Get groundwater reservoir output variables

        Returns:
            variables: output variables

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

    def set_initial_conditions(self):
        """Initialize PRMSSnow snowpack variables."""

        # I have set metadata on each var for var_init: [zero, False] to init those vars

        # (should these be labeled as such in metadata?)
        self.deninv = one / self.den_init.copy()
        self.denmaxinv = one / self.den_max.copy()
        self.pkwater_equiv = self.snowpack_init.copy()
        self.pkwater_equiv_ante = None

        sd = self.ndeplval / 11
        # self.snarea_curve_2d = reshape(self.snarea_curve, (/11, sd/))

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
        self.pk_ice[wh_pkweq_gt_zero] = self.pkwater_equiv[wh_pkweq_gt_zero]
        self.freeh2o[wh_pkweq_gt_zero] = (
            self.pk_ice[wh_pkweq_gt_zero] * self.freeh2o_cap[wh_pkweq_gt_zero]
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

        # self.snowcov_area[wh_pkweq_gt_zero] = self.sca_deplcrv(self.snarea_curve_2d(1:11, self.hru_deplcrv(chru)), &
        #                                           self.frac_swe[wh_pkweq_gt_zero])

        self.pss = self.pkwater_equiv
        self.pst = self.pkwater_equiv

        # now get from restart
        # if restart then get these variables from the restart file.
        # JLM: this seems like an overkill list and dosent shed light on what states have memory.
        # JLM: could there be a diagnostic part of the advance?
        restart_variables = [
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
        ]
        return

    def advance(self) -> None:
        """Advance the snow pack
        Returns:
            None

        """
        self.pkwater_ante = self.pkwater_equiv
        self._itime_step += 1

        for key, value in self._input_variables_dict.items():
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

        newsnow = np.zeros([self.nhru], dtype=int)
        net_snow_gt_zero = self.net_snow > zero
        wh_net_snow_gt_zero = np.where(net_snow_gt_zero)
        newsnow[wh_net_snow_gt_zero] = 1

        net_rain_gt_zero = self.net_rain > zero
        wh_net_snow_gt_zero_and_net_rain_gt_zero = np.where(
            net_rain_gt_zero & net_rain_gt_zero
        )
        pptmix = np.zeros([self.nhru])

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
