from typing import Literal

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import nan, zero
from ..parameters import Parameters
from .prms_groundwater import PRMSGroundwater


class PRMSGroundwaterNoDprst(PRMSGroundwater):
    """PRMS groundwater reservoir.

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
        soil_to_gw: Portion of excess flow to the capillary reservoir that
            drains to the associated GWR for each HRU
        ssr_to_gw: Drainage from the gravity-reservoir to the associated GWR
            for each HRU
        dprst_seep_hru: Seepage from surface-depression storage to associated
            GWR for each HRU
        budget_type: one of ["defer", None, "warn", "error"] with "defer" being
            the default and defering to control.options["budget_type"] when
            available. When control.options["budget_type"] is not avaiable,
            budget_type is set to "warn".
        calc_method: one of ["fortran", "numba", "numpy"]. None defaults to
            "numba".
        verbose: Print extra information or not?

    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        soil_to_gw: adaptable,
        ssr_to_gw: adaptable,
        budget_type: Literal["defer", None, "warn", "error"] = "defer",
        calc_method: Literal["fortran", "numba", "numpy"] = None,
        verbose: bool = None,
    ) -> None:
        self._dprst_flag = False

        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
            soil_to_gw=soil_to_gw,
            ssr_to_gw=ssr_to_gw,
            dprst_seep_hru=None,
            budget_type=budget_type,
            calc_method=calc_method,
            verbose=verbose,
        )

        self.name = "PRMSGroundwaterNoDprst"
        self._set_budget()

        return

    @staticmethod
    def get_dimensions() -> tuple:
        return ("nhru",)

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "hru_area",
            "hru_in_to_cf",
            "gwflow_coef",
            "gwsink_coef",
            "gwstor_init",
            "gwstor_min",
        )

    @staticmethod
    def get_inputs() -> tuple:
        return (
            "soil_to_gw",
            "ssr_to_gw",
        )

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "soil_to_gw",
                "ssr_to_gw",
            ],
            "outputs": [
                "gwres_flow",
            ],
            "storage_changes": [
                "gwres_stor_change",
            ],
        }

    @staticmethod
    def get_init_values() -> dict:
        return {
            "gwres_flow": nan,
            "gwres_flow_vol": nan,
            "gwres_sink": nan,
            "gwres_stor": nan,
            "gwres_stor_old": nan,
            "gwres_stor_change": nan,
        }

    def _calculate(self, simulation_time):
        zero_array = self.gwres_stor * zero
        self._simulation_time = simulation_time
        (
            self.gwres_stor[:],
            self.gwres_flow[:],
            self.gwres_sink[:],
            self.gwres_stor_change[:],
            self.gwres_flow_vol[:],
        ) = self._calculate_gw(
            self.hru_area,
            self.soil_to_gw,
            self.ssr_to_gw,
            zero_array,
            self.gwres_stor,
            self.gwflow_coef,
            self.gwsink_coef,
            self.gwres_stor_old,
            self.hru_in_to_cf,
        )
        return
