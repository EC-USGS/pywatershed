from typing import Literal
from warnings import warn

import numpy as np

from ..base.adapter import adaptable
from ..base.conservative_process import ConservativeProcess
from ..base.control import Control
from ..constants import nan, numba_num_threads
from ..parameters import Parameters

try:
    from ..prms_groundwater_f import calc_groundwater as _calculate_fortran

    has_prmsgroundwater_f = True
except ImportError:
    has_prmsgroundwater_f = False


class PRMSGroundwater(ConservativeProcess):
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
        budget_type: one of [None, "warn", "error"]
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
        dprst_seep_hru: adaptable,
        budget_type: Literal[None, "warn", "error"] = None,
        calc_method: Literal["fortran", "numba", "numpy"] = None,
        verbose: bool = None,
    ) -> None:
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "PRMSGroundwater"

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
            "dprst_seep_hru",
        )

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "soil_to_gw",
                "ssr_to_gw",
                "dprst_seep_hru",
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

    def _set_initial_conditions(self):
        # initialize groundwater reservoir storage
        self.gwres_stor[:] = self.gwstor_init.copy()
        self.gwres_stor_old[:] = self.gwstor_init.copy()
        return

    def _init_calc_method(self):
        if self._calc_method is None:
            self._calc_method = "numba"

        avail_methods = ["numpy", "numba", "fortran"]
        fortran_msg = ""
        if self._calc_method == "fortran" and not has_prmsgroundwater_f:
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

        if self._calc_method.lower() == "numba":
            import numba as nb

            numba_msg = f"{self.name} jit compiling with numba "
            nb_parallel = (numba_num_threads is not None) and (
                numba_num_threads > 1
            )
            if nb_parallel:
                numba_msg += f"and using {numba_num_threads} threads"
            print(numba_msg, flush=True)

            self._calculate_gw = nb.njit(
                nb.types.UniTuple(nb.float64[:], 5)(
                    nb.types.Array(nb.types.float64, 1, "C", readonly=True),
                    nb.float64[:],
                    nb.float64[:],
                    nb.float64[:],
                    nb.float64[:],
                    nb.types.Array(nb.types.float64, 1, "C", readonly=True),
                    nb.types.Array(nb.types.float64, 1, "C", readonly=True),
                    nb.float64[:],
                    nb.types.Array(nb.types.float64, 1, "C", readonly=True),
                ),
                fastmath=True,
                parallel=False,
            )(self._calculate_numpy)

        elif self._calc_method.lower() == "fortran":
            self._calculate_gw = _calculate_fortran

        else:
            self._calculate_gw = self._calculate_numpy

        return

    def _advance_variables(self) -> None:
        self.gwres_stor_old[:] = self.gwres_stor
        return

    def _calculate(self, simulation_time):
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
            self.dprst_seep_hru,
            self.gwres_stor,
            self.gwflow_coef,
            self.gwsink_coef,
            self.gwres_stor_old,
            self.hru_in_to_cf,
        )
        return

    @staticmethod
    def _calculate_numpy(
        gwarea,
        soil_to_gw,
        ssr_to_gw,
        dprst_seep_hru,
        gwres_stor,
        gwflow_coef,
        gwsink_coef,
        gwres_stor_old,
        hru_in_to_cf,
    ):
        soil_to_gw_vol = soil_to_gw * gwarea
        ssr_to_gw_vol = ssr_to_gw * gwarea
        dprst_seep_hru_vol = dprst_seep_hru * gwarea

        # todo: what about route order

        _gwres_stor = gwres_stor * gwarea
        _gwres_stor += soil_to_gw_vol + ssr_to_gw_vol + dprst_seep_hru_vol

        _gwres_flow = _gwres_stor * gwflow_coef
        _gwres_stor -= _gwres_flow

        _gwres_sink = _gwres_stor * gwsink_coef
        idx = np.where(_gwres_sink > _gwres_stor)
        _gwres_sink[idx] = _gwres_stor[idx]
        _gwres_stor -= _gwres_sink

        # convert most units back to self variables
        # output variables
        gwres_stor = _gwres_stor / gwarea
        # for some stupid reason this is left in acre-inches
        gwres_flow = _gwres_flow / gwarea
        gwres_sink = _gwres_sink / gwarea

        gwres_stor_change = gwres_stor - gwres_stor_old
        gwres_flow_vol = gwres_flow * hru_in_to_cf

        return (
            gwres_stor,
            gwres_flow,
            gwres_sink,
            gwres_stor_change,
            gwres_flow_vol,
        )
