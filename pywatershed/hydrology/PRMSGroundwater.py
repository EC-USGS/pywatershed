import numpy as np

from pywatershed.base.storageUnit import StorageUnit

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import nan, numba_num_threads

try:
    from ..PRMSGroundwater_f import calc_groundwater as _calculate_fortran

    has_prmsgroundwater_f = True
except ImportError:
    has_prmsgroundwater_f = False


class PRMSGroundwater(StorageUnit):
    """PRMS groundwater reservoir."""

    def __init__(
        self,
        control: Control,
        soil_to_gw: adaptable,
        ssr_to_gw: adaptable,
        dprst_seep_hru: adaptable,
        budget_type: str = None,
        calc_method: str = None,
        verbose: bool = False,
        load_n_time_batches: int = 1,
    ):
        super().__init__(
            control=control,
            verbose=verbose,
            load_n_time_batches=load_n_time_batches,
        )
        self.name = "PRMSGroundwater"

        self._calc_method = str(calc_method)

        self._set_inputs(locals())
        self._set_budget(budget_type)

        if calc_method == "numba":
            # read-only arrays dont have numba signatures
            self._hru_area = self.hru_area.copy()
            self._gwflow_coef = self.gwflow_coef.copy()
            self._gwsink_coef = self.gwsink_coef.copy()

        return

    @staticmethod
    def get_dimensions() -> tuple:
        """Get groundwater reservoir dimensions

        Returns:
            tuple of dimension names

        """
        return ("nhru",)

    @staticmethod
    def get_parameters() -> tuple:
        """Get groundwater reservoir parameters

        Returns:
            tuple of parameter name

        """
        return (
            "hru_area",
            "gwflow_coef",
            "gwsink_coef",
            "gwstor_init",
            "gwstor_min",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get groundwater reservoir input variables

        Returns:
            tuple: input variable name

        """
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
        """Get groundwater initial values

        Returns:
            dict: initial values for named variables
        """
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

    def _advance_variables(self) -> None:
        """Advance the groundwater reservoir variables
        Returns:
            None
        """
        self.gwres_stor_old[:] = self.gwres_stor
        return

    def _calculate(self, simulation_time):
        """Calculate groundwater reservoir terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

        if self._calc_method.lower() == "numba":
            import numba as nb

            if not hasattr(self, "_calculate_numba"):
                numba_msg = f"{self.name} jit compiling with numba "
                nb_parallel = (numba_num_threads is not None) and (
                    numba_num_threads > 1
                )
                if nb_parallel:
                    numba_msg += f"and using {numba_num_threads} threads"
                print(numba_msg, flush=True)

                self._calculate_numba = nb.njit(
                    nb.types.UniTuple(nb.float64[:], 5)(
                        nb.float64[:],
                        nb.float64[:],
                        nb.float64[:],
                        nb.float64[:],
                        nb.float64[:],
                        nb.float64[:],
                        nb.float64[:],
                        nb.float64[:],
                        nb.float64[:],
                    ),
                    parallel=False,
                )(self._calculate_numpy)

            (
                self.gwres_stor[:],
                self.gwres_flow[:],
                self.gwres_sink[:],
                self.gwres_stor_change[:],
                self.gwres_flow_vol[:],
            ) = self._calculate_numba(
                self._hru_area,
                self.soil_to_gw,
                self.ssr_to_gw,
                self.dprst_seep_hru,
                self.gwres_stor,
                self._gwflow_coef,
                self._gwsink_coef,
                self.gwres_stor_old,
                self.params.hru_in_to_cf,
            )

        elif self._calc_method.lower() == "fortran":
            (
                self.gwres_stor[:],
                self.gwres_flow[:],
                self.gwres_sink[:],
                self.gwres_stor_change[:],
                self.gwres_flow_vol[:],
            ) = _calculate_fortran(
                self.hru_area,
                self.soil_to_gw,
                self.ssr_to_gw,
                self.dprst_seep_hru,
                self.gwres_stor,
                self.gwflow_coef,
                self.gwsink_coef,
                self.gwres_stor_old,
                self.params.hru_in_to_cf,
            )

        elif self._calc_method.lower() in ["none", "numpy"]:
            (
                self.gwres_stor[:],
                self.gwres_flow[:],
                self.gwres_sink[:],
                self.gwres_stor_change[:],
                self.gwres_flow_vol[:],
            ) = self._calculate_numpy(
                self.hru_area,
                self.soil_to_gw,
                self.ssr_to_gw,
                self.dprst_seep_hru,
                self.gwres_stor,
                self.gwflow_coef,
                self.gwsink_coef,
                self.gwres_stor_old,
                self.params.hru_in_to_cf,
            )

        else:
            msg = f"Invalid calc_method={self._calc_method} for {self.name}"
            raise ValueError(msg)

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
