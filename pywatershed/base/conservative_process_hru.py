from typing import Literal

import numpy as np

from ..base import meta
from ..constants import HruType, zero
from ..parameters import Parameters
from .conservative_process import ConservativeProcess
from .control import Control


class ConservativeProcessHru(ConservativeProcess):
    """Base class for conservative processes on un-routed HRUs, gridded or not.

    A subclass of ConservativeProcess.
    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        budget_type: Literal["defer", None, "warn", "error"] = "defer",
        metadata_patches: dict[dict] = None,
        metadata_patch_conflicts: Literal["left", "warn", "error"] = "error",
    ):
        if not hasattr(self, "name"):
            self.name = "ConservativeProcess"

        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
            metadata_patches=metadata_patches,
            metadata_patch_conflicts=metadata_patch_conflicts,
        )

        self._set_active_locations()
        self._set_active_mask()
        self._active_mask_variables()

    def _set_active_locations(self):
        """Set the variables _active_locations and _nactive_locations."""
        if self.hru_type.min() == zero:
            self._active_hrus = np.where(self.hru_type != 0)[0]
            self._nactive_hrus = len(self._active_hrus)
        else:
            self._active_hrus = np.arange(self.nhru)
            self._nactive_hrus = self.nhru

    def _set_active_mask(self):
        """Set the _active_mask variable from HruType

        To be overridden for processes that maintain TimeseriesArrays.
        """
        if self.hru_type.min() == HruType.INACTIVE.value:
            self._active_mask = self.hru_type != HruType.INACTIVE.value

    def _active_mask_variables(self):
        """Set all variables to missing values outside of _active_mask."""
        if not hasattr(self, "_active_mask"):
            return
        # Implemented for 1-D variables not TimeseriesArrays
        for var_name in self.get_variables():
            var_dim_name = list(meta.get_dimensions(var_name).values())[0][0]
            if var_dim_name != "nhru":
                return
            self[var_name][:] = np.where(
                self._active_mask, self[var_name], np.nan
            )

        return
