from typing import Union

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.parameters import PrmsParameters

from ..base.adapter import Adapter, adapter_factory
from ..base.control import Control
from ..constants import zero

adaptable = Union[str, np.ndarray, Adapter]

FT2_PER_ACRE = 43560.0
INCHES_PER_FOOT = 12.0
SECS_PER_HOUR = 3600.0


class PRMSChannel(StorageUnit):
    """PRMS channel flow

    Args:
        params: parameter object

    """

    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        sroff: adaptable,
        ssres_flow: adaptable,
        gwres_flow: adaptable,
        verbose: bool = False,
    ) -> "PRMSChannel":

        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
        )

        self.name = "PRMSChannel"

        self._input_variables_dict = {}
        self._input_variables_dict["sroff"] = adapter_factory(
            sroff,
            "sroff",
        )
        self._input_variables_dict["ssres_flow"] = adapter_factory(
            ssres_flow,
            "ssres_flow",
        )
        self._input_variables_dict["gwres_flow"] = adapter_factory(
            gwres_flow, "gwres_flow"
        )

        # process channel data
        self._process_channel_data()

        return

    def _process_channel_data(self):
        # convert prms data to zero-based
        self.hru_segment -= 1

        # calculate connectivity

    def set_initial_conditions(self):
        # initialize channel segment storage
        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get channel segment parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nhru",
            "nsegment",
            "hru_area",
            "hru_segment",
            "mann_n",
            "seg_depth",
            "seg_length",
            "seg_slope",
            "segment_type",
            "tosegment",
            "tosegment_nhm",
            "x_coef",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get channel segment input variables

        Returns:
            variables: input variables

        """
        return (
            "sroff",
            "ssres_flow",
            "gwres_flow",
        )

    @staticmethod
    def get_variables() -> tuple:
        """Get channel segment output variables

        Returns:
            variables: output variables

        """
        return ("seg_lateral_inflow",)

    @staticmethod
    def get_init_values() -> dict:
        """Get channel segment initial values

        Returns:
            dict: initial values for named variables
        """
        # No GW res values need initialized prior to calculation.
        return {
            "seg_lateral_inflow": zero,
        }

    def _advance_variables(self) -> None:
        """Advance the channel segment variables
        Returns:
            None
        """
        return

    def calculate(self, simulation_time):
        """Calculate channel segment terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None

        """

        self._simulation_time = simulation_time

        hru_area = self.hru_area
        # FT2_PER_ACRE / INCHES_PER_FOOT / Timestep_seconds
        time_step_seconds = self.control.time_step / np.timedelta64(1, "s")
        in_to_cfs = hru_area * FT2_PER_ACRE
        in_to_cfs /= INCHES_PER_FOOT
        in_to_cfs /= time_step_seconds

        # calculate lateral flow term
        self.seg_lateral_inflow[:] = 0.0
        for ihru in range(self.nhru):
            irch = self.hru_segment[ihru]
            if irch < 0:
                continue
            lateral_inflow = (
                self.sroff[ihru]
                + self.ssres_flow[ihru]
                + self.gwres_flow[ihru]
            ) * in_to_cfs[ihru]
            self.seg_lateral_inflow[irch] += lateral_inflow

        return
