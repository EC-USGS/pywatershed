import numpy as np

from ..base.data_model import DatasetDict
from ..constants import HruType, ACTIVE
from ..parameters import Parameters

import xarray as xr


# NOTES:
# * A preprocess needs to return new parametr objects or files for
#   PRMSRunoff, PRMSSoilzone, and MAYBE PRMSGroundwater
# * This combines functionality from cascade.f90 and basin.f90
# * Neglecting/commenting basin calculations below.
# * Is it strage that hru_route_order is calculated here?

# Is this a discretization edit?


def calc_hru_route_order(parameters: Parameters):
    """Calculate the HRU routing order.

    This is taken from basin.f90 with all its error trapping checks.

    Args:
      parameters: A Parameters object for the domain which includes hru_type

    """
    nhru = parameters.dims["nhru"]
    hru_type = parameters.parameters["hru_type"]
    hru_route_order = np.zeros(nhru, dtype=np.int32)

    nlake = 0
    if nlake > 0:
        numlake_hrus = 0  # to verify we have all the lakes
        numlakes_check = 0
        lake_hru_id = parameters.parameters["hru_lake_id"]

    frozen_flag = HruType.INACTIVE.value  # until further notice

    jj = -1
    for ii in range(nhru):

        if hru_type[ii] == HruType.INACTIVE.value:
            continue

        # ?? need to fix for lakes with multiple HRUs and PRMS lake routing ??
        if hru_type[ii] == HruType.LAKE.value:
            numlake_hrus = numlake_hrus + 1
            lakeid = lake_hru_id[ii]
            if lakeid > 0:
                if lakeid > numlakes_check:
                    numlakes_check = lakeid
            else:
                msg = f"ERROR, hru_type = 2 for HRU: {ii} and lake_hru_id = 0"
                raise ValueError(msg)

            if Nlake == 0:
                msg = (
                    f"ERROR, hru_type = 2 for HRU: {ii} "
                    "and dimension nlake = 0"
                )
                raise ValueError(msg)

        else:
            if nlake > 0 and lake_hru_id[ii] > 0:
                msg = (
                    f"ERROR, HRU: {ii} specifed to be a lake by lake_hru_id "
                    "but hru_type not equal 2"
                )
                raise ValueError(msg)

            if frozen_flag == ACTIVE:
                if hru_type[ii] == HruType.SWALE.value:
                    msg = (
                        "ERROR, a swale HRU cannot be frozen for CFGI, "
                        f"HRU: {ii}"
                    )
                    raise ValueError(msg)

        # <<<
        jj = jj + 1
        hru_route_order[jj] = ii

        if hru_type[ii] == HruType.LAKE.value:
            continue

    # <
    if nlake > 0:
        if numlakes_check != nlake:
            msg = (
                "ERROR, number of lakes specified in lake_hru_id"
                f"does not equal dimension nlake: {nlake} number of "
                "lakes: numlakes_check. For PRMS lake routing each lake "
                "must be a single HRU."
            )
            raise ValueError(msg)

    new_params = parameters.to_xr_ds()
    new_params["hru_route_order"] = xr.Variable("nhru", hru_route_order)
    return Parameters.from_dataset_dict(DatasetDict.from_ds(new_params))
