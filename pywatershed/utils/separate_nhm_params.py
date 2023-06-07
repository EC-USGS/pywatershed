import numpy as np

import pywatershed

from ..constants import fileish
from ..parameters import PrmsParameters

# these are prarameters that are provided as scalars which we will
# force expand to their full dimensions
params_expand_scalar_to_dims = {
    "obsout_segment": "nsegment",
    "seg_humidity": "nsegment",
    "width_m": "nsegment",
}

var_meta_to_attrs = [
    "default",
    "desc",
    "help",
    "maximum",
    "minimum",
    "units",
]


def shape(val):
    if isinstance(val, np.ndarray):
        return val.shape
    elif isinstance(val, (int, float)):
        return 1
    else:
        return len(val)


def separate_domain_params_to_ncdf(
    prms_param_file: fileish,
    domain_name: str,
    out_dir: fileish,
    process_list: list = None,
    use_xr=True,
):
    """Separate PRMS parameter file into files for individual processes

    Args:
        prms_param_file: the native PRMS parameter file to separate
        domain_name: string name to include in the output files
        out_dir: the directory to which output files will be writen
        process_list: optional, the list of process classes desired for
            individual output files. If not specified, all process classes
            will be assumed.

    Returns:
        A dictionary of `process_class: file` pairs corresponding to the files
        written for (requested) processes in process_list.

    """
    nhm_processes = [
        pywatershed.PRMSSolarGeometry,
        pywatershed.PRMSAtmosphere,
        pywatershed.PRMSCanopy,
        pywatershed.PRMSSnow,
        pywatershed.PRMSRunoff,
        pywatershed.PRMSSoilzone,
        pywatershed.PRMSGroundwater,
        pywatershed.PRMSChannel,
    ]
    if process_list is None:
        process_list = nhm_processes

    prms_parameters = PrmsParameters.load(prms_param_file)

    written_files = {}

    for proc in process_list:
        proc_param_names = proc.get_parameters()
        proc_params = prms_parameters.subset(proc_param_names)
        nc_out_file = out_dir / f"parameters_{domain_name}_{proc.__name__}.nc"
        proc_params.to_netcdf(nc_out_file, use_xr=use_xr)
        written_files[proc] = nc_out_file

    return written_files
