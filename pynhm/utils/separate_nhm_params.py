import numpy as np

import pynhm

from ..base import meta
from ..constants import fileish
from ..utils import PrmsParameters

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
    domain_name: str,
    prms_param_file: fileish,
    out_dir: fileish,
    process_list=None,
):
    import xarray as xr

    nhm_processes = [
        pynhm.PRMSSolarGeometry,
        pynhm.PRMSAtmosphere,
        pynhm.PRMSCanopy,
        pynhm.PRMSSnow,
        pynhm.PRMSRunoff,
        pynhm.PRMSSoilzone,
        pynhm.PRMSGroundwater,
        pynhm.PRMSChannel,
    ]
    if process_list is None:
        process_list = nhm_processes

    param_dict = PrmsParameters.load(prms_param_file).parameters

    written_files = {}

    # Get the dimensions up front, use as needed later.
    meta_dim_names = list(meta.dimensions.keys())
    dim_meta = {}
    for dim in meta_dim_names:
        if dim in param_dict.keys():
            dim_meta[dim] = int(param_dict[dim])

    for proc in process_list:
        proc_param_names = proc.get_parameters()
        proc_params = {name: param_dict[name] for name in proc_param_names}

        # need all dims for process over inputs, variables & parameters
        proc_vars = proc.get_variables()
        proc_inputs = proc.get_inputs()
        proc_dims = set()
        for pp in [proc_params, proc_vars, proc_inputs]:
            for name in pp:
                mm = meta.find_variables(name)[name]
                if "dimensions" in mm.keys():
                    proc_dims = proc_dims.union(mm["dimensions"].values())
        proc_dims = list(proc_dims)

        ds = xr.Dataset(
            attrs=dict(
                description=(
                    f"NHM Parameters for {proc.__name__} in "
                    f"{domain_name} domain"
                ),
                domain_name=domain_name,
                nhm_process=proc.__name__,
            ),
        )
        nc_out_file = out_dir / f"parameters_{domain_name}_{proc.__name__}.nc"

        # dimensions
        dim_dict = {}
        for dim in meta_dim_names:
            if dim in proc_params.keys():
                dim_dict[dim] = proc_params[dim]
        # these are implied and not included
        implied_dims = ["ndoy", "nmonth", "scalar"]
        for dd in implied_dims:
            dim_dict[dd] = meta.dimensions[dd]["default"]

        # Add parameters with dimensions
        for param_name in proc_params.keys():
            if param_name in dim_dict.keys():
                continue

            param_meta = meta.find_variables(param_name)[param_name]
            param_vals = proc_params[param_name]

            param_shape = shape(param_vals)
            if isinstance(param_shape, int):
                param_shape = [param_shape]

            param_dims = {
                dim_name: dim_dict[dim_name]
                for dim_name in param_meta["dimensions"].values()
            }
            for pp in param_dims:
                # remove dims used for parameter data
                if pp in proc_dims:
                    proc_dims.remove(pp)

            # sometimes a scalar is allowed to represent a uniform values for
            # the full dimensions of a parameter. Going to handle those on a
            # case-by-case basis for now: just expand them to full size
            if param_name in params_expand_scalar_to_dims:
                for ii, dim_name in enumerate(param_dims.keys()):
                    full_dim = params_expand_scalar_to_dims[param_name]
                    if dim_name == full_dim:
                        param_shape = list(param_shape)
                        param_shape[ii] = dim_dict[full_dim]
                        param_shape = tuple(param_shape)
                        param_vals = np.zeros(param_shape) + param_vals

            assert tuple(param_dims.values()) == param_shape

            ds[param_name] = xr.DataArray(data=param_vals, dims=param_dims)
            ds[param_name].attrs = {
                kk: vv
                for kk, vv in param_meta.items()
                if kk in var_meta_to_attrs
            }

        # to file
        for vv in proc_dims:
            ds.attrs[vv] = dim_meta[vv]
        ds.to_netcdf(nc_out_file)
        ds.close()
        assert nc_out_file.exists()
        written_files[proc] = nc_out_file

    return written_files
