import numpy as np

import pywatershed

from ..base import meta
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
    domain_name: str,
    prms_param_file: fileish,
    out_dir: fileish,
    process_list=None,
):
    import xarray as xr

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
    param_dict = prms_parameters.parameters

    # TODO: PRMSParameters should merge the static metadata when this is done,
    #       it will very much transform the following code

    written_files = {}

    for proc in process_list:
        proc_param_names = proc.get_parameters()
        proc_params = {name: param_dict[name] for name in proc_param_names}

        # in some cases (e.g. PRMSAtmosphere) the parameters dont need
        # all the process dimensions
        # proc_dims = proc.get_dimensions()  ## too much
        proc_dims = set(
            [
                vv["dims"]
                for kk, vv in prms_parameters.metadata.items()
                if kk in proc_params.keys()
            ]
        )
        proc_dims = tuple([dd for tt in proc_dims for dd in tt])

        ds = xr.Dataset(
            attrs=dict(
                description=(
                    f"NHM Parameters for {proc.__name__} in " f"{domain_name} domain"
                ),
                domain_name=domain_name,
                nhm_process=proc.__name__,
            ),
        )
        nc_out_file = out_dir / f"parameters_{domain_name}_{proc.__name__}.nc"

        # dimension
        # Add parameters with dimensions
        for param_name in proc_params.keys():
            param_meta = meta.find_variables(param_name)[param_name]

            param_vals = proc_params[param_name]

            param_shape = shape(param_vals)
            if isinstance(param_shape, int):
                param_shape = [param_shape]

            param_dim_names = prms_parameters.metadata[param_name]["dims"]
            param_dims = {
                dim_name: prms_parameters.dims[dim_name] for dim_name in param_dim_names
            }

            proc_dims = list(proc_dims)
            for pp in param_dims:
                # remove dims used for parameter data
                if pp in proc_dims:
                    proc_dims.remove(pp)

            proc_dims = tuple(proc_dims)

            # sometimes a scalar is allowed to represent a uniform values for
            # the full dimensions of a parameter. Going to handle those on a
            # case-by-case basis for now: just expand them to full size
            if param_name in params_expand_scalar_to_dims:
                for ii, dim_name in enumerate(param_dims.keys()):
                    full_dim = params_expand_scalar_to_dims[param_name]
                    if dim_name == full_dim:
                        param_shape = list(param_shape)
                        param_shape[ii] = prms_parameters.dims[full_dim]
                        param_shape = tuple(param_shape)
                        param_vals = np.zeros(param_shape) + param_vals

            assert tuple(param_dims.values()) == param_shape

            ds[param_name] = xr.DataArray(data=param_vals, dims=param_dims)
            ds[param_name].attrs = {
                kk: vv for kk, vv in param_meta.items() if kk in var_meta_to_attrs
            }

        # to file
        for vv in proc_dims:
            ds.attrs[vv] = prms_parameters.dims[vv]
        ds.to_netcdf(nc_out_file)
        ds.close()
        assert nc_out_file.exists()
        written_files[proc] = nc_out_file

    return written_files
