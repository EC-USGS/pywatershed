import pathlib as pl
import shutil
import time

import numpy as np
import pywatershed as pws
import xarray as xr

sttime = time.time()

model_output_netcdf = True

work_dir = pl.Path("/Users/jmccreight/usgs/pywatershed2/test_data/drb_2yr")

out_dir = pl.Path("./custom_output")
shutil.rmtree(out_dir)  # CAREFUL HERE
out_dir.mkdir()
custom_output_file = out_dir / "model_custom_output.nc"

param_file = work_dir / "myparam.param"
params = pws.parameters.PrmsParameters.load(param_file)

control = pws.Control.load(work_dir / "control.test")

control.options = control.options | {
    "input_dir": work_dir,
    "budget_type": None,
    "verbose": False,
    "calc_method": "numba",
}

if model_output_netcdf:
    control.options = control.options | {
        "netcdf_output_var_names": [
            "hru_actet",
            "sroff_vol",
            "ssres_flow_vol",
            "gwres_flow_vol",
            "seg_outflow",
            "hru_streamflow_out",
        ],
        "netcdf_output_dir": out_dir,
    }


model = pws.Model(
    [
        pws.PRMSSolarGeometry,
        pws.PRMSAtmosphere,
        pws.PRMSCanopy,
        pws.PRMSSnow,
        pws.PRMSRunoff,
        pws.PRMSSoilzone,
        pws.PRMSGroundwater,
        pws.PRMSChannel,
    ],
    control=control,
    parameters=params,
)


# Custom model output at selected spatial locations for all times.
# Generally, i'd be careful with xarray performance, but just writing at the
# end should be fine.
# Could move to netcdf4 if performance is a concern.

# /////////////////////////////////
# specfications: what we want this to look like to the user

var_list = [
    "hru_actet",
    "seg_outflow",
]

# want seg_outflow just on poi_gages
# make it a tuple like the return of np.where
wh_gages = (params.parameters["poi_gage_segment"] - 1,)
spatial_subsets = {
    "poi_gages": {
        "coord_name": "nhm_seg",
        "indices": wh_gages,
        "new_coord": params.parameters["poi_gage_id"],
        "variables": ["seg_outflow"],
    },
}


# A novel, diagnostic variable
def sum_hru_flows(sroff_vol, ssres_flow_vol, gwres_flow_vol):
    return sroff_vol + ssres_flow_vol + gwres_flow_vol


diagnostic_var_dict = {
    "hru_streamflow_out": {
        "inputs": ["sroff_vol", "ssres_flow_vol", "gwres_flow_vol"],
        "function": sum_hru_flows,
        "like_var": "sroff_vol",
        "metadata": {"desc": "something or other", "units": "parsecs"},
    },
}

# TODO: specify subsets in time
# TODO: specify different output files

# /////////////////////////////////
# code starts here

out_subset_ds = xr.Dataset()

needed_vars = var_list + [
    var for key, val in diagnostic_var_dict.items() for var in val["inputs"]
]
needed_metadata = pws.meta.get_vars(needed_vars)
dims = set([dim for val in needed_metadata.values() for dim in val["dims"]])

subset_vars = [
    var for key, val in spatial_subsets.items() for var in val["variables"]
]

var_subset_key = {
    var: subkey
    for var in subset_vars
    for subkey in spatial_subsets.keys()
    if var in spatial_subsets[subkey]["variables"]
}

diagnostic_vars = list(diagnostic_var_dict.keys())

# solve the processes for each variable
var_proc = {
    var: proc_key
    for var in needed_vars
    for proc_key, proc_val in model.processes.items()
    if var in proc_val.get_variables()
}

time_coord = np.arange(
    control.start_time, control.end_time, dtype="datetime64[D]"
)
n_time_steps = len(time_coord)
out_subset_ds["time"] = xr.Variable(["time"], time_coord)
out_subset_ds = out_subset_ds.set_coords("time")

# annoying to have to hard-code this
dim_coord = {"nhru": "nhm_id", "nsegment": "nhm_seg"}


# declare memory for the outputs
for var in var_list + diagnostic_vars:
    # impostor approach
    orig_diag_var = None
    if var in diagnostic_vars:
        orig_diag_var = var
        var = diagnostic_var_dict[var]["like_var"]

    proc = model.processes[var_proc[var]]
    dim_name = needed_metadata[var]["dims"][0]
    dim_len = proc.params.dims[dim_name]
    coord_name = dim_coord[dim_name]
    coord_data = proc.params.coords[dim_coord[dim_name]]
    type = needed_metadata[var]["type"]

    var_meta = {
        kk: vv
        for kk, vv in needed_metadata[var].items()
        if kk in ["desc", "units"]
    }

    if orig_diag_var is not None:
        var = orig_diag_var
        del var_meta["desc"]
        if "metadata" in diagnostic_var_dict[var]:
            var_meta = diagnostic_var_dict[var]["metadata"]
        if "desc" not in var_meta.keys():
            var_meta["desc"] = "Custom output diagnostic variable"

    if var in subset_vars:
        subset_key = var_subset_key[var]
        subset_info = spatial_subsets[subset_key]
        dim_name = f"n{subset_key}"
        coord_name = subset_key
        dim_len = len(subset_info["indices"][0])
        coord_data = subset_info["new_coord"]

    if coord_name not in list(out_subset_ds.variables):
        out_subset_ds[coord_name] = xr.DataArray(coord_data, dims=[dim_name])
        out_subset_ds = out_subset_ds.set_coords(coord_name)

    out_subset_ds[var] = xr.Variable(
        ["time", dim_name],
        np.full(
            [n_time_steps, dim_len],
            pws.constants.fill_values_dict[np.dtype(type)],
            type,
        ),
    )

    out_subset_ds[var].attrs = var_meta


for istep in range(n_time_steps):
    model.advance()
    model.calculate()

    if model_output_netcdf:
        model.output()

    for var in var_list:
        proc = model.processes[var_proc[var]]
        if var not in subset_vars:
            out_subset_ds[var][istep, :] = proc[var]
        else:
            indices = spatial_subsets[var_subset_key[var]]["indices"]
            out_subset_ds[var][istep, :] = proc[var][indices]

    for diag_key, diag_val in diagnostic_var_dict.items():
        input_dict = {}
        for ii in diag_val["inputs"]:
            proc = model.processes[var_proc[ii]]
            input_dict[ii] = proc[ii]

        out_subset_ds[diag_key][istep, :] = diag_val["function"](**input_dict)


out_subset_ds.to_netcdf(custom_output_file)

del model
del out_subset_ds

if model_output_netcdf:
    out_subset_ds = xr.open_dataset(custom_output_file)

    for vv in var_list:
        default_output_file = out_dir / f"{vv}.nc"
        print("checking variable: ", vv)
        answer = xr.open_dataset(default_output_file)[vv]
        result = out_subset_ds[vv]

        if vv in subset_vars:
            indices = spatial_subsets[var_subset_key[vv]]["indices"]
            answer = answer[:, indices[0]]

        np.testing.assert_allclose(answer, result)

    for diag_key, diag_val in diagnostic_var_dict.items():
        print("checking diagnostic variable: ", diag_key)
        input_dict = {}
        for ii in diag_val["inputs"]:
            default_output_file = out_dir / f"{ii}.nc"
            input_dict[ii] = xr.open_dataset(default_output_file)[ii]

        answer = diag_val["function"](**input_dict)
        result = out_subset_ds[diag_key]

        np.testing.assert_allclose(answer, result)
