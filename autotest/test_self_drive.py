import pathlib as pl

import pytest
import pywatershed as pws
import xarray as xr

n_time_steps = 50


@pytest.fixture(scope="function")
def process_list(simulation):
    control = pws.Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )

    if (
        "dprst_flag" in control.options.keys()
        and control.options["dprst_flag"]
    ):
        Runoff = pws.PRMSRunoff
        Soilzone = pws.PRMSSoilzone
        Groundwater = pws.PRMSGroundwater

    else:
        Runoff = pws.PRMSRunoffNoDprst
        Soilzone = pws.PRMSSoilzoneNoDprst
        Groundwater = pws.PRMSGroundwaterNoDprst

    process_list = [
        pws.PRMSSolarGeometry,
        pws.PRMSAtmosphere,
        pws.PRMSCanopy,
        pws.PRMSSnow,
        Runoff,
        Soilzone,
        Groundwater,
    ]

    if control.options["streamflow_module"] != "strmflow":
        process_list += [pws.PRMSChannel]
    return process_list


# TODO: instead of individual processes, could test the last N processes
#       where N in [2,6]


@pytest.mark.domain
def test_drive_indiv_process(simulation, process_list, tmp_path):
    """Output of a full pywatershed drives indiv process models separately

    The results from the full model should be consistent with the results from
    the individual models, else there is likely something wrong with the
    full model.
    """

    # Run a full pws model to use its output to drive individual processes
    output_dir = pl.Path(tmp_path) / "output"

    control = pws.Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )

    control.edit_n_time_steps(n_time_steps)
    control.options["budget_type"] = "warn"
    control.options["calc_method"] = "numba"
    control.options["input_dir"] = simulation["dir"]
    del control.options["netcdf_output_var_names"]
    del control.options["netcdf_output_dir"]

    param_file = simulation["dir"] / control.options["parameter_file"]
    params = pws.parameters.PrmsParameters.load(param_file)

    model = pws.Model(
        process_list,
        control=control,
        parameters=params,
    )

    model.initialize_netcdf(output_dir=output_dir)

    model.run(finalize=True)
    del model, params, control

    # run individual process models
    for proc in process_list:
        if proc in [pws.PRMSSolarGeometry, pws.PRMSAtmosphere]:
            # These are not driven by outputs of above, only external outputs
            # or known/static inputs
            continue

        print(proc.__name__)
        proc_model_output_dir = tmp_path / proc.__name__
        proc_model_output_dir.mkdir()

        control = pws.Control.load_prms(
            simulation["control_file"], warn_unused_options=False
        )
        param_file = simulation["dir"] / control.options["parameter_file"]
        params = pws.parameters.PrmsParameters.load(param_file)

        control.edit_n_time_steps(n_time_steps)
        control.options["budget_type"] = "warn"
        control.options["calc_method"] = "numba"
        control.options["input_dir"] = output_dir
        control.options["netcdf_output_dir"] = proc_model_output_dir
        del control.options["netcdf_output_var_names"]

        proc_model = pws.Model(
            [proc],
            control=control,
            parameters=params,
        )
        proc_model.initialize_netcdf()
        proc_model.run(finalize=True)
        del proc_model, params, control

        # compare output netcdf files
        for vv in proc.get_variables():
            results_file = proc_model_output_dir / f"{vv}.nc"
            if not results_file.exists():
                print(f"results file not found: {results_file}")
                continue
            res = xr.open_dataset(proc_model_output_dir / f"{vv}.nc")
            ans = xr.open_dataset(output_dir / f"{vv}.nc")

            # Leaving the commented to diagnose what PRMSRunoff later.
            try:
                xr.testing.assert_allclose(res, ans)
            except AssertionError:
                print(vv, abs(res - ans).max())
                print(vv, (abs(res - ans) / ans).max())

            del res, ans

    return
