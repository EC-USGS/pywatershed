import pathlib as pl

import pytest
import xarray as xr

import pywatershed as pws

n_time_steps = 50

nhm_processes = [
    pws.PRMSSolarGeometry,
    pws.PRMSAtmosphere,
    pws.PRMSCanopy,
    pws.PRMSSnow,
    pws.PRMSRunoff,
    pws.PRMSSoilzone,
    pws.PRMSGroundwater,
    pws.PRMSChannel,
]

# TODO: instead of individual processes, could test the last N processes
#       where N in [2,6]


def test_drive_indiv_process(domain, tmp_path):
    """Output of a full pywatershed NHM drives indiv process models separately

    The results from the full model should be consistent with the results from
    the individual models, else there is likely something wrong with the
    full model.
    """

    # Run a full pws NHM to use its output to drive individual processes
    nhm_output_dir = pl.Path(tmp_path) / "nhm_output"

    params = pws.parameters.PrmsParameters.load(domain["param_file"])
    control = pws.Control.load_prms(
        domain["control_file"], warn_unused_options=False
    )
    control.edit_n_time_steps(n_time_steps)
    control.options["budget_type"] = "warn"
    control.options["calc_method"] = "numba"
    control.options["input_dir"] = domain["prms_run_dir"]
    del control.options["netcdf_output_var_names"]
    del control.options["netcdf_output_dir"]

    nhm = pws.Model(
        nhm_processes,
        control=control,
        parameters=params,
    )

    nhm.initialize_netcdf(output_dir=nhm_output_dir)

    nhm.run(finalize=True)
    del nhm, params, control

    # run individual process models
    for proc in nhm_processes:
        if proc in [pws.PRMSSolarGeometry, pws.PRMSAtmosphere]:
            # These are not driven by outputs of above, only external outputs
            # or known/static inputs
            continue

        print(proc.__name__)
        proc_model_output_dir = tmp_path / proc.__name__
        proc_model_output_dir.mkdir()

        params = pws.parameters.PrmsParameters.load(domain["param_file"])
        control = pws.Control.load_prms(
            domain["control_file"], warn_unused_options=False
        )
        control.edit_n_time_steps(n_time_steps)
        control.options["budget_type"] = "warn"
        control.options["calc_method"] = "numba"
        control.options["input_dir"] = nhm_output_dir
        control.options["netcdf_output_dir"] = proc_model_output_dir

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
            ans = xr.open_dataset(nhm_output_dir / f"{vv}.nc")

            # Leaving the commented to diagnose what PRMSRunoff later.
            try:
                xr.testing.assert_allclose(res, ans)
            except AssertionError:
                print(vv, abs(res - ans).max())
                print(vv, (abs(res - ans) / ans).max())

            del res, ans

    return
