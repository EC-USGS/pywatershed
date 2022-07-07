# This requires flopy to be installed
# which is not in the requirements.

import os
import pathlib as pl
import shutil

from flopy import run_model


def test_exe_available(exe):
    assert exe.is_file(), f"'{exe}'...does not exist"
    assert os.access(exe, os.X_OK)


def test_prms_run(simulations, exe):
    ws = simulations["ws"]
    control_file = simulations["control_file"]
    print(f"\n\n\n{'*' * 70}\n{'*' * 70}")
    print(
        f"run_domains.py: Running '{control_file}' in {ws}\n\n",
        flush=True,
    )

    # delete the existing output dir and re-create it
    output_dir = pl.Path(ws) / "output"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    success, buff = run_model(
        exe,
        control_file,
        model_ws=ws,
        cargs=[
            "-MAXDATALNLEN",
            "60000",
        ],
        normal_msg="Normal completion of PRMS",
    )
    assert success, f"could not run prms model in '{ws}'"

    print(f"run_domains.py: End of domain {ws}\n", flush=True)
