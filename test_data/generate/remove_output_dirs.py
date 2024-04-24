import pathlib as pl
import shutil


def test_remove_output_dirs(simulation):
    output_dir = pl.Path(simulation["output_dir"])
    if output_dir.exists():
        shutil.rmtree(output_dir)

    # this shouldnt matter but it is necessary to exist for PRMS to write to it
    output_dir.mkdir(parents=True)

    # remove the corresponding test_data_version file
    # this has to be synched with autotest/conftest.py
    dom_dir = simulation["ws"]
    control = simulation["control_file"]
    sim_key = f"{dom_dir.name}:{control.with_suffix('').name}"
    sim_tag = sim_key.replace(":", "_")
    test_data_version_file = pl.Path(f"../.test_data_version_{sim_tag}.txt")
    if test_data_version_file.exists():
        test_data_version_file.unlink()

    return
