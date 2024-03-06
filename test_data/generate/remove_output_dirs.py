import pathlib as pl
import shutil


def test_remove_output_dirs(simulation):
    output_dir = pl.Path(simulation["output_dir"])
    if output_dir.exists():
        shutil.rmtree(output_dir)

    # this shouldnt matter but it is necessary to exist for PRMS to write to it
    output_dir.mkdir(parents=True)

    return
