import pathlib as pl
from typing import Union

# A module for path/file utilities


def path_rel_to_yml(
    file_in_yml: Union[pl.Path, str], yml: Union[pl.Path, str]
):
    """Resolve a path from a yaml file

    Given a yaml file (yml) and a file specified within that yaml file,
    if the file is an absolute path, return it as a pathlib.Path object,
    otherwise resolve the file path relative to the location of the yaml file.

    Args:
        file_in_yml: a str or pathlib.Path from within a yaml file
        yml: the path of the yaml file.

    Return:
        pathlib.Path object with resolved/absolute path
    """
    yml_pl = pl.Path(yml)
    file_pl = pl.Path(file_in_yml)
    if not file_pl.is_absolute():
        file_pl = (yml_pl.parent / file_pl).resolve()
    return file_pl


def assert_exists(path):
    assert pl.Path(path).exists()
    return
