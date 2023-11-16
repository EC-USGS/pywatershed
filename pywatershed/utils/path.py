import pathlib as pl
from typing import Union

# A module for path/file utilities


def path_rel_to_yaml(
    file_in_yaml: Union[pl.Path, str], yaml: Union[pl.Path, str]
):
    """Resolve a path from a yaml file

    Given a yaml file (yaml) and a file specified within that yaml file,
    if the file is an absolute path, return it as a pathlib.Path object,
    otherwise resolve the file path relative to the location of the yaml file.

    Args:
        file_in_yaml: a str or pathlib.Path from within a yaml file
        yaml: the path of the yaml file.

    Return:
        pathlib.Path object with resolved/absolute path
    """
    yaml_pl = pl.Path(yaml)
    file_pl = pl.Path(file_in_yaml)
    if not file_pl.is_absolute():
        file_pl = (yaml_pl.parent / file_pl).resolve()
    return file_pl


def assert_exists(path):
    assert pl.Path(path).exists()
    return


def dict_pl_to_str(the_dict):
    """Convert dict items of class  pathlib.Path to strings, recursively"""
    for key, val in the_dict.items():
        if isinstance(val, dict):
            the_dict[key] = dict_pl_to_str(val)
        elif isinstance(val, pl.Path):
            the_dict[key] = str(val)

    return the_dict
