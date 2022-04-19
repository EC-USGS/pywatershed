import pathlib as pl
from typing import Union, Iterable

from ..constants import __pynhm_root__

import yaml

fileish = Union[str, pl.Path]


def load_yaml_file(the_file):
    with pl.Path(the_file).open("r") as file_stream:
        data = yaml.load(file_stream, Loader=yaml.Loader)
    return data


class Meta:
    def __init__(
        self,
        dimensions: fileish = (
            __pynhm_root__ / "static/metadata/dimensions.yaml"
        ),
        control: fileish = (__pynhm_root__ / "static/metadata/control.yaml"),
        parameters: fileish = (
            __pynhm_root__ / "static/metadata/parameters.yaml"
        ),
        variables: fileish = (
            __pynhm_root__ / "static/metadata/variables.yaml"
        ),
    ):

        self.name = "Meta"

        self.dimensions = load_yaml_file(dimensions)
        self.control = load_yaml_file(control)
        self.parameters = load_yaml_file(parameters)
        self.variables = load_yaml_file(variables)

        return

    def get_var_subclass(self, subclass_name: str) -> dict:
        """Get the variable metadata for the supplied subclass."""
        return {
            key: value
            for key, value in self.variables.items()
            if subclass_name in value["modules"]
        }

    def get_inputs_subclass(self, subclass_name: str) -> dict:
        """Get the inputs metadata for the supplied subclass."""
        return {
            key: value
            for key, value in self.variables.items()
            if ("input_to" in value) and (subclass_name in value["input_to"])
        }

    def _get_meta_in_list(self, meta: str, the_list: Iterable) -> dict:
        return {
            key: value
            for key, value in getattr(self, meta).items()
            if key in the_list
        }

    def get_dims(self, var_list: Iterable) -> dict:
        """Get the variable metadata for the supplied subclass."""
        return self._get_meta_in_list("dimensions", var_list)

    def get_control(self, var_list: Iterable) -> dict:
        """Get the variable metadata for the supplied subclass."""
        return self._get_meta_in_list("control", var_list)

    def get_params(self, inputs_list: Iterable) -> dict:
        """Get the variable metadata for the supplied subclass."""
        return self._get_meta_in_list("parameters", inputs_list)

    def get_vars(self, var_list: Iterable) -> dict:
        """Get the variable metadata for the supplied subclass."""
        return self._get_meta_in_list("variables", var_list)
