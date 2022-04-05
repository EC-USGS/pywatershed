import pathlib as pl
from typing import Union

from ..utils.netcdf_utils import NetCdfRead

fileish = Union[str, pl.Path]


class BoundaryConditions:
    """Boundary Condition class"""

    def __init__(self) -> "BoundaryConditions":
        self.variables = {}
        self.current = {}

    def add_boundary(
        self,
        name: fileish,
    ) -> None:
        ncf = NetCdfRead(name)
        for variable in ncf.variables:
            self.variables[variable] = ncf

    def advance(
        self,
    ) -> None:
        for variable, ncf in self.variables.items():
            self.current[variable] = ncf.advance(variable)

    def set_pointers(self, component: object) -> None:
        for key in dir(component):
            if not key.startswith("_") and key in self.variables.keys():
                setattr(component, key, self.current[key])
