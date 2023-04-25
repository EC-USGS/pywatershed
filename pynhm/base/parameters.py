import numpy as np

from .data_model import DatasetDict


class Parameters(DatasetDict):
    def __init__(
        self,
        dims: dict = {},
        coords: dict = {},
        data_vars: dict = {},
        metadata: dict = {},
        encoding: dict = {},
        validate: bool = True,
    ) -> "Parameters":
        """Parameter class

        This is a subclass of data_model.DatasetDict, it has all the same
        methods. New methods map to DatasetDict as follows:
            parameters: dd.variables
            get_param_values: get values from dd.variables
            get_dim_values: get values from dd.dims

        """
        super().__init__(
            dims=dims,
            coords=coords,
            data_vars=data_vars,
            metadata=metadata,
            encoding=encoding,
            validate=validate,
        )
        return

    @property
    def parameters(self) -> dict:
        return self.variables

    def get_param_values(
        self,
        keys: list | str = None,
    ) -> dict | np.ndarray:
        """Get the values of the parameters (coords or data_vars) by keys

        Also see:
            subset() method is a Parameter object is desired.
        """
        if not isinstance(keys, list):
            return self.variables[keys]
        else:
            return {kk: self.variables[kk] for kk in keys}

    def get_dim_values(
        self,
        keys: list | tuple | str = None,
    ) -> dict | np.ndarray:
        """Get the values of the dimensions by keys."""
        if not isinstance(keys, (list, tuple)):
            return self.dims[keys]
        else:
            return {kk: self.dims[kk] for kk in keys}
