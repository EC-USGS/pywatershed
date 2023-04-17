from .data_model import DatasetDict
from ..constants import listish


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

    def get_parameters(
        self,
        keys: listish = None,
        copy=False,
        keep_global: bool = False,
        keep_global_metadata: bool = None,
        keep_global_encoding: bool = None,
        process=None,
    ) -> "Parameters":
        return self.subset(
            keys=keys,
            copy=copy,
            keep_global=keep_global,
            keep_global_metadata=keep_global_metadata,
            keep_global_encoding=keep_global_encoding,
            process=process,
        )
