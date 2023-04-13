from .data_model import DatasetDict


class Parameters(DatasetDict):
    @property
    def parameters(self) -> dict:
        return self._data_vars
