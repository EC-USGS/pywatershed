import numpy as np


class DataAccess:
    def __init__(self):
        self.name = "DataAccess"
        self._coords = []
        self._variables = []
        self._potential_variables = []
        return

    @property
    def coords(self) -> list:
        return self._coords

    @property
    def variables(self) -> list:
        return self._variables

    def _has_potential_state(self, state_name, fail=True) -> bool:
        if state_name in self._potential_variables:
            return True
        return False

    def _has_state(self, state_name, fail=True) -> bool:
        if state_name in self._variables + self._coords:
            return True
        return False

    def _has_coord(self, coord_name, fail=True) -> bool:
        if coord_name in self._coords:
            return True
        return False

    def __setitem__(self, name: str, value: np.ndarray) -> None:
        if not isinstance(value, np.ndarray):
            raise TypeError(
                f"Value to set must be of type np.ndarray not '{type(value)}'"
            )
        if self._has_potential_state(name):
            setattr(self, name, value)
            self._variables = list(set(self._variables + [name]))
        elif self._has_coord(name):
            # only set coord if it is None.
            if hasattr(self, name):
                msg = f"Coordinate variable can not be re-set: '{name}'"
                raise KeyError(msg)
            setattr(self, name, value)
        else:
            msg = (
                f"{self.name} does not have potential state or coordinate "
                f"'{name}'"
            )
            raise KeyError(msg)
        return None

    def __getitem__(self, name: str) -> np.ndarray:
        if self._has_state(name) or self._has_coord(name):
            return getattr(self, name)
        else:
            msg = f"'{name}' not in state or coordinate variables"
            raise KeyError(msg)

    def __delitem__(self, name: str):
        if self._has_state(name):
            if name in self._variables:
                _ = self._variables.remove(name)
                return getattr(self, name)
            else:
                raise KeyError("Can not delete coordinate '{name}'")
        elif self._has_coord(name):
            msg = f"Can not delete coordinate '{name}'"
            raise KeyError(msg)
        else:
            msg = f"'{name}' not in state or coordinate variables"
            raise KeyError(msg)
        return None
