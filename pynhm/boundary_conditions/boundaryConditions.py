import pandas as pd
import numpy as np

class BoundaryConditions:
    def __init__(self) -> "BoundaryConditions":
        self.boundary_conditions = {}
        self.current = {}

    def add_boundary_dataframe(self, boundary_name: str, boundary_condition: pd.DataFrame,):
        self.boundary_conditions[boundary_name] = boundary_condition.to_numpy()

    def add_boundary_recarray(self, boundary_name: str, boundary_condition: np.recarray,):
        self.boundary_conditions[boundary_name] = boundary_condition

    def advance(self, itime_step: int,) -> None:
        for idx, (name, values) in enumerate(self.boundary_conditions.items()):
            arr = np.zeros(len(values.dtype.names[1:]), dtype=float)
            for jdx, key in enumerate(values.dtype.names[1:]):
                arr[jdx] = values[key][itime_step]
            self.current[name] = arr

    def set_pointers(self, component: object) -> None:
        keys = tuple(self.current.keys())
        for key in dir(component):
            if not key.startswith("_") and key in keys:
                setattr(component, key, self.current[key])


