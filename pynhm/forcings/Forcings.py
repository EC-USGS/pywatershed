import datetime

try:
    import pandas as pd
    from pandas import DataFrame
except ModuleNotFoundError:
    pd = None
    DataFrame = None


class Forcings:
    def __init__(
        self,
        name: str = None,
        data: DataFrame = None,
        verbose: int = 0,
    ):
        self.verbose = verbose

        if name is None or data is None:
            self.forcings = {}
        else:
            self.forcings[name] = data
        self.current = {}
        self.itime_step = None
        self.current_date = None

    def add_forcing(
        self,
        name: str,
        data: DataFrame,
    ):
        if name in self.forcings.keys():
            raise KeyError(f"'{name}' key already exists in Forcings object")
        self.forcings[name] = data

    def advance(
        self,
        itime_step: int,
        current_date: datetime,
    ):
        self.itime_step = itime_step
        self.current_date = current_date
        dt = pd.to_datetime(current_date)
        for key, df in self.forcings.items():
            self.current[key] = df.iloc[df.index.get_loc(dt, method="nearest")]
