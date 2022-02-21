import pathlib as pl

import pandas as pd


class CsvFile:
    def __init__(
        self,
        name: (pl.Path, str) = None,
        convert: bool = False,
    ):
        self.paths = {}
        if name is not None:
            self._add_path(name)
        self.convert = convert

    def add_path(
        self,
        name: (pl.Path, str),
    ):
        self._add_path(name)

    def to_dataframe(self):
        templist = []
        for key, path in self.paths.items():
            if path.exists():
                try:
                    tdf = pd.read_csv(
                        path,
                        parse_dates=["Date"],
                        index_col=["Date"],
                        dtype=float,
                    )
                except:
                    print(f"pandas could not parse...'{path}'")
                    continue
            column_base = path.stem

            column_names = [
                f"{column_base}_{hru_id.strip()}" for hru_id in tdf.columns
            ]
            tdf.columns = [column_names]
            tdf.index.names = ["date"]
            if self.convert:
                raise NotImplementedError("conversion not implemented")
            templist.append(tdf)

        # concatenate individual dataframes
        return pd.concat([v for v in templist], axis=1)

    def _add_path(
        self,
        name: (pl.Path, str),
    ):
        if isinstance(name, str):
            name = pl.Path(name)
        elif not isinstance(name, pl.Path):
            raise TypeError(f"{name} must be a string or pathlib.Path object")
        self.paths[name.name] = name
