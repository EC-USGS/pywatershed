import datetime as dt
import pathlib as pl
from typing import Union

import netCDF4 as nc4
import numpy as np
import pandas as pd

from .cbh_metadata import cbh_metadata

fileish = Union[str, pl.PosixPath]


class CsvFile:
    def __init__(
        self,
        name: fileish = None,
        convert: bool = False,
    ) -> "CsvFile":
        self.paths = {}
        if name is not None:
            self._add_path(name)
        self.convert = convert
        self._data = None

    @property
    def hru_ids(self) -> list:
        """Get a list of hru ids

        Returns:
            hru_ids: list of hru_ids

        """
        if self._data is None:
            self._get_data()

        hru_ids = []
        for name in self._data.dtype.names[1:]:
            hru_id = name.split("_")[-1]
            if hru_id not in hru_ids:
                hru_ids.append(int(hru_id))
        return hru_ids

    @property
    def variable_names(self) -> list:
        """Get a list of unique variables

        Returns:
            variables: list of variables

        """
        if self._data is None:
            self._get_data()

        variables = []
        for name in self._data.dtype.names[1:]:
            hru_id = name.split("_")[-1]
            variable = name.replace(f"_{hru_id}", "")
            if variable not in variables:
                variables.append(variable)
        return variables

    @property
    def data(self) -> np.recarray:
        """Get csv output data as a numpy recarray

        Returns:
            data : numpy recarray containing all of the csv data

        """
        if self._data is None:
            self._get_data()
        return self._data

    def add_path(
        self,
        name: fileish,
    ) -> None:
        """Add a csv output file path to the object

        Args:
            name: path for csv output file

        Returns:
            None

        """

        self._add_path(name)

    def to_dataframe(self) -> pd.DataFrame:
        """Get the csv output data as a pandas dataframe

        Returns:
            df: csv output data as a pandas dataframe

        """
        if self._data is None:
            self._get_data()

        df = pd.DataFrame(self._data).set_index("date")
        return df

    def to_netcdf(
        self,
        name: fileish,
        global_atts: dict = None,
        clobber: bool = True,
        zlib: bool = True,
        complevel: int = 4,
        chunk_sizes: dict = {"time": 30, "hruid": 0},
    ) -> None:
        """Output the csv output data to a netcdf file

        Args:
            name: path for netcdf output file
            clobber: boolean indicating if an existing netcdf file should
                be overwritten
            zlib: boolean indicating if the data should be compressed
                (default is True)
            complevel: compression level (default is 4)
            chunk_sizes: dictionary defining chunk sizes for the data

        Returns:
            None

        """
        if self._data is None:
            self._get_data()

        hru_ids = self.hru_ids
        variables = self.variable_names
        ntimes = self._data.shape[0]
        nhrus = len(hru_ids)

        ds = nc4.Dataset(name, "w", clobber=clobber)
        ds.setncattr("Description", "PRMS output data")
        if global_atts is not None:
            for key, val in global_atts.items():
                ds.setncattr(key, val)

        # Dimensions
        # None for the len argument gives an unlimited dim
        ds.createDimension("time", ntimes)
        ds.createDimension("hru_id", nhrus)

        # Dim Variables
        time = ds.createVariable("datetime", "f4", ("time",))
        start_date = self._data["date"][0].strftime("%Y-%m-%d %H:%M:%S")
        time_units = f"days since {start_date}"
        time[:] = nc4.date2num(
            self._data["date"].astype(dt.datetime),
            units=time_units,
            calendar="standard",
        )

        hruid = ds.createVariable("hru_id", "i4", ("hru_id"))
        hruid[:] = np.array(hru_ids, dtype=int)

        # Variables
        for vv in variables:
            vvtype = "f4"
            var = ds.createVariable(
                vv,
                vvtype,
                ("time", "hru_id"),
                fill_value=nc4.default_fillvals[vvtype],  # JLM: sus
                zlib=zlib,
                complevel=complevel,
                chunksizes=tuple(chunk_sizes.values()),
            )
            arr = np.zeros((ntimes, nhrus), dtype=np.float32)
            for idx, hru in enumerate(hru_ids):
                key = f"{vv}_{hru}"
                arr[:, idx] = self._data[key][:]
            ds.variables[vv][:, :] = arr

        ds.close()
        print(f"Wrote netcdf file: {name}")
        return

    def _add_path(
        self,
        name: (pl.Path, str),
    ):
        if isinstance(name, str):
            name = pl.Path(name)
        elif not isinstance(name, pl.Path):
            raise TypeError(f"{name} must be a string or pathlib.Path object")
        self.paths[name.name] = name

    def _get_data(self) -> None:
        """Read csv data into a single numpy recarray

        Returns:
            None

        """
        str2date = lambda x: dt.datetime.strptime(
            x.decode("utf-8"), "%Y-%m-%d"
        )
        all_data = []
        ntimes = 0
        dtype = [("date", dt.datetime)]
        for key, path in self.paths.items():
            if path.exists():
                try:
                    arr = np.genfromtxt(
                        path,
                        dtype=None,
                        names=True,
                        delimiter=",",
                        converters={0: str2date},
                    )
                except:
                    raise IOError(f"numpy could not parse...'{path}'")

            ntimes = max(arr.shape[0], ntimes)
            column_base = path.stem
            column_names = [
                f"{column_base}_{hru_id.strip()}"
                for hru_id in arr.dtype.names[1:]
            ]
            arr.dtype.names = ["date"] + column_names

            # add additional column names to the dtype
            for name in column_names:
                dtype.append((name, np.float32))

            all_data.append(arr)

        self._data = np.zeros(ntimes, dtype=dtype)
        for idx, arr in enumerate(all_data):
            if idx == 0:
                i0 = 0
            else:
                i0 = 1
            for name in arr.dtype.names[i0:]:
                self._data[name][:] = arr[name][:]
