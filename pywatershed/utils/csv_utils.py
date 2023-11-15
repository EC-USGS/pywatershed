import datetime as dt
import pathlib as pl
from typing import Union

import netCDF4 as nc4
import numpy as np
import pandas as pd

from ..base import meta

fileish = Union[str, pl.PosixPath, dict]


class CsvFile:
    """CSV file object
    path: a string, pathlib.Path or dict. The key of the dict can be used
          to rename the variable in the recarray and upon output to Netcdf.
          The value of the dict should be a string or a pathlib.Path. Only
          dicts of len 1 are allowed currently.

    """

    def __init__(
        self,
        path: fileish = None,
        convert: bool = False,
    ) -> "CsvFile":
        self.paths = {}
        if path is not None:
            self._add_path(path)
        self.convert = convert
        self._variables = None
        self._coordinates = None
        self._data = None
        self.meta = meta

    @property
    def nhm_id(self) -> list:
        """Get a list of nhm ids

        Returns:
            nhm_ids: list of  nhm ids

        """
        self._lazy_data_evaluation()

        key = "nhm_id"
        if key in self._coordinates.keys():
            nhm_ids = self._coordinates[key]
        else:
            nhm_ids = None
        return nhm_ids

    @property
    def nhm_seg(self) -> list:
        """Get a list of nhm segments

        Returns:
            nhm_segs: list of nhm segments

        """
        self._lazy_data_evaluation()

        key = "nhm_seg"
        if key in self._coordinates.keys():
            nhm_segs = self._coordinates[key]
        else:
            nhm_segs = None
        return nhm_segs

    @property
    def variable_names(self) -> list:
        """Get a list of unique variables

        Returns:
            variables: list of variables

        """
        self._lazy_data_evaluation()
        return self._variables

    @property
    def data(self) -> np.recarray:
        """Get csv output data as a numpy recarray

        Returns:
            data : numpy recarray containing all of the csv data

        """
        self._lazy_data_evaluation()
        return self._data

    def add_path(
        self,
        path: fileish,
    ) -> None:
        """Add a csv output file path to the object

        Args:
            name: path for csv output file

        Returns:
            None

        """

        self._add_path(path)

    def to_dataframe(self) -> pd.DataFrame:
        """Get the csv output data as a pandas dataframe

        Returns:
            df: csv output data as a pandas dataframe

        """
        self._lazy_data_evaluation()
        df = pd.DataFrame(self._data).set_index("date")
        return df

    def to_netcdf(
        self,
        name: fileish,
        global_atts: dict = None,
        clobber: bool = True,
        zlib: bool = True,
        complevel: int = 4,
        chunk_sizes: dict = {"time": 30, "nhm_id": 0, "nhm_seg": 0},
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
        self._lazy_data_evaluation()

        ds = nc4.Dataset(name, "w", clobber=clobber)
        ds.setncattr("Description", "PRMS output data")
        if global_atts is not None:
            for key, val in global_atts.items():
                ds.setncattr(key, val)

        # Dimensions
        # None for the len argument gives an unlimited dim
        ntimes = self._data.shape[0]
        ds.createDimension("time", ntimes)
        for key, value in self._coordinates.items():
            ds.createDimension(key, len(value))

        # Dim Variables
        time = ds.createVariable("time", "f4", ("time",))
        start_date = self._data["date"][0].strftime("%Y-%m-%d %H:%M:%S")
        time_units = f"days since {start_date}"
        time.units = time_units
        time[:] = nc4.date2num(
            self._data["date"].astype(dt.datetime),
            units=time_units,
            calendar="standard",
        )

        for key, value in self._coordinates.items():
            coord_id = ds.createVariable(key, "i4", (key))
            coord_id[:] = np.array(value, dtype=int)

        dimensions = self.meta.get_dimensions(self.variable_names)

        # Variables
        for variable_name in self.variable_names:
            if self.meta.is_available(variable_name):
                meta_dict = self.meta.find_variables(variable_name)
                variable_type = meta.meta_netcdf_type(meta_dict[variable_name])
                dimension = dimensions[variable_name]
                if "nsegment" in dimension:
                    dim_name = "nhm_seg"
                else:
                    dim_name = "nhm_id"
                dtype = meta.meta_numpy_type(meta_dict[variable_name])
            else:
                variable_type = "f4"
                dim_name = "nhm_id"
                dtype = np.float32

            ids = self._coordinates[dim_name]
            nids = len(ids)

            dims = ("time", dim_name)
            chunk_sizes_var = [chunk_sizes[vv] for vv in dims]

            _ = ds.createVariable(
                variable_name,
                variable_type,
                dims,
                fill_value=nc4.default_fillvals[variable_type],  # JLM: sus
                zlib=zlib,
                complevel=complevel,
                chunksizes=tuple(chunk_sizes_var),
            )
            # add additional meta data
            if self.meta.is_available(variable_name):
                var_meta = self.meta.find_variables(variable_name)[
                    variable_name
                ]
                for key, val in var_meta.items():
                    if isinstance(val, dict):
                        continue
                    ds.variables[variable_name].setncattr(key, val)

            arr = np.zeros((ntimes, nids), dtype=dtype)
            for idx, on_id in enumerate(ids):
                key = f"{variable_name}_{on_id}"
                arr[:, idx] = self._data[key][:]
            ds.variables[variable_name][:, :] = arr

        ds.close()
        print(f"Wrote netcdf file: {name}")
        return

    def _add_path(
        self,
        path: fileish,
    ):
        if isinstance(path, (str, pl.Path)):
            path = pl.Path(path)
            self.paths[path.stem] = path
        elif isinstance(path, dict):
            if len(path) > 1:
                raise ValueError("Only dicts of len 1 allowed currently")
            for key, val in path.items():
                self.paths[key] = pl.Path(val)
        elif not isinstance(path, pl.Path):
            raise TypeError("path must be a string or pathlib.Path object")

    def _lazy_data_evaluation(self):
        if self._data is None:
            self._get_data()

    def _get_data(self) -> None:
        """Read csv data into a single numpy recarray

        Returns:
            None

        """

        def str2date(x):
            return dt.datetime.strptime(x.decode("utf-8"), "%Y-%m-%d")

        all_data = []
        ntimes = 0
        dtype = [("date", dt.datetime)]
        for variable_name, path in self.paths.items():
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

            if self._variables is None:
                self._variables = [variable_name]
            else:
                self._variables.append(variable_name)

            # determine the variable type
            if self.meta.is_available(variable_name):
                variable_type = meta.meta_numpy_type(
                    self.meta.find_variables(variable_name)[variable_name]
                )
                if (
                    "nsegment"
                    in self.meta.get_dimensions(variable_name)[variable_name]
                ):
                    coordinate_name = "nhm_seg"
                else:
                    coordinate_name = "nhm_id"
            else:
                variable_type = np.float32
                coordinate_name = "nhm_id"

            # set coordinates
            if self._coordinates is None:
                self._coordinates = {}
            if coordinate_name not in list(self._coordinates.keys()):
                self._coordinates[coordinate_name] = [
                    idx for idx in arr.dtype.names[1:]
                ]

            column_names = [
                f"{variable_name}_{idx.strip()}" for idx in arr.dtype.names[1:]
            ]
            arr.dtype.names = ["date"] + column_names

            # add additional column names to the dtype
            for name in column_names:
                dtype.append((name, variable_type))

            all_data.append(arr)

            # reset ntimes, if necessary
            ntimes = max(arr.shape[0], ntimes)

        self._data = np.zeros(ntimes, dtype=dtype)
        for idx, arr in enumerate(all_data):
            if idx == 0:
                i0 = 0
            else:
                i0 = 1
            for name in arr.dtype.names[i0:]:
                self._data[name][:] = arr[name][:]
