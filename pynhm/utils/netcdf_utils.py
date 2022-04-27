import pathlib as pl
from typing import Union

import netCDF4 as nc4
import numpy as np

fileish = Union[str, pl.Path]
listish = Union[list, tuple]
arrayish = Union[list, tuple, np.ndarray]
ATOL = np.finfo(np.float32).eps


class NetCdfRead:
    def __init__(
        self,
        name: fileish,
        nc_read_vars: list = None,
    ) -> "NetCdfRead":

        self._nc_file = name
        self._nc_read_vars = nc_read_vars
        self._open_nc_file()
        self._itime_step = {}
        for variable in self.variables:
            self._itime_step[variable] = 0

    def __del__(self):
        self.close()

    def close(self):
        if self.dataset.isopen():
            self.dataset.close()

    def _open_nc_file(self):
        self.dataset = nc4.Dataset(self._nc_file, "r")
        self.ds_var_list = list(self.dataset.variables.keys())
        if self._nc_read_vars is None:
            self._nc_read_vars = self.ds_var_list
        self.ds_var_chunking = {
            vv: self.dataset.variables[vv].chunking()
            for vv in self.ds_var_list
        }
        # Set dimension variables which are not chunked
        self._datetime = (
            nc4.num2date(
                self.dataset.variables["datetime"][:],
                units=self.dataset.variables["datetime"].units,
                calendar="standard",
                only_use_cftime_datetimes=False,
            )
            .filled()
            .astype("datetime64[s]")
            # JLM: the global time type as in cbh_utils, define somewhere
        )
        self._ntimes = self._datetime.shape[0]
        spatial_id_name = (
            "nhm_id" if "nhm_id" in self.ds_var_list else "hru_id"
        )
        self._spatial_id = self.dataset.variables[spatial_id_name][:]

        self._variables = [
            name
            for name in self.ds_var_list
            if name != "datetime" and name != spatial_id_name
        ]

    @property
    def ntimes(
        self,
    ) -> int:
        """Get number of times in the netcdf file

        Returns:
            ntimes: number of times in the NetCDF file

        """
        return self._ntimes

    @property
    def date_times(
        self,
    ) -> np.ndarray:
        """Get the datetimes in the NetCDF file

        Returns:
            data_times: numpy array of datetimes in the NetCDF file

        """
        return self._datetime

    @property
    def nhru(
        self,
    ) -> int:
        """Get number of HRUs in the NetCDF file

        Returns:
            nhru: number of HRUs in the NetCDF file

        """
        return self._spatial_id.shape[0]

    @property
    def spatial_ids(
        self,
    ) -> np.ndarray:
        """Get the spatial IDs in the NetCDF file

        Returns:
            arr: numpy array with the spatial IDs in the NetCDF file

        """
        return self._spatial_id

    @property
    def variables(self):
        """Get a list of variable names

        Returns:
            variables: list of variable names, excluding the datetime and
                nhru variables

        """
        return self._variables

    def get_data(
        self,
        variable: str,
        itime_step: int = None,
    ) -> np.ndarray:
        """Get data for a variable

        Args:
            variable: variable name
            itime_step: time step to return. If itime_step is None all of the
              data for a variable is returned

        Returns:
            arr: numpy array with the data for a variable

        """
        if variable not in self._nc_read_vars:
            raise ValueError(
                f"'{variable}' not in list of available variables"
            )
        if itime_step is None:
            return self.dataset[variable][:, :]
        else:
            if itime_step >= self._ntimes:
                raise ValueError(
                    f"requested time step {itime_step} but only "
                    + f"{self._ntimes} time steps are available."
                )
            return self.dataset[variable][itime_step, :]

    def advance(
        self,
        variable: str,
    ) -> np.ndarray:
        """Get the data for a variable for the next time step

        Args:
            variable: variable name

        Returns:
            arr: numpy array with the data for a variable for the current
                time step

        """
        arr = self.get_data(
            variable,
            itime_step=self._itime_step[variable],
        )
        self._itime_step[variable] += 1
        return arr


class NetCdfWrite:
    """Output the csv output data to a netcdf file

    Args:
        name: path for netcdf output file
        clobber: boolean indicating if an existing netcdf file should
            be overwritten
        zlib: boolean indicating if the data should be compressed
            (default is True)
        complevel: compression level (default is 4)
        chunk_sizes: dictionary defining chunk sizes for the data
    """

    def __init__(
        self,
        name: fileish,
        hru_ids: arrayish,
        variables: listish,
        var_meta: dict,
        time_units: str = "days since 1970-01-01 00:00:00",
        clobber: bool = True,
        zlib: bool = True,
        complevel: int = 4,
        chunk_sizes: dict = {"time": 30, "hruid": 0},
    ) -> "NetCdfWrite":

        self.dataset = nc4.Dataset(name, "w", clobber=clobber)
        self.dataset.setncattr("Description", "PYNHM output data")

        if isinstance(hru_ids, (list, tuple)):
            nhrus = len(hru_ids)
        elif isinstance(hru_ids, np.ndarray):
            nhrus = hru_ids.shape[0]
        self.nhrus = nhrus
        self.hru_ids = hru_ids

        # Dimensions
        # None for the len argument gives an unlimited dim
        self.dataset.createDimension("time", None)
        self.dataset.createDimension("hru_id", self.nhrus)

        # Dim Variables
        self.time = self.dataset.createVariable("datetime", "f4", ("time",))
        self.time.units = time_units

        self.hruid = self.dataset.createVariable("hru_id", "i4", ("hru_id"))
        self.hruid[:] = np.array(hru_ids, dtype=int)

        self.variables = {}
        for vv in variables:
            variabletype = "f4"
            self.variables[vv] = self.dataset.createVariable(
                vv,
                variabletype,
                ("time", "hru_id"),
                fill_value=nc4.default_fillvals[variabletype],
                zlib=zlib,
                complevel=complevel,
                chunksizes=tuple(chunk_sizes.values()),
            )
            for key, val in var_meta.items():
                if isinstance(val, dict):
                    continue
                self.variables[vv].setncattr(key, val)

    def __del__(self):
        self.close()

    def close(self):
        if self.dataset.isopen():
            self.dataset.close()

    def add_simulation_time(self, itime_step: int, simulation_time: float):
        # var = self.variables["datetime"]
        # var[itime_step] = simulation_time
        self.time[itime_step] = simulation_time

    def add_data(
        self, name: str, itime_step: int, current: np.ndarray
    ) -> None:
        """Add data for a time step to a NetCDF variable

        Args:
            name:
            itime_step:

        Returns:

        """
        if name not in self.variables.keys():
            raise KeyError(f"{name} not a valid variable name")
        var = self.variables[name]
        var[itime_step, :] = current[:]

    def add_all_data(self, name: str, current: np.ndarray) -> None:
        """Add data to a NetCDF variable

        Args:
            name:

        Returns:

        """
        if name not in self.variables.keys():
            raise KeyError(f"{name} not a valid variable name")
        var = self.variables[name]
        var[:, :] = current[:, :]


class NetCdfCompare:
    """Compare variables in two NetCDF files

    Args:
        base_pth: path to base NetCDF file
        compare_pth: path to NetCDF file that will be compared to
            the base NetCDF file
        verbose: boolean determining if information should be written
            to stdout
    """

    def __init__(
        self,
        base_pth: fileish,
        compare_pth: fileish,
        verbose: bool = False,
    ) -> None:
        self._base = NetCdfRead(base_pth)
        self._compare = NetCdfRead(compare_pth)
        self._verbose = verbose
        if not self.__validate_comparison:
            raise KeyError(
                f"Variables in {compare_pth} do not exist in {base_pth}"
            )

    def __del__(self):
        if self._base.dataset.isopen():
            self._base.dataset.close()
        if self._compare.dataset.isopen():
            self._compare.dataset.close()

    def compare(
        self,
        atol: float = ATOL,
    ) -> (bool, dict):
        """

        Args:
            atol: absolute tolerance to use in numpy allclose (default
                is single precision machine precisions ~1e-7)

        Returns:
            success: boolean indicating if comparisons for all variables
                in the NetCDF file passed
            failures: dictionary with the maximum absolute difference for
                each variable in the NetCDF file.

        """
        return self.__compare(atol=atol)

    def __compare(
        self,
        atol: float = ATOL,
    ):
        success = True
        failures = {}
        for variable in self._compare.variables:
            base_arr = self._base.get_data(variable)
            compare_arr = self._compare.get_data(variable)
            if not np.allclose(base_arr, compare_arr, atol=atol):
                success = False
                failures[variable] = np.max(np.abs(compare_arr - base_arr))
        return success, failures

    @property
    def __validate_comparison(self) -> bool:
        """Validate the variables in the comparison NetCDF file relative
        to the base NetCDF file


        Returns:
            valid: boolean indicating if all variables in the comparison
                NetCDF file are in the base NetCDF file

        """
        valid = False
        for variable in self._compare.variables:
            if variable in self._base.variables:
                valid = True
                break
        return valid
