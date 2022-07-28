import pathlib as pl
from typing import Union

import netCDF4 as nc4
import numpy as np

from ..base.accessor import Accessor
from ..base.meta import meta_dimensions, meta_netcdf_type
from ..utils.time_utils import datetime_doy

fileish = Union[str, pl.Path]
listish = Union[list, tuple]
arrayish = Union[list, tuple, np.ndarray]
ATOL = np.finfo(np.float32).eps

# JLM TODO: start_time: np.datetime64 = None ?
#     Still need some mechanism for initial itime_step for datetime coordinate
#     variables


class NetCdfRead(Accessor):
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

        if "datetime" in self.dataset.variables:
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

        if "doy" in self.dataset.variables:
            self._doy = self.dataset.variables["doy"][:].data
            self._ntimes = self._doy.shape[0]

        spatial_id_names = []
        if "nhm_id" in self.ds_var_list:
            spatial_id_names.append("nhm_id")
        elif "hru_id" in self.ds_var_list:
            spatial_id_names.append("hru_id")
        if "nhm_seg" in self.ds_var_list:
            spatial_id_names.append("nhm_seg")
        elif "hru_seg" in self.ds_var_list:
            spatial_id_names.append("hru_seg")

        # set default spatial id
        if len(spatial_id_names) < 1:
            spatial_id_names.append("hru_id")

        # set spatial id dictionary
        self._spatial_ids = {}
        for spatial_id_name in spatial_id_names:
            self._spatial_ids[spatial_id_name] = self.dataset.variables[
                spatial_id_name
            ][:]

        self._variables = [
            name
            for name in self.ds_var_list
            if name != "datetime" and name not in spatial_id_names
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
            nhru_shape: number of HRUs in the NetCDF file

        """
        hru_key = None
        if "nhm_id" in self._spatial_ids.keys():
            hru_key = "nhm_id"
        elif "hru_id" in self._spatial_ids.keys():
            hru_key = "hru_id"
        if hru_key is not None:
            nhru_shape = self._spatial_ids[hru_key].shape[0]
        else:
            nhru_shape = 0
        return nhru_shape

    @property
    def nsegment(
        self,
    ) -> int:
        """Get number of segments in the NetCDF file

        Returns:
            seg_shape: number of segments in the NetCDF file

        """
        seg_key = None
        if "nhm_seg" in self._spatial_ids.keys():
            seg_key = "nhm_seg"
        if seg_key is not None:
            seg_shape = self._spatial_ids[seg_key].shape[0]
        else:
            seg_shape = 0
        return seg_shape

    @property
    def spatial_ids(
        self,
    ) -> np.ndarray:
        """Get the spatial IDs in the NetCDF file

        Returns:
            arr: numpy array with the spatial IDs in the NetCDF file

        """
        return self._spatial_ids

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
        self, variable: str, current_time: np.datetime64 = None
    ) -> np.ndarray:
        """Get the data for a variable for the next time step

        Args:
            variable: variable name

        Returns:
            arr: numpy array with the data for a variable for the current
                time step

        """
        if "datetime" in self.dataset.variables:
            arr = self.get_data(
                variable,
                itime_step=self._itime_step[variable],
            )

        if "doy" in self.dataset.variables:
            arr = self.get_data(
                variable,
                itime_step=datetime_doy(current_time) - 1,
            )

        self._itime_step[variable] += 1

        return arr


class NetCdfWrite(Accessor):
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
        coordinates: dict,
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

        variable_dimensions = {}
        nhru_coordinate = False
        nsegment_coordinate = False
        for var_name in variables:
            dimension_name = meta_dimensions(var_meta[var_name])
            variable_dimensions[var_name] = dimension_name
            if (
                "nhru" in dimension_name
                or "ngw" in dimension_name
                or "nssr" in dimension_name
            ):
                nhru_coordinate = True
            if "nsegment" in dimension_name:
                nsegment_coordinate = True

        if nhru_coordinate:
            hru_ids = coordinates["nhm_id"]
            if isinstance(hru_ids, (list, tuple)):
                nhrus = len(hru_ids)
            elif isinstance(hru_ids, np.ndarray):
                nhrus = hru_ids.shape[0]
            self.nhrus = nhrus
            self.hru_ids = hru_ids

        if nsegment_coordinate:
            hru_segments = coordinates["nhm_seg"]
            if isinstance(hru_segments, (list, tuple)):
                nsegments = len(hru_segments)
            elif isinstance(hru_segments, np.ndarray):
                nsegments = hru_segments.shape[0]
            self.nsegments = nsegments
            self.hru_segments = hru_segments

        # Dimensions

        # Time is an implied dimension in the netcdf file for most variables
        # time/datetime is necessary if an alternative time
        # dimenison does not appear in even one variable
        ndays_time_vars = [
            "soltab_potsw",
            "soltab_horad_potsw",
            "soltab_sunhrs",
        ]

        for var_name in variables:
            if var_name in ndays_time_vars:
                continue
            # None for the len argument gives an unlimited dim
            self.dataset.createDimension("time", None)
            self.datetime = self.dataset.createVariable(
                "datetime", "f4", ("time",)
            )
            self.datetime.units = time_units
            break

        # similarly, if alternative time dimenions exist.. define them
        for var_name in variables:
            if var_name not in ndays_time_vars:
                continue
            self.dataset.createDimension("ndays", 366)
            self.doy = self.dataset.createVariable("doy", "i4", ("ndays",))
            self.doy.units = "Day of year"
            break

        if nhru_coordinate:
            self.dataset.createDimension("nhm_id", self.nhrus)
        if nsegment_coordinate:
            self.dataset.createDimension("nhm_seg", self.nsegments)

        if nhru_coordinate:
            self.hruid = self.dataset.createVariable(
                "nhm_id", "i4", ("nhm_id")
            )
            self.hruid[:] = np.array(self.hru_ids, dtype=int)
        if nsegment_coordinate:
            self.segid = self.dataset.createVariable(
                "nhm_seg", "i4", ("nhm_seg")
            )
            self.segid[:] = np.array(self.hru_segments, dtype=int)

        self.variables = {}
        for var_name in variables:
            variabletype = meta_netcdf_type(var_meta[var_name])
            if "nsegment" in variable_dimensions[var_name]:
                spatial_coordinate = "nhm_seg"
            else:
                spatial_coordinate = "nhm_id"

            if var_name in ndays_time_vars:
                time_dim = "ndays"
            else:
                time_dim = "time"

            self.variables[var_name] = self.dataset.createVariable(
                var_name,
                variabletype,
                (time_dim, spatial_coordinate),
                fill_value=nc4.default_fillvals[variabletype],
                zlib=zlib,
                complevel=complevel,
                chunksizes=tuple(chunk_sizes.values()),
            )
            for key, val in var_meta[var_name].items():
                if isinstance(val, dict):
                    continue
                self.variables[var_name].setncattr(key, val)

        return

    def __del__(self):
        self.close()
        return

    def close(self):
        if self.dataset.isopen():
            self.dataset.close()
            return

    def add_simulation_time(self, itime_step: int, simulation_time: float):
        self.datetime[itime_step] = nc4.date2num(
            simulation_time, self.datetime.units
        )
        return

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
        return

    def add_all_data(
        self,
        name: str,
        data: np.ndarray,
        time_data: np.ndarray,
        time_coord: str = "datetime",
    ) -> None:
        """Add data to a NetCDF variable

        Args:
            name:

        Returns:

        """
        self[time_coord][:] = time_data
        if name not in self.variables.keys():
            raise KeyError(f"{name} not a valid variable name")
        self.variables[name][:, :] = data[:, :]

        return


class NetCdfCompare(Accessor):
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
            if not np.allclose(compare_arr, base_arr, atol=atol):
                success = False
                diff = np.abs(compare_arr - base_arr)
                failures[variable] = (
                    np.max(diff),
                    atol,
                    np.argmax(np.max(diff, axis=0)),
                )
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
