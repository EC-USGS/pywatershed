import datetime as dt
import pathlib as pl
from math import ceil
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

# JLM TODO: the implied time dimension seems like a bad idea, it should be
#    an argument.


class NetCdfRead(Accessor):
    """NetCDF file reader (for input/forcing data)

    Args:
      name: the netcdf file path to open
      start_time: optional np.datetime64 which is the start of the simulation
      end_time: optional np.datetime64 which is the end of the simulation
      nc_read_vars: a subset of available variables to read from the file
      load_n_times: optional integer for the length of time to load into memory
        from file see load_n_time_batches for more details. default is None. A
        value of -1 gives the same result as load_n_time_batches = 1.
      load_n_time_batches: optional integer for the number to batches (time
        partitions) of the input data to use. Default value is 1 whic loads all
        the input data on the first time advance. This may not work well for
        domains large in space and/or time as it will require large amounts of
        memory. These load options exist to optimize IO patterns against memory
        usage: specify EXACTLY ONE of load_n_times or load_ntime_batches,
        whichever is more convenient for you. If both _load_times and
        _load_n_time_batches are None, then no time batching is used. This has
        proven an inefficient pattern. The time batching is not implemented for
        DOY (cyclic) variables, only for variables with time dimension "time".
    """

    def __init__(
        self,
        name: fileish,
        start_time: np.datetime64 = None,
        end_time: np.datetime64 = None,
        nc_read_vars: list = None,
        load_n_times: int = None,
        load_n_time_batches: int = 1,
    ) -> "NetCdfRead":
        self.name = "NetCdfRead"
        self._nc_file = name
        self._nc_read_vars = nc_read_vars
        self._start_time = start_time
        self._end_time = end_time

        if (load_n_times is not None) and (load_n_time_batches is not None):
            msg = "Can only specify one of load_n_times or load_ntime_batches"
            raise ValueError(msg)

        self._load_n_times = load_n_times
        self._load_n_time_batches = load_n_time_batches

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

        if "time" in self.dataset.variables:
            self._time = (
                nc4.num2date(
                    self.dataset.variables["time"][:],
                    units=self.dataset.variables["time"].units,
                    calendar="standard",
                    only_use_cftime_datetimes=False,
                )
                .filled()
                .astype("datetime64[s]")
                # JLM: the global time type as in cbh_utils, define somewhere
            )

            if self._start_time is None:
                self._start_index = 0
            else:
                wh_start = np.where(self._time == self._start_time)
                self._start_index = wh_start[0][0]

            if self._end_time is None:
                self._end_index = self._time.shape[0] - 1
            else:
                wh_end = np.where(self._time == self._end_time)
                self._end_index = wh_end[0][0]

            self._time = self._time[self._start_index : (self._end_index + 1)]
            self._ntimes = self._end_index - self._start_index + 1

            # time batching
            # data are actually loaded in get_data
            if self._load_n_time_batches is not None:
                # use ceil because we want exactly the requested # of batches
                self._load_n_times = ceil(
                    self._ntimes / self._load_n_time_batches
                )
                self._data_loaded = {}

            elif self._load_n_times is not None:
                # Use ceil to account for the remainder batch
                if self._load_n_times == -1:
                    self._load_n_times = self._ntimes
                self._load_time_batches = ceil(
                    self._ntimes / self._load_n_times
                )
                self._data_loaded = {}

            # Note that if neither _load variables is specified, then no time
            # batching is used

        if "doy" in self.dataset.variables:
            self._doy = self.dataset.variables["doy"][:].data
            self._ntimes = self._doy.shape[0]
            self._start_index = 0
            self._end_index = 365

        spatial_id_names = []
        if "nhm_id" in self.ds_var_list:
            spatial_id_names.append("nhm_id")
        elif "hru_id" in self.ds_var_list:
            spatial_id_names.append("hru_id")
        elif "grand_id" in self.ds_var_list:
            spatial_id_names.append("grand_id")

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
            if name != "time" and name not in spatial_id_names
        ]

        return

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
    def times(
        self,
    ) -> np.ndarray:
        """Get the times in the NetCDF file

        Returns:
            data_times: numpy array of datetimes in the NetCDF file

        """
        if hasattr(self, "_time"):
            return self._time
        elif hasattr(self, "_doy"):
            return self._doy
        else:
            raise KeyError(f"{self.name} has neither _time nor _doy")

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
            variables: list of variable names, excluding the time and
                nhru variables

        """
        return self._variables

    def all_time(self, variable):
        return self.get_data(variable)

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
            return self.dataset[variable][
                self._start_index : (self._end_index + 1), :
            ]

        else:
            if itime_step >= self._ntimes:
                raise ValueError(
                    f"requested time step {itime_step} but only "
                    + f"{self._ntimes} time steps are available."
                )

            if hasattr(self, "_data_loaded"):
                # load when needed, at the beginning of each batch
                batch_index = itime_step % self._load_n_times
                if batch_index == 0:
                    ith_batch = itime_step // self._load_n_times
                    # print(
                    #     f"load batch "
                    #     f"#{ith_batch}/{self._load_n_time_batches-1}: "
                    #     f"{variable}"
                    # )
                    start_ind = self._start_index + (
                        ith_batch * self._load_n_times
                    )
                    end_ind = start_ind + self._load_n_times
                    self._data_loaded[variable] = self.dataset[variable][
                        start_ind:end_ind, :
                    ]

                return self._data_loaded[variable][batch_index, :]

            else:
                # no time batching
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
        if "time" in self.dataset.variables:
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
        global_attrs: dict = {},
        time_units: str = "days since 1970-01-01 00:00:00",
        clobber: bool = True,
        zlib: bool = True,
        complevel: int = 4,
        chunk_sizes: dict = {"time": 1, "hruid": 0},
    ):
        if isinstance(variables, dict):
            group_variables = []
            for group, vars in variables.items():
                for var_name in vars:
                    if group is None:
                        group_variables += [var_name]
                    else:
                        group_variables += [f"{group}/{var_name}"]

            v2 = [
                var_name
                for group, vars in variables.items()
                for var_name in vars
            ]
            variables = v2
        else:
            group_variables = variables

        self.dataset = nc4.Dataset(name, "w", clobber=clobber)
        self.dataset.setncattr("Description", "pywatershed output data")
        for att_key, att_val in global_attrs.items():
            self.dataset.setncattr(att_key, att_val)

        variable_dimensions = {}
        nhru_coordinate = False
        nsegment_coordinate = False
        one_coordinate = False
        nreservoirs_coordinate = False
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
            if "one" in dimension_name:
                one_coordinate = True
            if "nreservoirs" in dimension_name:
                nreservoirs_coordinate = True

        if nhru_coordinate:
            hru_ids = coordinates["nhm_id"]
            if isinstance(hru_ids, (list, tuple)):
                nhrus = len(hru_ids)
            elif isinstance(hru_ids, np.ndarray):
                nhrus = hru_ids.shape[0]
            self.nhrus = nhrus
            self.hru_ids = hru_ids

        if nsegment_coordinate:
            segment_ids = coordinates["nhm_seg"]
            if isinstance(segment_ids, (list, tuple)):
                nsegments = len(segment_ids)
            elif isinstance(segment_ids, np.ndarray):
                nsegments = segment_ids.shape[0]
            self.nsegments = nsegments
            self.segment_ids = segment_ids

        if one_coordinate:
            self.one_ids = coordinates["one"]

        if nreservoirs_coordinate:
            self.nreservoirs = len(coordinates["grand_id"])

        # Dimensions

        # Time is an implied dimension in the netcdf file for most variables
        # time is necessary if an alternative time
        # dimenison does not appear in even one variable
        doy_time_vars = [
            "soltab_potsw",
            "soltab_horad_potsw",
            "soltab_sunhrs",
        ]

        for var_name in variables:
            if var_name in doy_time_vars:
                continue
            # None for the len argument gives an unlimited dim
            self.dataset.createDimension("time", None)
            self.time = self.dataset.createVariable("time", "f4", ("time",))
            self.time.units = time_units
            break

        # similarly, if alternative time dimenions exist.. define them
        for var_name in variables:
            if var_name not in doy_time_vars:
                continue
            self.dataset.createDimension("doy", 366)
            self.doy = self.dataset.createVariable("doy", "i4", ("doy",))
            self.doy.units = "Day of year"
            break

        if nhru_coordinate:
            self.dataset.createDimension("nhm_id", self.nhrus)
        if nsegment_coordinate:
            self.dataset.createDimension("nhm_seg", self.nsegments)
        if one_coordinate:
            self.dataset.createDimension("one", 1)
        if nreservoirs_coordinate:
            self.dataset.createDimension("grand_id", self.nreservoirs)

        if nhru_coordinate:
            self.hruid = self.dataset.createVariable(
                "nhm_id", "i4", ("nhm_id")
            )
            self.hruid[:] = np.array(self.hru_ids, dtype=int)
        if nsegment_coordinate:
            self.segid = self.dataset.createVariable(
                "nhm_seg", "i4", ("nhm_seg")
            )
            self.segid[:] = np.array(self.segment_ids, dtype=int)
        if one_coordinate:
            self.oneid = self.dataset.createVariable("one", "i4", ("one"))
            self.oneid[:] = np.array(self.one_ids, dtype=int)
        if nreservoirs_coordinate:
            self.grandid = self.dataset.createVariable(
                "grand_id", "i4", ("grand_id")
            )
            self.grandid[:] = coordinates["grand_id"]

        self.variables = {}
        for var_name, group_var_name in zip(variables, group_variables):
            variabletype = meta_netcdf_type(var_meta[var_name])
            if len(
                set(["nhru", "ngw", "nssr"]).intersection(
                    set(variable_dimensions[var_name])
                )
            ):
                spatial_coordinate = "nhm_id"
            elif "nsegment" in variable_dimensions[var_name]:
                spatial_coordinate = "nhm_seg"
            elif "one" in variable_dimensions[var_name]:
                spatial_coordinate = "one"
            elif "nreservoirs" in variable_dimensions[var_name]:
                spatial_coordinate = "grand_id"
            else:
                msg = (
                    "Undefined spatial coordinate name in "
                    f"{variable_dimensions[var_name]}"
                )
                raise ValueError(msg)

            if var_name in doy_time_vars:
                time_dim = "doy"
            else:
                time_dim = "time"

            self.variables[var_name] = self.dataset.createVariable(
                group_var_name,
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
        self.time[itime_step] = nc4.date2num(simulation_time, self.time.units)
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
        time_coord: str = "time",
    ) -> None:
        """Add data to a NetCDF variable

        Args:
            name:

        Returns:

        """
        if name not in self.variables.keys():
            raise KeyError(f"{name} not a valid variable name")

        if time_coord == "time":
            start_date = (
                time_data[0].astype(dt.datetime).strftime("%Y-%m-%d %H:%M:%S")
            )
            self[time_coord].units = f"days since {start_date}"
            self[time_coord][:] = nc4.date2num(
                time_data.astype(dt.datetime),
                units=self[time_coord].units,
                calendar="standard",
            )
        else:
            # currently just doy
            self[time_coord][:] = time_data

        self.variables[name][:, :] = data[:, :]

        return
