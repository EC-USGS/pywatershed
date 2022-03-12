# JLM: is this stuff in parker's XML?
# https://github.com/paknorton/pyPRMS/blob/f41a7911af0b34c575005ed5d133c36111ccccf3/pyPRMS/xml/variables.xml
# in dimensions... sort of, not exactly, there are sort of variations...
# JLM: specify fill values? or adopt global fill values somewhere. this probably depends on what the model uses
# the netcdf output is currently taking what the defaults are for netCDF4

# This might live in a kindof constantsmodule or something
# some of these values are also defined in xarray or netcdf4 i think
fill_value_f4 = 9.96921e36

cbh_metadata = {
    "datetime": {
        "type": "f4",
        "long_name": "time",
        "standard_name": "time",
        "calendar": "standard",  # Depends, revisit
        "units": "days since 1979-01-01 00:00:00",  # Depends, may not be correct
    },
    "hru_ind": {
        "type": "i4",
        "long_name": "Hydrologic Response Unit (HRU) index",
        "cf_role": "timeseries_id",
    },
    "nhm_id": {
        "type": "i4",
        "long_name": "NHM Hydrologic Response Unit (HRU) ID",
        "cf_role": "timeseries_id",
    },
    "tmax": {
        "type": "f4",
        "_FillValue": fill_value_f4,
        "long_name": "Maximum daily air temperature",
        "units": "degree_fahrenheit",
        "standard_name": "maximum_daily_air_temperature",
    },
    "tmax_adj": {
        "type": "f4",
        "_FillValue": fill_value_f4,
        "long_name": "Maximum daily air temperature from parameter adjustments",
        "units": "degree_fahrenheit",
        "standard_name": "maximum_daily_air_temperature",
    },
    "tmin": {
        "type": "f4",
        "_FillValue": fill_value_f4,
        "long_name": "Minimum daily air temperature",
        "units": "degree_fahrenheit",
        "standard_name": "minimum_daily_air_temperature",
    },
    "tmin_adj": {
        "type": "f4",
        "_FillValue": fill_value_f4,
        "long_name": "Minimum daily air temperature from parameter adjustments",
        "units": "degree_fahrenheit",
        "standard_name": "minimum_daily_air_temperature",
    },
    "prcp": {
        "type": "f4",
        "_FillValue": fill_value_f4,
        "long_name": "daily total precipitation",
        "units": "in",
        "standard_name": "daily_total_precipitation",
    },
    "prcp_adj": {
        "type": "f4",
        "_FillValue": fill_value_f4,
        "long_name": "daily total precipitation from parameter adjustments",
        "units": "in",
        "standard_name": "daily_total_precipitation",
    },
    "rainfall_adj": {
        "type": "f4",
        "_FillValue": fill_value_f4,
        "long_name": "daily rainfall from parameter adjustments",
        "units": "in",
        "standard_name": "daily_rainfall",
    },
    "rhavg": {
        "type": "f4",
        "_FillValue": fill_value_f4,
        "long_name": "Daily mean relative humidity",
        "units": "percent",
        "standard_name": "rhavg",
    },
    "snowfall_adj": {
        "type": "f4",
        "_FillValue": fill_value_f4,
        "long_name": "daily snowfall from parameter adjustments",
        "units": "in",
        "standard_name": "daily_snowfall",
    },
}
