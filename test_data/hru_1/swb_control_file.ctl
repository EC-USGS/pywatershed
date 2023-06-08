GRID   1    1    1864520       2231406        5000
BASE_PROJECTION_DEFINITION +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23.0 +lon_0=-96.0 +x_0=0.0 +y_0=0.0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs

%% Define methods
-----------------

INTERCEPTION_METHOD              BUCKET
EVAPOTRANSPIRATION_METHOD        HARGREAVES
RUNOFF_METHOD                    CURVE_NUMBER
SOIL_MOISTURE_METHOD             THORNTHWAITE
PRECIPITATION_METHOD             TABULAR
GROWING_DEGREE_DAY_METHOD        BASKERVILLE_EMIN
FOG_METHOD                       NONE
FLOW_ROUTING_METHOD              NONE
IRRIGATION_METHOD                NONE  #NONE turns off all irrigation  #FAO-56 uses that method
ROOTING_DEPTH_METHOD             NONE
CROP_COEFFICIENT_METHOD          NONE  #Used to calculate actual ET with soil moisture
DIRECT_RECHARGE_METHOD           NONE
SOIL_STORAGE_MAX_METHOD          CALCULATED
AVAILABLE_WATER_CONTENT_METHOD   GRIDDED


%% define location, projection, and conversions for weather data
----------------------------------------------------------------


INITIAL_CONTINUOUS_FROZEN_GROUND_INDEX CONSTANT 100.0
UPPER_LIMIT_CFGI 83.
LOWER_LIMIT_CFGI 55.

%% specify location and projection for input GIS grids
------------------------------------------------------

HYDROLOGIC_SOILS_GROUP CONSTANT 2

LAND_USE CONSTANT 22

AVAILABLE_WATER_CONTENT CONSTANT 2.5


%% specify location and names for all lookup tables
---------------------------------------------------

LAND_USE_LOOKUP_TABLE landuse_nlcd.txt
WEATHER_DATA_LOOKUP_TABLE hru_1_weather.txt

%% initial conditions for soil moisture and snow storage amounts
%% may be specified as grids, but using a constant amount and
%% allowing the model to "spin up" for a year is also acceptable.

INITIAL_PERCENT_SOIL_MOISTURE CONSTANT 100.0
INITIAL_SNOW_COVER_STORAGE CONSTANT 2.0


%% OUTPUT CONTROL SECTION:
OUTPUT DISABLE snow_storage
OUTPUT ENABLE tmin tmax
OUTPUT DISABLE crop_et soil_storage delta_soil_storage reference_ET0 
OUTPUT DISABLE runon interception soil_storage delta_soil_storage

OUTPUT ENABLE gross_precipitation 
OUTPUT ENABLE runoff_outside rejected_net_infiltration 
OUTPUT ENABLE runoff actual_et rainfall snowmelt 

DUMP_VARIABLES COORDINATES 1866006 2234241 

%% start and end date may be any valid dates in SWB version 2.0
%% remember to allow for adequate model spin up; running the
%% model for just a month or two will give questionable results

START_DATE 01/01/1979
END_DATE 12/31/2019
