from numpy import nanmax, float_power
from ..constants import DEGREES_TO_RADIANS, RADIANS_TO_DEGREES
from .solar_and_meteorological_functions import relative_earth_sun_distance__D_r, \
                                                solar_declination__delta,         \
                                                sunrise_sunset_angle__omega_s,    \
                                                extraterrestrial_radiation__Ra,   \
                                                equivalent_evaporation

def et0_hargreaves_samani(extraterrestrial_radiation__Ra, air_temp_min, air_temp_max, air_temp_mean, 
                            et_slope=0.0023, et_constant=17.8, et_exponent=0.5):
    """
    Return the daily reference evapotranspriation in millimeters given extraterrestrial
    radiation (expressed as the number of millimeters of water that could be evaporated by
    the applied radiation) and the minimum and maximum daily air temperatures in Celsius degrees.

    Implemented as equation 4 in Hagreaves and Samani (1985) and as equation 50 in Allen and others (1998).

    """
    air_temp_delta = air_temp_max - air_temp_min

    et0 = nanmax((0.0,
                  et_slope * extraterrestrial_radiation__Ra * (air_temp_mean + et_constant)          \
                       * float_power(air_temp_delta,et_exponent)) 
             )

    return et0


def calculate_et0_hargreaves_samani(day_of_year, number_of_days_in_year, latitude,
                                        air_temp_min, air_temp_max, air_temp_mean):
    """
    Return the calculated reference evapotranspiration in millimeters per day, 
    given the day number, latitude, and min, max, and mean air temperatures in Celsius.

    day_of_year                 Day number of the current solar year.
    number_of_days_in_year      Number of days in the solar year.
    latitude                    Latitude (in degrees) at which calculation is to be made.
    air_temp_min                Minimum daily air temperature, in degrees Celsius
    air_temp_max                Maximum daily air temperature, in degrees Celsius
    air_temp_mean               Mean daily air temperature, in degrees Celsius
    """    
    latitude_radians = latitude * DEGREES_TO_RADIANS

    d_r     = relative_earth_sun_distance__D_r(day_of_year, number_of_days_in_year)
    delta   = solar_declination__delta(day_of_year, number_of_days_in_year)
    omega_s = sunrise_sunset_angle__omega_s(latitude_radians, delta)

    Ra = equivalent_evaporation(extraterrestrial_radiation__Ra(latitude_radians, delta, omega_s, d_r))
    ref_ET = et0_hargreaves_samani(Ra, air_temp_min, air_temp_max, air_temp_mean)
    print(f'ref_ET: {ref_ET}')
    return ref_ET



def references():
    """
    Allen, R.G., Pereira, L.S., Raes, D., and Smith, M., 1998, Crop evapotranspiration-Guidelines
        for computing crop water requirements-FAO Irrigation and drainage paper 56: Food and Agriculture
        Organization of the United Nations, Rome, 333 p.

    Hargreaves, G.H., and Samani, Z.A., 1985, Reference Crop Evapotranspiration from Ambient Air
        Temperature: American Society of Agricultural Engineers 85–2517, 12 p.,
        accessed December 21, 2015, at http://www.zohrabsamani.com/papers/Hargreaves_Samani_85.pdf.
    """