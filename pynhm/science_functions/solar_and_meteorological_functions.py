from numpy import arccos, sin, cos, tan, pi


def daylight_hours( omega_s ):
    """
    Calculate the number of daylight hours at a location.

    omega_s   Sunset hour angle in Radians.

    Implementation follows equation 34, Allen and others (1998).
    """
    return 24.0 / pi * omega_s



def sunrise_sunset_angle__omega_s(latitude, delta):
    """
    Return the sunset angle in radians, given the latitude and solar declination, in radians.

    latitude        Latitude, in Radians
    delta           Solar declination, in Radians

    Implementation follows equation 25, Allen and others (1998).
    """
    omega_s = arccos( - tan(latitude) * tan(delta) )
    return omega_s



def day_angle__Gamma(day_of_year, number_of_days_in_year):
    """
    Return the day angle in Radians.  

    day _of_year             Integer day of the year (January 1 = 1)
    number_of_days_in_year   Number of days in the current year

    Implementation follows equation 1.2.2 in Iqbal (1983)
    """
    day_angle = 2.0 * pi * ( day_of_year - 1.0 ) / number_of_days_in_year
    return day_angle



def solar_declination__delta(day_of_year, number_of_days_in_year):
    """
    Return the solar declination for a given day of the year in Radians.

    day _of_year             Integer day of the year (January 1 = 1)
    number_of_days_in_year   Number of days in the current year

    Implementation follows equation 1.3.1 in Iqbal (1983).
    Iqbal (1983) reports maximum error of 0.0006 radians; if the last two terms are omitted,
      the reported accuracy drops to 0.0035 radians.
    """
    Gamma = day_angle__Gamma( day_of_year, number_of_days_in_year )

    delta =   0.006918                                        \
             - 0.399912 * cos( Gamma )                        \
             + 0.070257 * sin( Gamma )                        \
             - 0.006758 * cos( 2.0 * Gamma )                  \
             + 0.000907 * sin( 2.0 * Gamma )                  \
             - 0.002697 * cos( 3.0 * Gamma )                  \
             + 0.00148  * sin( 3.0 * Gamma )

    return delta


def relative_earth_sun_distance__D_r(day_of_year, number_of_days_in_year):
    """
    Return the inverse relative Earth-Sun distance (unitless) for a given day of the year.

    day_of_year              Integer day of the year (January 1 = 1)
    number_of_days_in_year   Number of days in the current year

    Implementation follows equation 23, Allen and others (1998). 
    See also Equation 1.2.3 in Iqbal, Muhammad (1983-09-28).
    """
    d_r = 1.0 + 0.033                                      \
            * cos( 2.0 * pi * day_of_year )                \
            / number_of_days_in_year      

    return d_r


def extraterrestrial_radiation__Ra(latitude, delta, omega_s, d_r):
    """
    Return the extraterrestrial solar radiation for a point above earth, in MJ / m^2 / day

    latitude        latitude in Radians
    delta           solar declination in Radians
    omega_s         sunset hour angle in Radians
    d_r             Inverse relative sun-earth distance (unitless)

    Implemented as equation 21, Allen and others (1998).
    """
    Gsc = 0.0820   # MJ / m^2 / min

    part_a = omega_s * sin( latitude ) * sin( delta )
    part_b = cos( latitude ) * cos( delta ) * sin( omega_s )

    Ra = 24.0 * 60.0 * Gsc * d_r * ( part_a + part_b ) / pi
    return Ra



def equivalent_evaporation(radiation_energy):
    """
    Returns the equivalent depth of water (in millimeters) that would be evaporated
    per day for a given amount of solar radiation (expressed in MJ/m^2-day).

    radiation_energy        Solar radiation energy expressed in MJ per sq. meter per day

    Implementation follows equation 20, Allen and others (1998).
    """

    return radiation_energy * 0.408


def references():
    """    
    Allen, R.G., Pereira, L.S., Raes, D., and Smith, M., 1998, Crop evapotranspiration-Guidelines for 
        computing crop water requirements-FAO Irrigation and drainage paper 56: Food and Agriculture Organization
        of the United Nations, Rome, 333 p.

    Iqbal, M., 1983, An Introduction To Solar Radiation: Academic Press Canada, Ontario, Canada, 617 p.

    """
    pass