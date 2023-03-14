import numpy as np

TM_SLOPE_TERM_IN      = 0.478769194198665
TM_SLOPE_TERM_MM      = 0.539815721014123
TM_EXPONENT_TERM      = -1.03678439421169

def thornthwaite_mather_soil_moisture_inches(max_soil_moisture, apwl):
    """
    Return the current soil moisture (in inches) given the current soil moisture maximum (in inches),
    as well as the current accumulated potential water loss (APWL). APWL is the running sum of the difference
    between the potential ET and the actual ET. Equation and constants come from an R equation fitting exercise 
    applied to the original soil-mosture retention tables included in Thornthwaite and Mather (1957).

    max_soil_moisture       Maximum moisture content of a soil at field capacity, in inches.
    apwl                    Running sum of daily difference between PET and AET, in inches.
    """

    soil_moisture = np.where( max_soil_moisture > 0.,
        max_soil_moisture * 10**(-TM_SLOPE_TERM_IN*apwl*max_soil_moisture**TM_EXPONENT_TERM),
        0.0
    )
    return(soil_moisture)


def thornthwaite_mather_soil_moisture_millimeters(max_soil_moisture, apwl):
    """
    Return the current soil moisture (in millimeters) given the current soil moisture maximum (in millimeters),
    as well as the current accumulated potential water loss (APWL). APWL is the running sum of the difference
    between the potential ET and the actual ET. Equation and constants come from an R equation fitting exercise 
    applied to the original soil-mosture retention tables included in Thornthwaite and Mather (1957).

    max_soil_moisture       Maximum moisture content of a soil at field capacity, in millimeters.
    apwl                    Running sum of daily difference between PET and AET, in millimeters.
    """

    soil_moisture = np.where( max_soil_moisture > 0.,
        max_soil_moisture * 10**(-TM_SLOPE_TERM_MM*apwl*max_soil_moisture**TM_EXPONENT_TERM),
        0.0
    )
    return(soil_moisture)


def thornthwaite_mather_accumulated_potential_water_loss_inches(max_soil_moisture, soil_moisture):
    """
    Return the accumulated potential water loss (in inches) given the current soil moisture (in inches),
    as well as maximum soil moisture (in inches). APWL is the running sum of the difference
    between the potential ET and the actual ET. Equation and constants come from an R equation fitting exercise 
    applied to the original soil-mosture retention tables included in Thornthwaite and Mather (1957).    
    """
    apwl = np.where(max_soil_moisture > 0.,
        - (np.log(soil_moisture / max_soil_moisture))/(np.log(10)*TM_SLOPE_TERM_IN*max_soil_moisture**TM_EXPONENT_TERM),
        0.0
        )
    return(apwl)



def thornthwaite_mather_accumulated_potential_water_loss_millimeters(max_soil_moisture, soil_moisture):
    """
    Return the accumulated potential water loss (in millimeters) given the current soil moisture (in millimeters),
    as well as maximum soil moisture (in millimeters). APWL is the running sum of the difference
    between the potential ET and the actual ET. Equation and constants come from an R equation fitting exercise 
    applied to the original soil-mosture retention tables included in Thornthwaite and Mather (1957).    
    """
    apwl = np.where(max_soil_moisture > 0.,
        - (np.log(soil_moisture) - np.log(max_soil_moisture))/(np.log(10)*TM_SLOPE_TERM_MM*max_soil_moisture**TM_EXPONENT_TERM),
        0.0
        )
    return(apwl)


def actual_et_references():
    """
    Thornthwaite, C.W., and Mather, J.R., 1957, Instructions and tables for computing potential evapotranspiration
        and the water balance: Publications in Climatology, v. 10, no. 3, p. 1-104.

    
    """
    pass