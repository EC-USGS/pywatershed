import numpy as np


def calculate_cn_S_inches(curve_number):
    """
    Return the curve number storage (S) term, in inches. Equation 2-4, Cronshey and others (1986).
    
    """

    S_inches = ( 1000.0 / curve_number ) - 10.0
    return S_inches



def calculate_cn_S_millimeters(curve_number):
    """
    Return the curve number storage (S) term, in millimeters. Equation 2-4, Cronshey and others (1986),
    with constants multiplied by 25.4 (mm to inches).
    
    """

    S_mm = ( 25400.0 / curve_number ) - 254.0
    return S_mm



def calculate_cn_runoff(inflow, storage_S, initial_abstraction_Ia = 0.05):
    """
    Return the runoff value given the inflow (precip), storage, and initial abstraction. 
    Equation 2-3, Cronshey and others (1986).
    """
    Ia = initial_abstraction_Ia
    runoff = np.where(inflow > Ia,
                            np.float_power(inflow - Ia * storage_S, 2.0) / (inflow + (1.0 - Ia) * storage_S),
                            0.0
                     )

    return runoff



def calculate_cn_alternative_S_0_05(storage_S):
    """
    Return the curve number storage term, assuming that the initial abstraction is 0.05, rather than 0.2.
    Equation 8, Woodward and others (2003).
    """
    return 1.33 * (storage_S^1.15)



def calculate_cn_arc2_to_arc1(curve_number_arc2):
    """
    Return a curve number corresponding to antecedant runoff condition 1, given a
    curve number corresponding to antecedant runoff condition 2.
    Implemented as equation 3.145 of "SCS Curve Number Methodology", Mishra and Singh (2003),
    and as equation 15 in Ponce and Hawkins (1996).
    
    Resulting curve numbers are clipped to the range 30-100.
    """
    return np.clip((curve_number_arc2 / (2.281 - 0.01281 * curve_number_arc2 )),      
                30.0,
                100.0
               )


def calculate_cn_arc2_to_arc3(curve_number_arc2):
    """
    Return a curve number corresponding to antecedant runoff condition 3, given a
    curve number corresponding to antecedant runoff condition 2.
    Implemented as equation 3.146 of "SCS Curve Number Methodology", Mishra and Singh (2003), 
    and as equation 16 in Ponce and Hawkins (1996).

    Resulting curve numbers are clipped to the range 30-100.
    """
    return np.clip((curve_number_arc2 / (0.427 - 0.00573 * curve_number_arc2 )),      
                30.0,
                100.0
               )

def adjust_curve_number(curve_number, inflow_5_day_sum, is_growing_season=False, cfgi=0, cfgi_ll=55, cfgi_ul=85):

    ARC_DRY_GROWING = 1.40
    ARC_DRY_DORMANT = 0.50
    ARC_WET_GROWING = 2.10
    ARC_WET_DORMANT = 1.10

    if( cfgi > cfgi_ll ):

        p_er = calculate_probability_of_enhanced_runoff(cfgi, cfgi_ll, cfgi_ul)

        curve_number_adj = curve_number * (1. - p_er) + curve_number_arc3 * p_er

    elif( is_growing_season ):

        curve_number_adj = np.where(inflow_5_day_sum > ARC_WET_GROWING, curve_number_arc3,
                             np.where(inflow_5_day_sum < ARC_DRY_GROWING, curve_number_arc1, curve_number))

    else:

        curve_number_adj = np.where(inflow_5_day_sum > ARC_WET_DORMANT, curve_number_arc3,
                             np.where(inflow_5_day_sum < ARC_DRY_DORMANT, curve_number_arc1, curve_number))

    return curve_number_adj



def cn_references():
    """
    Cronshey, R., McCuen, R., Miller, N., Rawls, W., Robbins, S., and Woodward, D., 1986, Urban Hydrology
        for Small Watersheds - Technical release 55: US Dept. of Agriculture, Soil Conservation Service, 
        Engineering Division, accessed at http://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/16/stelprdb1044171.pdf.

    Mishra, S.K., and Singh, V.P., 2003, Soil Conservation Service Curve Number (SCS-CN) Methodology: Water Science
         and Technology Library, Springer Netherlands, Dordrecht, 534 p.

    Ponce, V.M., and Hawkins, R.H., 1996, Runoff Curve Number: Has It Reached Maturity? Journal of Hydrologic
        Engineering, v. 1, no. 1, p. 11-19.

    Woodward, D.E., Hawkins, R.H., Jiang, R., Hjelmfelt, J., Van Mullem, J.A., and Quan, Q.D., 2003, 
        Runoff Curve Number Method: Examination of the Initial Abstraction Ratio, 
        in World Water and Environmental Resources Congress 2003,
        American Society of Civil Engineers, Philadelphia, Pennsylvania, p. 1-10.

    """
    pass

# 