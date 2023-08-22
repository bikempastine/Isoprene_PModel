
"""
Module containing the J, Jmax, and f(T) functions which take inputs from the Pmodel 
and returns parameters required to run the Energetic Status Isoprene model.

Functions:
    find_J(a_j,ci,gammastar) : light limited electron flux
    find_Jv(vcmax, ci, gammastar, kmm): Rubisco limited electron flux
    SS(T) : Estimates Isoprene synthase activity
    IspS_function(T_celsius, T_standard) : Function referred to as f(T) in the Energetic Status Model

Dependencies: 
    import math as math
"""

# import required package
import math

def find_J(a_j, ci, gammastar):
    """
    Function finds J (light limited electron flux)
    Equation comes from Morfopolous PhD Chapter 2, Equation 2.16 for Aj

    Args:
        a_j (float): light limmited gross assimilation (umol m-2 s-1)
        ci (float): internal CO2 concentration (Pa)
        gammastar (float): compensation point in the absance of dark resperation (Pa)
    
    Returns:
        J (float): light limited electron flux (electron flux used for C and O2 assimilation) (umols m-2 s-1)
    """
    J = 4 * a_j * ((ci + 2*gammastar) / (ci - gammastar))
    return J

def find_Jv(vcmax, ci, gammastar, kmm):
    """
    Function finds Jv (electron flux for Rubisco limmited assimilation)
    Formula comes from Morfopolous et al., 2014 (new physiologist)

    Args: #need to get the units
        vcmax (float): maximum carboxylation capacity (umols m-2 s-1)
        ci (float): light limmited gross assimilation (ubar)
        gammastar (float): compensation point in the absance of dark resperation (ubar)
        kmm (float): Michalis Menten coefficient for C, when there is competition for Rubisco by O2 and C (ubar)
    
    Returns:
        Jv (float): Rubisco limited electron flux (umols m-2 s-1)
    """    
   
    #input comes in in Pa , assume Vcmax is correct for now
    # convert Pa to ubar
    ci_ubar = ci *10
    gammastar_ubar = gammastar *10
    kmm_ubar = kmm *10


    Jv = 4 * vcmax * ((ci_ubar + 2*gammastar_ubar)/(ci_ubar + kmm_ubar)) 
    return Jv


def SS(T):
    """
    Function referred to as SS in the Energetic Status Model
    Estimates Isoprene synthase activity as introduced by Niinemets et al (1999)

    Args:
        T_celsius (float): The temperature input in Celsius

    Returns:
        float: SS: isoprene synthase activity (umol isoprene (g isoprene synthase)-1 s-1)
    """
     
    Tk = T + 273.15  # Convert Celsius to Kelvin (absolute temperature)
    R = 8.314  # Gas constant (Jmol-1K-1)
    C = 35.478  # Unitless scaling constant
    dHa = 83.129 * 1000  # Activation energy (Jmol-1)
    dHd = 284.6 * 1000  # Deactivation energy (Jmol-1)
    dS = 0.8875 * 1000  # Entropy term (JK-1mol-1)

    isoprene_synthase = math.exp(C - (dHa / (R * Tk))) / (1 + math.exp((dS * Tk - dHd) / (R * Tk)))
    return isoprene_synthase


def IspS_function(T_celsius, T_standard):
    """
    Function referred to as f(T) in the Energetic Status Model
    Estimates Isoprene synthase activity as introduced by Niinemets et al (1999)

    Args:
        T_celsius (float): The temperature input in Celsius
        T_standard (float): The standard temperature in Celsius to evaluate SSS at in Celcius

    Returns:
        float: SS/SSS or f(T) (unitless)
    """

    fT = SS(T_celsius) / SS(T_standard)  # Unitless

    return fT