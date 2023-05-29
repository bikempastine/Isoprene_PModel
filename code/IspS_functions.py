# import required package
import math

# Define the isoprene synthase function
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

# Define the full SS function as it will be called in the Energetic Status Model
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
