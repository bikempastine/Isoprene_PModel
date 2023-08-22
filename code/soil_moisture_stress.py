"""
Function which calculates the soil moisture stress parameter beta using the pmodel and preprepared data from SPLASH 

Functions

Dependencies: 
    import math as math
    import netCDF4
    import numpy as np
    from pyrealm import pmodel

"""

import math
import netCDF4
import numpy as np
from pyrealm import pmodel

def calculate_beta(SOILM_PATH, PET_PATH, AET_PATH, START_YEAR, END_YEAR, MEAN_ALPHA_MODE, BETA_SAVE_PATH):
    """
    Find beta (empirical soil moisture stress) for each gridcell 0.5/0.5 degree grid 

    Inputs: 
     - SOILM_PATH: path to SPLASH_v2.0_sm_lim_monthly_1901-2020.nc
     - PET_PATH: path to SPLASH_v2.0_pet_monthly_1901-2020.nc
     - AET_PATH: path to SPLASH_v2.0_aet_monthly_1901-2020

     - JAN_START: the year your time series starts
     - JAN_END: the year your time series ends (inclusive)

     - MEAN_ALPA_MODE: The method to calculate the mean alpha
        - 'year': calculate a different mean alpha for each year as described by the Stocker et al., 2020
        - 'total': calculate a single mean alpha for the whole time series as suppested by Dr. David Sandoval Calle

    Outputs:
        - Saves Beta as a .npy for each gridcell and for each month in the timeseries (upsidedown and not masked)

    """

    # Get the month number
    def find_jan_no(year):
        """Finds the month number since January 1901 of the January of the year provided"""
        month = 12 * (year - 1901)
        return month

    def set_nan(data):
        """Sets fill values to nan"""
        masked = np.where((data < 0), np.nan, data)
        return masked

    # Get the start and end month number
    start = find_jan_no(START_YEAR)
    end = find_jan_no(END_YEAR + 1)

    # Soilm
    soilm_ncdf = netCDF4.Dataset(SOILM_PATH)  # fraction of available water content
    soilm = soilm_ncdf['sm_lim'][start:end]  # subset soilm to the correct months
    soilm_ncdf.close()

    soilm_clean = set_nan(soilm)  # set the fill values to NAN

    # AET/PET
    pet_netcdf = netCDF4.Dataset(PET_PATH)  # potential evapotranspiration
    aet_netcdf = netCDF4.Dataset(AET_PATH)  # actual evapotranspiration

    pet = pet_netcdf['pet'][start:end]
    aet = aet_netcdf['aet'][start:end]

    pet_netcdf.close()
    aet_netcdf.close()

    pet_clean = set_nan(pet)
    aet_clean = set_nan(aet)

    if MEAN_ALPHA_MODE == 'total':
        # Average over time
        average_pet = np.nanmean(pet_clean, axis=0)
        average_aet = np.nanmean(aet_clean, axis=0)

        # AET/PET = mean alpha
        mean_alpha_total = np.divide(average_aet, average_pet, out=np.zeros_like(average_aet), where=(average_pet != 0))
        mean_alpha = np.tile(mean_alpha_total, (soilm_clean.shape[0], 1, 1)) # need to get this into the same shape as soilm to run the pmodel

    elif MEAN_ALPHA_MODE == 'year':
        # AET/PET = mean alpha
        mean_alpha = np.divide(aet_clean, pet_clean, out=np.zeros_like(aet_clean), where=(pet_clean != 0))
    else:
        print("Please input a value for MEAN_ALPHA_MODE that is either 'total' or 'year' depending on the output you need")
        return

    # Get Beta
    beta = pmodel.calc_soilmstress(soilm_clean, mean_alpha)

    # Save beta to the file spesified 
    np.save(BETA_SAVE_PATH, beta)

    # print a statement to the console
    
    return (
        print("BETA found from:", START_YEAR, 'to', END_YEAR, '\nResults saved to the path:', 
        BETA_SAVE_PATH, 
        "\nMean alpha was found by the method:", MEAN_ALPHA_MODE)
    )


if __name__ == "__main__":
    # Input file paths
    SOILM_PATH = "../data/Splash/SPLASH_v2.0_sm_lim_monthly_1901-2020.nc"
    PET_PATH = "../data/Splash/SPLASH_v2.0_pet_monthly_1901-2020.nc"
    AET_PATH = "../data/Splash/SPLASH_v2.0_aet_monthly_1901-2020.nc"
    BETA_SAVE_PATH = "../results/Soil_Moisture_Stress_Beta.npy"

    # Start time inputs
    START_YEAR = 2005
    END_YEAR = 2016

    # mean_alpha mode
    MEAN_ALPHA_MODE = 'year' #'total' 

    # Calculate beta
    calculate_beta(SOILM_PATH, PET_PATH, AET_PATH, START_YEAR, END_YEAR, MEAN_ALPHA_MODE, BETA_SAVE_PATH)

