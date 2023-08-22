"""
Process MEGAN isoprene monthly nc files, and concatenate them together to produce a single file with the entire time series
"""

import numpy as np
import netCDF4
import subprocess

# Define the range of years
start_year = 2005
end_year = 2014

# Name the output file
output_file = "clean_data/MEGAN_2005_2014.nc"

def process_megan_data(start_year, end_year):
    """
    Process MEGAN data for a range of years and save as netCDF files.
    
    Args:
        start_year (int): Starting year.
        end_year (int): Ending year.
    """
    
    # Loop over the years
    for i in range(start_year, end_year + 1):
        # Open MEGAN data for the current year
        megan_file = f"../data/MEGAN/ISOP_MEGAN-ECMWF-v2-MONTHLY_{i}.nc"
        MEGAN_cdf = netCDF4.Dataset(megan_file)
        MEGAN_cdf.set_auto_mask(False)

        # Read MEGAN data
        megan_lon_flip = MEGAN_cdf['Longitude'][:]
        megan_flip = MEGAN_cdf['Flux'][:]

        # Perform necessary operations on the data

        # Create a new netCDF file for the current year
        new_file = f"Megan_iso_{i}.nc"
        new_cdf = netCDF4.Dataset(new_file, 'w')

        # Create dimensions in the new netCDF file
        lat_dim = new_cdf.createDimension('lat', megan.shape[1])
        lon_dim = new_cdf.createDimension('lon', megan.shape[2])
        time_dim = new_cdf.createDimension('time', None)

        # Create variables in the new netCDF file
        lat_var = new_cdf.createVariable('lat', 'f4', ('lat',))
        lon_var = new_cdf.createVariable('lon', 'f4', ('lon',))
        time_var = new_cdf.createVariable('time', 'f8', ('time',))
        megan_iso = new_cdf.createVariable('megan_iso', 'f4', ('time', 'lat', 'lon',))

        # Assign values to variables in the new netCDF file
        lat_var[:] = MEGAN_cdf.variables['Latitude'][:] * -1
        lon_var[:] = megan_lon
        megan_iso[:] = np.flip(megan, axis=1)

        # Set attributes for the time variable
        time_var.units = 'months since 2005-01-01 00:00:00'
        time_var.standard_name = 'time'
        time_var.calendar = '360_day'

        # Assign the time values to the time variable
        time_var[:] = range(0 + (12*(i-start_year)), 12 + (12*(i-start_year)))

        # Close the current netCDF file
        new_cdf.close()
        MEGAN_cdf.close()

def merge_megan_data(start_year, end_year, output_file):
    """
    Merge MEGAN data for a range of years into a single netCDF file.
    
    Args:
        start_year (int): Starting year.
        end_year (int): Ending year.
        output_file (str): Output netCDF file name.
    """
    
    # Define the Bash command for concatenating the files
    bash_cat = f"cdo cat {' '.join([f'Megan_iso_{i}.nc' for i in range(start_year, end_year + 1)])} {output_file}"

    # Execute the Bash command
    subprocess.run(bash_cat, shell=True)
    
    # Define the Bash command for cleaning up the temporary files
    bash_clean = 'rm Megan_iso_*.nc'

    # Execute the Bash command
    subprocess.run(bash_clean, shell=True)


# Process and save MEGAN data
process_megan_data(start_year, end_year)

# Merge MEGAN data into a single file
merge_megan_data(start_year, end_year, output_file)
