"""
Save CO2 and elevation into the Pmodel_data file
"""

# Name the file paths
INPUT_FILE_MAIN = "$HOME/full_averaged/Pmodel_4_of_6_data.nc"
INPUT_FILE_CO2 = "$HOME/co2_data/2020monthly_in_situ_co2_mlo_clean.csv"
INPUT_FILE_ELEV = "$HOME/elevation/ASurf_WFDE5_CRU_v2.0.nc"
OUTPUT_FILE_PATH = "$HOME/full_averaged/Pmodel_data.nc"
START_YEAR = 2005
END_YEAR = 2016

# Load in the required packages
import math
import netCDF4
import numpy as np
import pandas as pd

# Load the Pmodel_data file with all the rest of the variables
ds = netCDF4.Dataset(INPUT_FILE_MAIN,  "r")
ds.set_auto_mask(True)

# Load the elevation data
el = netCDF4.Dataset(INPUT_FILE_ELEV)
el.set_auto_mask(True)

# Load the CO2 data
df = pd.read_csv(INPUT_FILE_CO2) #load in cleaned data: cleaning is removing info at top, and putting colnames in one cell
co2_df = pd.DataFrame({'year': df.iloc[:,0], 'month': df.iloc[:,1], 'ppm':df.iloc[:,4]}) # grab the columns needed, and rename the columns.
filtered_co2_df = co2_df[(co2_df['year'] >= START_YEAR) & (co2_df['year'] <= END_YEAR)] #filter between referace years
co2_values = filtered_co2_df.iloc[:, 2] # Extracting the values of interest from filtered_co2_np

# Extract an example variable from the main netCDF file that will be used to get the correct shape
temp = ds["Tair"][:] 
ds.close()

# Prepare the elevation data ro be added to the main NetCDF
elev_one = el["ASurf"][:] #elevation only has one time point
elev = np.tile(elev_one, (temp.shape[0],1,1)) #tile across all time points
el.close()

# Prepare the CO2 data by creating a new array with the same shape as ppfd, filled with co2_values
co2 = np.tile(co2_values[:, np.newaxis, np.newaxis], (1, temp.shape[1], temp.shape[2]))

# Save these variables to a new netcdf file along with the rest of the variables
with netCDF4.Dataset(INPUT_FILE_MAIN, "r") as input_dataset: # Open the input NetCDF file in read mode
    # Create a new NetCDF file in write mode
    with netCDF4.Dataset(OUTPUT_FILE_PATH, "w") as output_dataset:

        # Copy all dimensions from the input dataset to the output dataset
        for dim_name, dim in input_dataset.dimensions.items():
            output_dataset.createDimension(dim_name, len(dim))

        # Copy all variables and their data from the input dataset to the output dataset
        for var_name, var in input_dataset.variables.items():
            output_dataset.createVariable(var_name, 'float64', input_dataset[var_name].dimensions)
            output_dataset.variables[var_name][:] = input_dataset[var_name][:]

        # Create a new variable in the output dataset for elev
        output_dataset.createVariable("Elev", elev.dtype, input_dataset['Tair'].dimensions)
        output_dataset.variables["Elev"][:] = elev
        
        # Create a new variable in the output dataset for co2
        output_dataset.createVariable("CO2", co2.dtype, input_dataset['Tair'].dimensions)
        output_dataset.variables["CO2"][:] = co2

        # Save the changes and close the output dataset
        output_dataset.sync()


