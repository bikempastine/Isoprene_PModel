"""
The isoprene model, getting and saving critical variables
"""
# Load in the required packages
import math
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import datetime
from dateutil.relativedelta import relativedelta
import sys

## DEFINE VARIABLES
#INDEX = int(sys.argv[1])
INDEX = 8
START_DATE = datetime.datetime(2005, 1, 1)
REFERANCE_TEMP = 30 #reference temperature to run the IspS_function at
INPUT_FILE = "../temporary/Pmodel_data.nc"
OUTPUT_FILE = "../temporary/Iso_results_month_" + str(INDEX)+".nc"

# Load the functions for J, Jv, IspS
import isoprene_functions as Iso_fun

# Load in the pyrealm package
from pyrealm import pmodel
from pyrealm import hygro

# Load the  dataset containing the input variables.
ds = netCDF4.Dataset(INPUT_FILE)
ds.set_auto_mask(True)

# Extract the six variables for all months
temp = ds["Tair"][:] #Celcius
ppfd = ds["SWdown"][:] 
elev = ds["Elev"][:]
vap = ds["vap"][:]
co2 = ds["CO2"][:]
time = ds["time"][:]
fapar_wrong = ds["FPAR"][:] # this needs to be reshaped as there is an issue with the longitudes
ds.close()

np.max(vap)
np.min(vap)
np.max(temp)
np.min(temp)

# Get VPD from VAP
vpd = hygro.convert_vp_to_vpd(vap*0.1 , temp) #vap in hPa, convert to hPA 

# Fix the fapar longitude problem
halfway = np.shape(fapar_wrong)[2]//2
first_half = fapar_wrong[:,:,:halfway]
second_half = fapar_wrong[:,:,halfway:]
fapar = np.concatenate((second_half,first_half), axis=2)
fapar[fapar > 1] = np.nan #fapar is a fraction

# Manipulation of the inputs to fit the assumptions of the pmodel
patm = pmodel.calc_patm(elev) # Convert elevation into  atmospheric pressure
temp[temp < -25] = np.nan # Remove temperature values below -25 degrees
vpd[vpd < 0] = 0 #Clip Vapour pressure defficit so that VPD below zero is set to zero

# Run the PModel for the index defined above
env = pmodel.PModelEnvironment(tc=temp[INDEX,:,:], co2 = co2[INDEX,:,:], patm=patm[INDEX,:,:], vpd=vpd[INDEX,:,:], )
model = pmodel.PModel(env, method_jmaxlim = 'smith19') #run with smith19 to be consistent with J calculation
model.estimate_productivity(fapar[INDEX,:,:], ppfd[INDEX,:,:])

# Get the nessisary variables for J and Jv
gammastar = env.gammastar
ci = model.optchi.ci
vcmax = model.vcmax
kmm = env.kmm
a_j = (model.vcmax / model.optchi.mjoc) * model.optchi.mj #getting ths from the code pyrealm/pyrealm/pmodel/pmodel.py  line 613 and 595

# Get J and Jv
J = Iso_fun.find_J(a_j, ci, gammastar)
Jv = Iso_fun.find_Jv(vcmax, ci, gammastar, kmm)

# Run the IspS function
IspS = np.vectorize(Iso_fun.IspS_function)(temp[INDEX,:,:], REFERANCE_TEMP)

# Produce objects that will be stored
GPP = model.gpp
Jmax = model.jmax
LUE = model.lue
Vcmax = model.vcmax

J_times_IspS_umol = J*IspS


# Save these variables to a new netcdf file along with the rest of the variables
with netCDF4.Dataset(INPUT_FILE, "r") as input_dataset: # Open the input N$
    # Create a new NetCDF file in write mode
    with netCDF4.Dataset(OUTPUT_FILE, "w") as output_dataset:

        dims_to_keep = ['lon', 'lat']  # from the original file we want to get the shape

        # Copy all dimensions from the input dataset to the output dataset
        for dim_name, dim in input_dataset.dimensions.items():
            if dim_name in dims_to_keep:
                output_dataset.createDimension(dim_name, len(dim))
        
        # Create the time dimension with a length of 1
        output_dataset.createDimension("time", 1)

        for var_name, var in input_dataset.variables.items():
            if var_name in dims_to_keep:
                output_dataset.createVariable(var_name, 'float64', var.dimensions)
                output_dataset.variables[var_name][:] = var[:]
                print(var.dimensions)

        time_var = output_dataset.createVariable("time", np.float64, ("time",))

        # Define the units and calendar attributes for the time variable
        time_var.units = "months since 2005-01-01"
        time_var.calendar = "360_day"

        # Get the current date based on the INDEX and set it as the time variable
        time_values = START_DATE + relativedelta(months = INDEX)

        time_numeric = netCDF4.date2num(time_values, units=time_var.units, calendar=time_var.calendar)

        # Assign the time values to the time variable
        time_var[:] = time_numeric

        # Create and assign data to variables in the output file
        output_dataset.createVariable("gpp", GPP.dtype, ("time", "lat", "lon"))[:] = GPP
        output_dataset.createVariable("jmax", Jmax.dtype, ("time", "lat", "lon"))[:] = Jmax
        output_dataset.createVariable("lue", LUE.dtype, ("time", "lat", "lon"))[:] = LUE
        output_dataset.createVariable("vcmax", Jv.dtype, ("time", "lat", "lon"))[:] = Vcmax
        output_dataset.createVariable("J_times_IspS", J_times_IspS_umol.dtype, ("time", "lat", "lon"))[:] = J_times_IspS_umol
        output_dataset.createVariable("J", J.dtype, ("time", "lat", "lon"))[:] = J
        output_dataset.createVariable("Jv", Jv.dtype, ("time", "lat", "lon"))[:] = Jv
        
        # Save the changes and close the output dataset
        output_dataset.sync()

ds = netCDF4.Dataset("iso_result_all.nc")
ds.set_auto_mask(True)

ds.variables
ds['time'][:]
var = ds["J"][:]
dar = ds["J_times_IspS"][:]

# Create a figure with two subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# First subplot
im1 = axs[0].imshow(var[0,:,:], origin="lower", extent=[-180, 180, -90, 90])
axs[0].set_title("Iso_2005")
fig.colorbar(im1, ax=axs[0], fraction=0.022, pad=0.03)

# Second subplot
im2 = axs[1].imshow(dar[12,:,:], origin="lower", extent=[-180, 180, -90, 90])
axs[1].set_title("Iso_2006")
fig.colorbar(im2, ax=axs[1], fraction=0.022, pad=0.03)

# Adjust the spacing between subplots
plt.tight_layout()

# Display the plot
plt.show()