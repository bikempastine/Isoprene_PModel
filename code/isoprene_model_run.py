"""
Test case for the isoprene model, getting and saving critical variables
"""

## DEFINE VARIABLES
REFERANCE_TEMP = 30 #reference temperature to run the IspS_function at
INPUT_FILE = "../data/pmodel_inputs.nc"
OUTPUT_FILE = "../results/Iso_test_results.nc"


# Load the functions for J, Jv, IspS
import isoprene_functions as Iso_fun

# Load in the pyrealm package
from pyrealm import pmodel

# Load in the required packages
import math
import netCDF4
import numpy as np
from matplotlib import pyplot as plt


# Load an example dataset containing the main variables.
ds = netCDF4.Dataset(INPUT_FILE)
ds.set_auto_mask(False)


# Extract the six variables for all months
temp = ds["temp"][:]
co2 = ds["CO2"][:]  # Note - spatially constant 
elev = ds["elevation"][:] # Note - temporally constant
vpd = ds["VPD"][:]
fapar = ds["fAPAR"][:]
ppfd = ds["ppfd"][:]

# Manipulation of the inputs to fit the assumptions of the pmodel
patm = pmodel.calc_patm(elev) # Convert elevation into  atmospheric pressure
temp[temp < -25] = np.nan # Remove temperature values below -25 degrees
vpd[vpd < 0] = 0 #Clip Vapour pressure defficit so that VPD below zero is set to zero

# Run the PModel
env = pmodel.PModelEnvironment(tc=temp, co2 = co2, patm=patm, vpd=vpd, )
model = pmodel.PModel(env, method_jmaxlim = 'smith19') #run with smith19 to be consistent with J calculation
model.estimate_productivity(fapar, ppfd)

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
IspS = np.vectorize(Iso_fun.IspS_function)(temp, REFERANCE_TEMP)

# Produce objects that will be stored
GPP = model.gpp
Jmax = model.jmax
LUE = model.lue
Vcmax = model.vcmax

J_times_IspS_umol = J*IspS

# Create a new netCDF file for storing the variables
output_file = netCDF4.Dataset(OUTPUT_FILE, "w")

# Create dimensions in the output netCDF file
output_file.createDimension('lat', ds.dimensions["latitude"].size)  
output_file.createDimension('lon', ds.dimensions["longitude"].size)  
output_file.createDimension('time', None)  # unlimited axis (can be appended to).

# Copy dimention variables to the output file from ds
output_file.createVariable("time", ds["Time"].dtype, ("time",))[:] = ds.variables['Time'][:]
output_file.createVariable("lat", ds["latitude"].dtype, ("lat",))[:] = ds.variables['latitude'][:]
output_file.createVariable("lon", ds["longitude"].dtype, ("lon",))[:] = ds.variables['longitude'][:]

# Create and assign data to variables in the output file
output_file.createVariable("gpp", GPP.dtype, ("time", "lat", "lon"))[:] = GPP
output_file.createVariable("jmax", Jmax.dtype, ("time", "lat", "lon"))[:] = Jmax
output_file.createVariable("lue", LUE.dtype, ("time", "lat", "lon"))[:] = LUE
output_file.createVariable("vcmax", Jv.dtype, ("time", "lat", "lon"))[:] = Vcmax
output_file.createVariable("J_times_IspS", J_times_IspS_umol.dtype, ("time", "lat", "lon"))[:] = J_times_IspS_umol
output_file.createVariable("J", J.dtype, ("time", "lat", "lon"))[:] = J
output_file.createVariable("Jv", Jv.dtype, ("time", "lat", "lon"))[:] = Jv


# Close the netCDF to prevent corruption
ds.close()
output_file.close()

# Print to the console that the run is complete
print("Done and saved in: ", OUTPUT_FILE )