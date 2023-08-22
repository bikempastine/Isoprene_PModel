"""
Test case for the isoprene model, getting and saving critical variables
"""

## DEFINE VARIABLES
REFERANCE_TEMP = 30 #reference temperature to run the IspS_function at
#INPUT_FILE = "../data/pmodel_inputs.nc"
#INPUT_FILE = "../data/Pmodel_data_full.nc"
OUTPUT_FILE = "../results/Iso_results.nc"
INPUT_FILE = "../temporary/Pmodel_data.nc"
INPUT_FILE = "../temporary/Pmodel_4_of_6_data.nc"
# Load the functions for J, Jv, IspS
import isoprene_functions as Iso_fun

# Load in the pyrealm package
from pyrealm import pmodel
from pyrealm import hygro

# Load in the required packages
import math
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


# Load an example dataset containing the main variables.
# fpr = netCDF4.Dataset("../data/Pmodel_data.nc")
# fpr.set_auto_mask(True)

ds = netCDF4.Dataset(INPUT_FILE)
ds.set_auto_mask(True)

ds.variables

# Extract the six variables for all months
temp = ds["Tair"][:] #Celcius
fapar_wrong = fpr["FPAR"][:] #Weird projectin, ocean is above 1
ppfd = ds["SWdown"][:] 
elev = ds["Elev"][:]
vap = ds["vap"][:]
co2 = ds["CO2"][:]
fapar_wrong = ds["FPAR"][:]
#ds.close()
fpr.close()

np.max(fapar_wrong)
# Get VPD from VAP
vpd = hygro.convert_vp_to_vpd(vap*0.1 , temp) #vap in hPa, convert to hPA 

halfway = np.shape(fapar_wrong)[2]//2
first_half = fapar_wrong[:,:,:halfway]
second_half = fapar_wrong[:,:,halfway:]

fapar = np.concatenate((second_half,first_half), axis=2)
fapar = fapar.astype(float) 
fapar[fapar > 1] = np.nan

np.nanmax(fapar)
fapar_wrong.dtype

# plot the results
x = fapar
im = plt.imshow(x[0, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar(im, fraction=0.022, pad=0.03)
plt.title("fpar");
plt.show()

# Manipulation of the inputs to fit the assumptions of the pmodel
patm = pmodel.calc_patm(elev) # Convert elevation into  atmospheric pressure
temp[temp < -25] = np.nan # Remove temperature values below -25 degrees
vpd[vpd < 0] = 0 #Clip Vapour pressure defficit so that VPD below zero is set to zero

# Fit the problem with fapar


# Run the PModel
env = pmodel.PModelEnvironment(tc=temp, co2 = co2, patm=patm, vpd=vpd, )
model = pmodel.PModel(env, method_jmaxlim = 'smith19') #run with smith19 to be consistent with J calculation
model.estimate_productivity(fapar, ppfd)

model.summarize()
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