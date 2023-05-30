"""Test case for the isoprene model, getting and saving critical variables"""

# Load the functions for J, Jv, IspS
import J_function as J_fun
import IspS_functions as IspS

# Load in the pyrealm package
from pyrealm import pmodel

# Load in the required packages
import netCDF4
import numpy as np
from matplotlib import pyplot as plt

## Define Variables
REFERANCE_TEMP = 30 

# Load an example dataset containing the main variables.
ds = netCDF4.Dataset("../data/pmodel_inputs.nc")
ds.set_auto_mask(False)

# Extract the six variables for all months
temp = ds["temp"][:]
co2 = ds["CO2"][:]  # Note - spatially constant 
elev = ds["elevation"][:] # Note - temporally constant
vpd = ds["VPD"][:]
fapar = ds["fAPAR"][:]
ppfd = ds["ppfd"][:]

ds.close()

# Manipulation to fit the assumptions of the pmodel
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

# Run the functions to find J and Jv
J = J_fun.find_J(a_j, ci, gammastar)
Jv = J_fun.find_Jv(vcmax, ci, gammastar, kmm)

# Run the IspS function
IspS = np.vectorize(IspS.IspS_function)(temp, REFERANCE_TEMP)

# Produce objects that will be stored together
J_times_IspS = J*IspS
GPP = model.gpp
Jmax = model.jmax
LUE = model.lue

# Create a new netCDF file for storing the variables
output_file = netCDF4.Dataset("../results/Iso_test_results.nc", "w")

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
output_file.createVariable("J_times_IspS", J_times_IspS.dtype, ("time", "lat", "lon"))[:] = J_times_IspS
output_file.createVariable("J", J.dtype, ("time", "lat", "lon"))[:] = J
output_file.createVariable("Jv", Jv.dtype, ("time", "lat", "lon"))[:] = Jv



# Little test
test_cdf = netCDF4.Dataset("../results/Iso_test_results.nc")
test_cdf.set_auto_mask(False)

#plot
ax = plt.imshow(test_cdf['Jv'][0,:,:], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()