import J_function as J_fun #import the functions from the J_function.py script
import netCDF4
from pyrealm import pmodel
import numpy as np
from matplotlib import pyplot as plt

# Load an example dataset containing the main variables.
ds = netCDF4.Dataset("../data/pmodel_inputs.nc")
ds.set_auto_mask(False)

# Extract the six variables for all months
temp = ds["temp"][:]
co2 = ds["CO2"][:]  # Note - spatially constant but mapped.
elev = ds["elevation"][:]  # Note - temporally constant but repeated #no need to find more just repeat 
vpd = ds["VPD"][:]
fapar = ds["fAPAR"][:]
ppfd = ds["ppfd"][:]

ds.close()

# Need to do some manipulation to fit the assumptions of the pmodel
patm = pmodel.calc_patm(elev) # Convert elevation into  atmospheric pressure
temp[temp < -25] = np.nan # Remove temperature values below -25 degrees
vpd[vpd < 0] = 0 #Clip Vapour pressure defficit so that VPD below zero is set to zero

### Calculate the photosynthetic environment
env = pmodel.PModelEnvironment(tc=temp, co2 = co2, patm=patm, vpd=vpd, )
env.summarize()

### Run the P model
model = pmodel.PModel(env)
model.estimate_productivity(fapar, ppfd)
model.summarize()


# Get the nessisary variables to run the models
gammastar = env.gammastar
ci = model.optchi.ci
vcmax = model.vcmax
kmm = env.kmm
a_j = (model.vcmax / model.optchi.mjoc) * model.optchi.mj #getting ths from the code pyrealm/pyrealm/pmodel/pmodel.py  line 613 and 595

gpp = model.gpp

### Run the functions to find J and Jv
J = J_fun.find_J(a_j, ci, gammastar)
Jv = J_fun.find_Jv(vcmax, ci, gammastar, kmm)

### Find out about J results
np.nanmean(J)

#month 1
ax = plt.imshow(J[0, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()

#month 2
ax = plt.imshow(J[1, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()

#difference between two months
ax = plt.imshow(J[0, :, :] - J[1, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()

#find out about JV
np.nanmean(Jv)

#month 1
ax = plt.imshow(Jv[0, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()

#month 2
ax = plt.imshow(Jv[1, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()

#difference between two months
ax = plt.imshow(Jv[0, :, :] - Jv[1, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()

# J - Jv
ax = plt.imshow(J[0, :, :] - Jv[0, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()

ax = plt.imshow(J[1, :, :] - Jv[1, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()


# Plot J and Jv against ppfd in europe

# Define the extent of Europe
europe_extent = [-30, 45, 25, 75]  # [min_lon, max_lon, min_lat, max_lat]


ds = netCDF4.Dataset("../data/pmodel_inputs.nc")
ds.set_auto_mask(False)

# Find the indices corresponding to Europe extent
lon_indices = np.where((ds['longitude'][:] >= europe_extent[0]) & (ds['longitude'][:] <= europe_extent[1]))[0]
lat_indices = np.where((ds['latitude'][:] >= europe_extent[2]) & (ds['latitude'][:] <= europe_extent[3]))[0]

ds.close()
# Extract the J and fapar values for Europe
J_europe = J[:, lat_indices, :][:, :, lon_indices]
Jv_europe = Jv[:, lat_indices, :][:, :, lon_indices]
ppfd_europe = ppfd[:, lat_indices, :][:, :, lon_indices]

# Reshape the arrays for plotting
J_europe_flat = J_europe[0,:,:].flatten()
ppfd_europe_flat = ppfd_europe[0,:,:].flatten()
Jv_europe_flat = Jv_europe[0,:,:].flatten()

# Plot J against fapar for Europe
plt.scatter(ppfd_europe_flat, J_europe_flat)
plt.xlabel('ppfd')
plt.ylabel('J')
plt.title('J vs ppfd (Europe)')
plt.show()

plt.scatter(ppfd_europe_flat, Jv_europe_flat)
plt.xlabel('ppfd')
plt.ylabel('Jv')
plt.title('Jv vs ppfd (Europe)')
plt.show()

