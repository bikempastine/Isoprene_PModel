import numpy as np
import netCDF4
import cartopy.crs as ccrs
from matplotlib import pyplot as plt

# Open the Model data
Model_cdf = netCDF4.Dataset("../data/model_clean.nc")
Model_cdf.set_auto_mask(False)
model = Model_cdf['model_data'][:]
Model_cdf.close()

# Open the OMI data
OMI_cdf = netCDF4.Dataset("../data/OMI_clean.nc")
OMI_cdf.set_auto_mask(False)
omi = OMI_cdf['omi_iso'][:]
OMI_cdf.close()

# filter out when Model is above 1 (I think here the error is driving the high F values)
model[model < 1] = np.nan
omi[model < 1] = np.nan

model[omi < 3e-05] = np.nan
omi[omi < 3e-05] = np.nan

# calculate the isoprene factor
isoprene_factor = omi/model
np.nanmax(isoprene_factor)

# visulalise the factor in  a histogram
f = isoprene_factor.flatten()
plt.hist(f, bins=5000)
plt.show()

# visualise in a scatter plot
plt.scatter(isoprene_factor, model, color='red',alpha=0.01, s=1)
plt.ylabel('model')
plt.xlabel('isoprene_factor')
plt.title('Scatter Plot of model vs isoprene_factor')
plt.tight_layout()
plt.show()

plt.scatter(omi, isoprene_factor, color='blue',alpha=0.01, s=1)
plt.xlabel('omi')
plt.ylabel('isoprene_factor')
plt.title('Scatter Plot of omi vs isoprene_factor')
plt.show()


