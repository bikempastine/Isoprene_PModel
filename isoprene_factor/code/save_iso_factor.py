import numpy as np
import netCDF4

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

# calculate the isoprene factor
isoprene_factor = omi/model

# Save the data to a npy file
np.save('../data/isoprene_factor_2005_2014.npy', isoprene_factor)