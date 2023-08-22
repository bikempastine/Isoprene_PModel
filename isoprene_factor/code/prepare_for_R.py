import numpy as np
import netCDF4
import xarray as xr
import pandas as pd


model_cdf = netCDF4.Dataset("../../results/iso_result_OMI_filtered.nc")
model_cdf.set_auto_mask(False)

J = model_cdf.variables['J'][:]
vcmax = model_cdf.variables['vcmax'][:]
lue = model_cdf.variables['lue'][:]
jmax = model_cdf.variables['jmax'][:]
gpp = model_cdf.variables['gpp'][:]

sm_continuous = np.load('../data/sm_stress_2005_2014.npy')
sm_categorical = np.where(sm_continuous == 1, 1, 2) 

F = np.load('../data/isoprene_factor_2005_2014.npy')
F = F* 1e6

fpar_cdf = netCDF4.Dataset("../../temporary/fpar_all.nc")
fpar_cdf.set_auto_mask(False)
fpar = fpar_cdf.variables['FPAR'][:]
fpar_lon_flip = fpar_cdf['lon'][:]
fpar[(fpar < 0)| (fpar > 1)] = np.nan

# Rearrange FPAR data to account for the longitude split
halfway = np.shape(fpar)[2] // 2
first_half = fpar[:, :, :halfway]
second_half = fpar[:, :, halfway:]
fpar_fixed = np.concatenate((second_half, first_half), axis=2)
fpar_cdf.close()

# Load the biome map
biome_cdf = netCDF4.Dataset("../data/clean_Type3_biome.nc")
biome = biome_cdf['type 3'][:]
biome_cdf.close()
biome[(biome < 1) | (biome > 9)] = np.nan
biome_rep = np.tile(biome[np.newaxis, :, :], (120, 1, 1))

# Load air temperature
temp_cdf = netCDF4.Dataset("../data/Tair_2005_2014.nc")
temp = temp_cdf['Tair'][:] - 273.15
temp[(temp > 1000)] = np.nan
temp_cdf.close()

# Load swdown
swdown_cdf = netCDF4.Dataset("../data/SWdown_2005_2014.nc")
swdown = swdown_cdf['SWdown'][:]
swdown[(swdown > 1000)] = np.nan
swdown_cdf.close()


# Create a new NetCDF file and write the features to it
new_nc = netCDF4.Dataset('../data/regression_data.nc', 'w')

# Copy over dimensions from the original file
for name, dimension in model_cdf.dimensions.items():
    new_nc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

# Create variables in the new file and set data

new_F = new_nc.createVariable('F', np.float32, model_cdf.variables['J'].dimensions)
new_F[:] = F

new_J = new_nc.createVariable('J', np.float32, model_cdf.variables['J'].dimensions)
new_J[:] = J

new_lue = new_nc.createVariable('lue', np.float32, model_cdf.variables['J'].dimensions)
new_lue[:] = lue

new_jmax = new_nc.createVariable('jmax', np.float32, model_cdf.variables['J'].dimensions)
new_jmax[:] = jmax

new_gpp = new_nc.createVariable('gpp', np.float32, model_cdf.variables['J'].dimensions)
new_gpp[:] = gpp

new_vcmax = new_nc.createVariable('vcmax', np.float32, model_cdf.variables['J'].dimensions)
new_vcmax[:] = vcmax

new_sm_continuous = new_nc.createVariable('sm_continuous', np.float32, model_cdf.variables['J'].dimensions)
new_sm_continuous[:] = sm_continuous

new_sm_categorical = new_nc.createVariable('sm_categorical', np.int, model_cdf.variables['J'].dimensions)
new_sm_categorical[:] = sm_categorical

new_biome = new_nc.createVariable('biome', np.int, model_cdf.variables['J'].dimensions)
new_biome[:] = biome_rep

new_temp = new_nc.createVariable('temp', np.float32, model_cdf.variables['J'].dimensions)
new_temp[:] = temp

new_swdown = new_nc.createVariable('swdown', np.float32, model_cdf.variables['J'].dimensions)
new_swdown[:] = swdown

new_fpar = new_nc.createVariable('fpar', np.float32, model_cdf.variables['J'].dimensions)
new_fpar[:] = fpar_fixed

new_nc.close()

# Close the original NetCDF file
model_cdf.close()


# Load your NetCDF data
ds = xr.open_dataset('../data/regression_data.nc')

# Convert the xarray dataset to a pandas DataFrame
df = ds.to_dataframe().reset_index()

# Filter out rows with NaN values in any column
df = df[df['biome'] > 0]

df_clean = df.dropna()

# Save the cleaned data to a CSV file (for R)
df_clean.to_csv('../data/cleaned_regress_for_R.csv', index=False)

y = np.log(df_clean['F'])
norm_y = (y - np.mean(y)) / np.std(y)
plt.hist(norm_y, bins=60, edgecolor='black')
plt.show()

import statsmodels.api as sm
import matplotlib.pyplot as plt

sm.qqplot(norm_y, line='s')
plt.title('Q-Q Plot')
plt.show()

