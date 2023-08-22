import numpy as np
import netCDF4
import cartopy.crs as ccrs
from matplotlib import pyplot as plt

# Open the Biome data
Biome_cdf = netCDF4.Dataset("../data/Type3_LC0.5.nc")
Biome_cdf.set_auto_mask(False)
biome_flip = Biome_cdf['type 3'][:]
Biome_cdf.close()


halfway = np.shape(biome_flip)[1] // 2
first_half = biome_flip[:, :halfway]
second_half = biome_flip[:, halfway:]
biome = np.concatenate((second_half, first_half), axis=1)


fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
img = ax.imshow(category_mask, origin='lower', cmap='RdBu_r', extent=[-180, 180, -90, 90],)
ax.coastlines()
ax.gridlines()
plt.colorbar(img, ax=ax, orientation='horizontal', label='sm')
plt.show()



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


# Calculate the average isoprene factor for each category in biome
unique_categories = np.unique(biome)
accumulated_factors = np.bincount(biome, weights=isoprene_factor[0,:,:])
category_counts = np.bincount(biome)

a = isoprene_factor[1,:,:]
average_isoprene_factors = []
unique_categories = np.unique(biome)

for category in unique_categories:
    category_mask = (biome == category)
    average_factor = np.nanmean(a[category_mask])
    average_isoprene_factors.append(average_factor)

# Avoid division by zero in case some categories have no data
average_isoprene_factor = np.divide(accumulated_factors, category_counts, out=np.zeros_like(accumulated_factors), where=category_counts != 0)

# Print the results
for category, average_factor in zip(unique_categories, average_isoprene_factor):
    print(f"Category {category}: Average Isoprene Factor = {average_factor:.2f}")