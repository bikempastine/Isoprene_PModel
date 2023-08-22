"""
Investigates the difference between the modelled (J * f(T)) and HCHO data for January 2010
"""

# Imports
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors

# Open HCHO data
HCHO_ncdf = netCDF4.Dataset("../results/HCHO_Jan_2010.nc")
HCHO_ncdf.set_auto_mask(False)
HCHO = HCHO_ncdf["HCHO_Iso"][0, :, :] # January, 2010, converted to kg/m2/month
HCHO = np.flipud(HCHO)  # need to flip because it comes upsidedown
HCHO_ncdf.close()

#get rid of outliers

def set_outliers_to_nan(arr):
    """Get rid of the top 2% of the data"""
    # Calculate the threshold value
    threshold = np.percentile(arr, 99.9)

    # Set values above the threshold to NaN
    arr[arr > threshold] = np.nan

    return arr

HCHO = set_outliers_to_nan(HCHO)

# Open model data
model_ncdf = netCDF4.Dataset("../results/Iso_test_results.nc")
model_ncdf.set_auto_mask(False)
model = model_ncdf["J_times_IspS"][0, :, :]
model_ncdf.close()

# Set values equal to 0 to NA : an ad hoc type of masking
def set_nan(data):
    masked = np.where((data == 0), np.nan, data)
    return masked

masked_HCHO = set_nan(HCHO)
masked_model = set_nan(model)

### Compare the two models directly, side by side 
# Create a new figure and axis
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Plot the first dataset on the first subplot
image1 = ax1.imshow(masked_model, origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='PuRd')

# Plot the second dataset on the second subplot
image2 = ax2.imshow(masked_HCHO, origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='PuRd')
ax2.coastlines()

# Add colorbar to the second subplot
cbar = plt.colorbar(image2, ax=[ax1, ax2], orientation='horizontal', pad=0.05)
cbar.set_label('micromols m-2 s-1')

# Set the map extent for both subplots
ax1.set_global()
ax2.set_global()

# Add titles and labels
ax1.set_title('Model')
ax2.set_title('HCHO')
ax1.set_xlabel('Longitude')
ax2.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')

plt.suptitle('Comparison of Datasets', fontsize=14)
plt.show()


###  Devide Model / HCHO
model_over_HCHO = masked_model / masked_HCHO

# Plotting the histogram
plt.hist(model_over_HCHO, bins=5, edgecolor='black')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Model / HCHO')
plt.show()

###  Subtract Model from HCHO 
diff = np.subtract(masked_HCHO, masked_model)

# Create a scatter plot
plt.scatter(x=masked_HCHO, y=diff, alpha=0.7)
plt.xlabel('masked_HCHO')
plt.ylabel('HCHO - Model')
plt.show()

# Plotting a histogram of diff
plt.hist(diff, bins=5, edgecolor='black')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('HCHO - model')
plt.show()

# Plot diff gloablly
# Create a new figure and axis
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot the data on the map
image = ax.imshow(diff, origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
ax.coastlines()

# Add a colorbar
cbar = plt.colorbar(image, ax=ax)
cbar.set_label('kg/ m2 / month')

# Set the map extent
ax.set_global()

# Add title and labels
plt.title('HCHO - model January 2010')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()


## Get the histogram values for diff
# Define the number of bins
num_bins = 10

# Calculate the histogram
diff_clean = diff[~np.isnan(diff)]
diff_clean_range = diff_clean[diff_clean < 0.0033]

counts, bins = np.histogram(diff_clean, bins=num_bins)
for i in range(len(counts)):
    bin_range = f"{bins[i]:.6f} - {bins[i + 1]:.6f}"
    print(f"Bin {i}: Range: {bin_range}, Count: {counts[i]}")

count = np.count_nonzero(diff_clean > 0)
print(count)


### Normalise between 0 and 1

def normalize_0_1(arr):
    """"Normalises arrays between 0 and 1"""
    min_value = np.nanmin(arr)
    max_value = np.nanmax(arr)
    normalized_array = (arr - min_value) / (max_value - min_value)
    return normalized_array


# Normalise the model isoprene data between 0 and 1 
normalised_model = normalize_0_1(masked_model)
np.nanmax(normalised_model)
np.nanmean(normalised_model)

# Normalise HCHO between 0 and 1
normalised_HCHO = normalize_0_1(masked_HCHO)
np.nanmax(normalised_HCHO)
np.nanmean(normalised_HCHO)



plt.scatter(normalised_model, normalised_HCHO)
plt.xlabel('model')
plt.ylabel('HCHO')
plt.show()

np.array(normalised_model).flatten()

np.polyfit(np.array(normalised_model).flatten(), np.array(normalised_HCHO).flatten, 1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Plot the first dataset on the first subplot
image1 = ax1.imshow(normalised_model, origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='PuRd')
ax1.coastlines()

# Plot the second dataset on the second subplot
image2 = ax2.imshow(normalised_HCHO, origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='PuRd')
ax2.coastlines()

# Add colorbar to the second subplot
cbar = plt.colorbar(image2, ax=[ax1, ax2], orientation='horizontal', pad=0.05)
cbar.set_label('normalised value')

# Set the map extent for both subplots
ax1.set_global()
ax2.set_global()

# Add titles and labels
ax1.set_title('Model')
ax2.set_title('HCHO')
ax1.set_xlabel('Longitude')
ax2.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')

plt.suptitle('Comparison of Datasets', fontsize=14)
plt.show()


relationship = normalised_model/normalised_HCHO

fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot the data on the map
image = ax.imshow(relationship, origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
ax.coastlines()

# Add a colorbar
cbar = plt.colorbar(image, ax=ax)
cbar.set_label('model/HCHO')

# Set the map extent
ax.set_global()

# Add title and labels
plt.title('model/HCHO January 2010')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

np.nanmean(relationship)
np.nanmax(relationship)
np.nanmin(relationship)