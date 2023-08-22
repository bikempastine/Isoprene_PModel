
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import matplotlib.colors as colors

# Land cover
lc = netCDF4.Dataset("../data/LandCover_half.nc")
lc.set_auto_mask(True)
lc.variables
land_cover = lc["Land Cover"][:]

# fix landcover 
halfway = np.shape(land_cover)[1]//2
first_half = land_cover[:,:halfway]
second_half = land_cover[:,halfway:]
land_cover = np.concatenate((second_half,first_half), axis=1)
#land_cover[land_cover == 0] = np.nan 
#land_cover[land_cover > 10] = np.nan 

land_cover_int = land_cover.astype(int)  # Convert land_cover to int

unique_land_cover = np.unique(land_cover_int)

# plot
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
image = ax.imshow(land_cover_int , origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()
plt.show()




#ds = netCDF4.Dataset("../results/iso_result_mon_mean.nc")
ds = netCDF4.Dataset("../results/iso_result_OMI_filterd_mon_mean.nc")
ds.set_auto_mask(True)

dar = ds["J_times_IspS"][:]
dar[dar > 1000] = np.nan

dar_mean = np.nanmean(dar, axis=0)
dar_mean.shape

# Subset everything above 30 degrees North
def calculate_mean_subset(latitude_threshold,latitude, data):
    lat_indices = np.where(latitude > latitude_threshold)[0]
    data_subset = data[:, lat_indices, :]
    data_mean_subset = np.nanmean(data_subset, axis=0)
    anomaly_subset = np.subtract(data_subset, data_mean_subset)
    norm_subset = np.divide(anomaly_subset, data_mean_subset)
    mean_data_subset = np.column_stack((np.arange(12), [np.nanmean(norm_subset[i, :, :]) for i in range(12)]))
    return mean_data_subset

mean_iso_subset_30 = calculate_mean_subset(30, latitude_iso, dar)

print(mean_iso_subset_30)

calculate_mean_subset(30, latitude_iso, dar)

# Normalize the data between the maximum and 0
anomaly = np.subtract(dar, dar_mean)
norm = np.divide(anomaly, dar_mean)

np.nanmean(norm[5,:,:])

np.nanmin(norm[5,:,:])
np.nanmax(norm[5,:,:])

# Create a figure with multiple subplots
fig, axs = plt.subplots(3, 4, figsize=(16, 12), subplot_kw={'projection': ccrs.PlateCarree()})

# Iterate over the range of plots you want to display (0 to 11)
for i in range(12):
    # Calculate the subplot indices based on the iteration index
    row = i // 4  # Row index
    col = i % 4   # Column index

    # Plot the normalized data in the current subplot
    im = axs[row, col].imshow(norm[i, :, :], origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='RdYlGn', vmin=-4, vmax=4, )

    axs[row, col].set_title("Iso_{:02d}".format(i))  # Format the title with leading zeros

    # Add coastlines
    axs[row, col].coastlines()

# Create a colorbar for the figure
fig.colorbar(im, ax=axs, fraction=0.022, pad=0.03, location='bottom')

# Adjust the spacing between subplots
plt.tight_layout()

# Move the colorbar below the subplots
fig.subplots_adjust(bottom=0.15)

# Display the plot
plt.show()
ds.close()

np.where(land_cover_int == 1)
norm[5, np.where(land_cover_int == 1)]
mean_norm_values = []
for lc in unique_land_cover:
    indices = np.where(land_cover_int == lc)  # Find indices of matching land_cover values
    mean_norm = np.mean(norm[indices])    # Calculate mean norm for the current land_cover value
    mean_norm_values.append(mean_norm)

print("Mean norm values for each land_cover group:")
for lc, mean_norm in zip(unique_land_cover, mean_norm_values):
    print(f"Land Cover {lc}: {mean_norm}")

#### Compare to the HCHO
HCHO = netCDF4.Dataset("../data/OMI_iso_estimate/mon_average_OMI.nc")
HCHO.set_auto_mask(False)

## Conveert isoprene values from the units kg/gridcell/month to kg/m2/month
def gridcell_to_m2(latitude):
    "Find the number of m2 in each grid cell in a 0.5 by 0.5 degree grid"
    half_degree_lat = 111111.1/2 # latitude lengths stay at about 111.1km per degree
    half_degree_lon = np.cos(np.deg2rad(latitude)) * (111111.1/2) # equation to get the length of each half degree of longitude (changes with latitude)
    meter_square_gridcell = half_degree_lat * half_degree_lon
    return meter_square_gridcell

latitudes = [float(lat) for lat in HCHO['lat'][:]] #get each latitude degree as a float 
no_m2_in_grid = [gridcell_to_m2(lat) for lat in latitudes] #gets the number of m2 blocks in each grid cell
tiled = np.tile(no_m2_in_grid, (len(HCHO['lon'][:]), 1)).T #repeat across each longitude as distance remains the same across each latitude degree

# Get the isoprene emmited per m2
isoprene_per_m2 = (HCHO["EMworldC5H8"][:])/(tiled) 

# Calculate the threshold value for the top 5%
threshold = np.percentile(isoprene_per_m2, 99.98)

isoprene_per_m2[isoprene_per_m2 == 0] = np.nan

# Set values above the threshold to NaN
isoprene_per_m2[isoprene_per_m2 > threshold] = np.nan

HCHO_mean = np.nanmean(isoprene_per_m2, axis=0)
HCHO_anomaly = np.subtract(isoprene_per_m2, HCHO_mean)
HCHO_norm = np.divide(HCHO_anomaly, HCHO_mean)


i = 6
np.nanmean(HCHO_norm[i,:,:])
np.nanmin(HCHO_norm[i,:,:])

# Subset everything above 30 degrees North
latitude_HCHO = HCHO.variables["lat"][:]


def HCHO_std_anomoly_subset(latitude_threshold, latitude, data):
    lat_indices = np.where(latitude > latitude_threshold)[0]
    data_subset = np.flipud(data)[:, lat_indices, :]
    data_mean_subset = np.nanmean(data_subset, axis=0)
    data_anomaly_subset = np.subtract(data_subset, data_mean_subset)
    data_norm_subset = np.divide(data_anomaly_subset, data_mean_subset)
    std_anomoly_subset = np.column_stack((np.arange(12), [np.nanmean(data_norm_subset[i, :, :]) for i in range(12)]))
    return std_anomoly_subset

mean_HCHO_subset_30 = HCHO_std_anomoly_subset(30, latitude_HCHO, isoprene_per_m2)

print(mean_HCHO_subset_30)
HCHO_std_anomoly_subset(0, latitude_HCHO, isoprene_per_m2)




# Create a figure and a single subplot
fig, ax = plt.subplots(figsize=(10, 4))

# Plot the first line with a red color
ax.plot(mean_HCHO_subset_30[:, 0], mean_HCHO_subset_30[:, 1], 'r-', label='HCHO')

# Plot the second line with a blue color
ax.plot(mean_iso_subset_30[:, 0], mean_iso_subset_30[:, 1], 'b-', label='Model')

# Set the title and legend
ax.set_title('Average anomaly in the Northern temperate zone (>30 degrees)')
ax.legend()

# Display the plot
plt.show()

np.nanmin(HCHO_norm)

# Create a figure with multiple subplots
fig, axs = plt.subplots(3, 4, figsize=(16, 12), subplot_kw={'projection': ccrs.PlateCarree()})

# Iterate over the range of plots you want to display (0 to 11)
for i in range(12):
    # Calculate the subplot indices based on the iteration index
    row = i // 4  # Row index
    col = i % 4   # Column index

    # Plot the normalized data in the current subplot
    im = axs[row, col].imshow(np.flipud(HCHO_norm[i, :, :]), origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='RdYlGn')#, vmin=-4, vmax=4, )

    axs[row, col].set_title("OMI_{:02d}".format(i))  # Format the title with leading zeros

    # Add coastlines
    axs[row, col].coastlines()

# Create a colorbar for the figure
fig.colorbar(im, ax=axs, fraction=0.022, pad=0.03, location='bottom')

# Adjust the spacing between subplots
plt.tight_layout()

# Move the colorbar below the subplots
fig.subplots_adjust(bottom=0.15)

# Display the plot
plt.show()


HCHO.close()