
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize

#ds = netCDF4.Dataset("../results/iso_result_mon_mean.nc")
ds = netCDF4.Dataset("../results/iso_result_OMI_filterd_mon_mean.nc")
ds.set_auto_mask(True)

dar = ds["J_times_IspS"][:]
dar[dar > 1000] = np.nan

# Calculate the absolute maximum value
dar_max = np.nanmax(dar)

# Normalize the data between the maximum and 0
norm = np.divide(dar, dar_max)

# Create a figure with multiple subplots
fig, axs = plt.subplots(3, 4, figsize=(16, 12))

# Iterate over the range of plots you want to display (0 to 11)
for i in range(12):
    # Calculate the subplot indices based on the iteration index
    row = i // 4  # Row index
    col = i % 4   # Column index


    # Plot the normalized data in the current subplot
    im = axs[row, col].imshow(norm[i, :, :], origin="lower", extent=[-180, 180, -90, 90])
    axs[row, col].set_title("Iso_{:02d}".format(i))  # Format the title with leading zeros

# Create a colorbar for the figure
fig.colorbar(im, ax=axs, fraction=0.022, pad=0.03,location='bottom')

# Adjust the spacing between subplots
plt.tight_layout()

# Move the colorbar below the subplots
fig.subplots_adjust(bottom=0.15)

# Display the plot
plt.show()

ds.close()

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

np.nanmean(isoprene_per_m2)
#isoprene_per_m2[isoprene_per_m2 == 0] = np.nan
#isoprene_per_m2[dar == np.nan] = np.nan

# Calculate the absolute maximum value
isoprene_per_m2_max = np.nanmax(isoprene_per_m2)

# Normalize the data between the maximum and 0
norm_HCHO = np.divide(isoprene_per_m2, isoprene_per_m2_max)
#norm_HCHO[np.isnan(np.flipud(dar))] = np.nan

# Create a figure with multiple subplots
fig, axs = plt.subplots(3, 4, figsize=(16, 12))

# Iterate over the range of plots you want to display (0 to 11)
for i in range(12):
    # Calculate the subplot indices based on the iteration index
    row = i // 4  # Row index
    col = i % 4   # Column index


    # Plot the normalized data in the current subplot
    im = axs[row, col].imshow(np.flipud(norm_HCHO[i, :, :]), origin="lower", extent=[-180, 180, -90, 90])
    axs[row, col].set_title("OMI_{:02d}".format(i))  # Format the title with leading zeros

# Create a colorbar for the figure
fig.colorbar(im, ax=axs, fraction=0.022, pad=0.03,location='bottom')

# Adjust the spacing between subplots
plt.tight_layout()

# Move the colorbar below the subplots
fig.subplots_adjust(bottom=0.15)

# Display the plot
plt.show()

HCHO.close()