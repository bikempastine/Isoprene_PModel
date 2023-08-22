import netCDF4
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import matplotlib.colors as colors

# MODEL RESULTS
ds = netCDF4.Dataset("../results/iso_result_OMI_filterd_mon_mean.nc")
ds.set_auto_mask(True)

iso = ds["J_times_IspS"][:]
iso[iso > 1000] = np.nan

iso_mean = np.nanmean(iso, axis=0)
iso_mean.shape

# Normalize the data between the maximum and 0
anomaly = np.subtract(iso, iso_mean)
norm = np.divide(anomaly, iso_mean)

# Compare to the HCHO
HCHO = netCDF4.Dataset("../data/OMI_iso_estimate/mon_average_OMI.nc")
HCHO.set_auto_mask(False)

# Convert isoprene values from the units kg/gridcell/month to kg/m2/month
def gridcell_to_m2(latitude):
    "Find the number of m2 in each grid cell in a 0.5 by 0.5 degree grid"
    half_degree_lat = 111111.1/2  # latitude lengths stay at about 111.1km per degree
    half_degree_lon = np.cos(np.deg2rad(latitude)) * (111111.1/2)  # equation to get the length of each half degree of longitude (changes with latitude)
    meter_square_gridcell = half_degree_lat * half_degree_lon
    return meter_square_gridcell

latitudes = [float(lat) for lat in HCHO['lat'][:]]  # get each latitude degree as a float
no_m2_in_grid = [gridcell_to_m2(lat) for lat in latitudes]  # gets the number of m2 blocks in each grid cell
tiled = np.tile(no_m2_in_grid, (len(HCHO['lon'][:]), 1)).T  # repeat across each longitude as the distance remains the same across each latitude degree

# Get the isoprene emitted per m2
isoprene_per_m2 = HCHO["EMworldC5H8"][:] / tiled

# Calculate the threshold value for the top 5%
threshold = np.percentile(isoprene_per_m2, 99.98)

isoprene_per_m2[isoprene_per_m2 == 0] = np.nan

# Set values above the threshold to NaN
isoprene_per_m2[isoprene_per_m2 > threshold] = np.nan

HCHO_mean = np.nanmean(isoprene_per_m2, axis=0)
HCHO_anomaly = np.subtract(isoprene_per_m2, HCHO_mean)
HCHO_norm = np.divide(HCHO_anomaly, HCHO_mean)

# # Create a figure with multiple subplots
fig, axs = plt.subplots(6, 4, figsize=(16, 24), subplot_kw={'projection': ccrs.PlateCarree()})

# Iterate over the range of plots you want to display (0 to 11)
for i in range(6):
    # Calculate the subplot indices based on the iteration index
    row = i  # Row index for the first set of plots
    col = 0   # Column index for the OMI plot

    # Plot the normalized OMI data in the current subplot
    im_omi = axs[row, col].imshow(np.flipud(HCHO_norm[i, :, :]), origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='RdYlGn', vmin=-4, vmax=4)

    axs[row, col].set_title("OMI_{:02d}".format(i))  # Format the title with leading zeros

    # Add coastlines
    axs[row, col].coastlines()

    # Calculate the subplot indices for the second set of plots
    row = i  # Row index for the second set of plots
    col = 1   # Column index for the isoprene plot

    # Plot the normalized isoprene data in the current subplot
    im_iso = axs[row, col].imshow(norm[i, :, :], origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='RdYlGn', vmin=-4, vmax=4)

    axs[row, col].set_title("Iso_{:02d}".format(i))  # Format the title with leading zeros

    # Add coastlines
    axs[row, col].coastlines()


for i in range(6, 12):
    # Calculate the subplot indices based on the iteration index
    row = i-6  # Row index for the first set of plots
    col = 2   # Column index for the OMI plot

    # Plot the normalized OMI data in the current subplot
    im_omi = axs[row, col].imshow(np.flipud(HCHO_norm[i, :, :]), origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='RdYlGn', vmin=-4, vmax=4)

    axs[row, col].set_title("OMI_{:02d}".format(i))  # Format the title with leading zeros

    # Add coastlines
    axs[row, col].coastlines()

    # Calculate the subplot indices for the second set of plots
    row = i-6  # Row index for the second set of plots
    col = 3   # Column index for the isoprene plot

    # Plot the normalized isoprene data in the current subplot
    im_iso = axs[row, col].imshow(norm[i, :, :], origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), cmap='RdYlGn', vmin=-4, vmax=4)

    axs[row, col].set_title("Iso_{:02d}".format(i))  # Format the title with leading zeros

    # Add coastlines
    axs[row, col].coastlines()

# Create a colorbar for the figure
cbar_ax = fig.add_axes([0.3, 0.05, 0.4, 0.02])  # Adjust the position of the colorbar
fig.colorbar(im_iso, cax=cbar_ax, orientation='horizontal')

# Adjust the spacing between subplots
fig.subplots_adjust(hspace=0.3, wspace=-0.7)

# Display the plot
plt.show()

ds.close()
HCHO.close()

