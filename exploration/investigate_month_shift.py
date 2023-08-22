"""
Plotting the average anomolies in isoprene emission using three geographical bounds for the HCHO derived data and the modlled emissions.
"""

# Import the packages
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
import calendar

# Define functions
def calculate_mean_subset_m(latitude_min, latitude_max, latitude, data):
    " Calculate the mean annomoly per month for the provided latitude range"

    # Find indices where latitude is within the specified range
    lat_indices = np.where((latitude > latitude_min) & (latitude < latitude_max))[0]

    # Subset the data based on latitude indices
    data_subset = data[:, lat_indices, :]

    # Calculate the mean of the data subset
    data_mean_subset = np.nanmean(data_subset, axis=0)

    # Calculate the anomaly by subtracting the mean from the subset
    anomaly_subset = np.subtract(data_subset, data_mean_subset)

    # (Value - Average)/Average = normalized values
    norm_subset = np.divide(anomaly_subset, data_mean_subset)

    # Return the normalized values next to the index of the month so they don't get mixed up
    mean_data_subset = np.column_stack(
        (np.arange(12), [np.nanmean(norm_subset[i, :, :]) for i in range(12)])
    )

    return mean_data_subset

def gridcell_to_m2(latitude):
    """Calculate the number of square meters in each grid cell in a 0.5 by 0.5 degree grid."""
    # Latitude lengths stay at about 111.1 km per degree
    half_degree_lat = 111111.1 / 2

    # Equation to get the length of each half degree of longitude (changes with latitude)
    half_degree_lon = np.cos(np.deg2rad(latitude)) * (111111.1 / 2)

    meter_square_gridcell = half_degree_lat * half_degree_lon
    return meter_square_gridcell

### MODEL DATA ###

# Import the Model data
model = netCDF4.Dataset("../results/iso_result_OMI_filterd_mon_mean.nc")
model.set_auto_mask(True)

# Get the isoprene estimate data
m_iso = model["J_times_IspS"][:]
m_lat = model["lat"][:] # get the latitude data
m_time = model["time"][:]
model.close()



m_iso[m_iso > 1000] = np.nan # set filled values to nan

indices = np.argsort(m_time)

# Get the means per month 
mean_iso_subset_NH= calculate_mean_subset_m(23.45, 90, m_lat, m_iso)[indices] #tropic of cancer up
mean_iso_subset_equator = calculate_mean_subset_m(-23.45, 23.45, m_lat, m_iso)[indices] # between tropics of cancer and capricorn
mean_iso_subset_SH = calculate_mean_subset_m(-90, -23.45, m_lat, m_iso)[indices] #tropic of capricorn down



### HCHO DATA ###

#### Compare to the HCHO
HCHO = netCDF4.Dataset("../data/OMI_iso_estimate/mon_average_OMI.nc")
HCHO.set_auto_mask(False)

# Convert HCHO from isoprene/grid call to per m2
# Convert latitude values to float
latitudes = [float(lat) for lat in HCHO['lat'][:]]

# Calculate the number of square meters in each grid cell
no_m2_in_grid = [gridcell_to_m2(lat) for lat in latitudes]

# Repeat the number of square meters across each longitude as the distance remains the same across each latitude degree
tiled = np.tile(no_m2_in_grid, (len(HCHO['lon'][:]), 1)).T

# Calculate the isoprene emitted per square meter
isoprene_per_m2 = HCHO["EMworldC5H8"][:] / tiled

# get the latitude data
h_lat = HCHO["lat"][:]
h_time = HCHO["time"][:]
HCHO.variables
HCHO.close()

# Filter out the outlying HCHO data
threshold = np.percentile(isoprene_per_m2, 99.98) # Filter out what is deemed to be outliers
isoprene_per_m2[isoprene_per_m2 > threshold] = np.nan # Set values above the threshold to NaN

# Filter out isoprene emissions of 0
isoprene_per_m2[isoprene_per_m2 == 0] = np.nan

# Get the means per month 
mean_HCHO_subset_NH= calculate_mean_subset_m(23.45, 90, h_lat, np.flipud(isoprene_per_m2))
mean_HCHO_subset_equator = calculate_mean_subset_m(-23.45, 23.45, h_lat, np.flipud(isoprene_per_m2))
mean_HCHO_subset_SH = calculate_mean_subset_m(-90, -23.45, h_lat, np.flipud(isoprene_per_m2))

a = np.concatenate(( mean_iso_subset_NH[:, 1][4:],  mean_iso_subset_NH[:, 1][:4]))

## Plot the figures
# Create a figure and three subplots arranged vertically
fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# Plot for the first subplot
axs[0].plot(mean_HCHO_subset_NH[:, 0], mean_HCHO_subset_NH[:, 1], 'r-', label='HCHO')
axs[0].plot(mean_HCHO_subset_NH[:, 0], a, 'b-', label='Model')
axs[0].set_title('Average anomaly in the Northern temperate zone', fontweight='bold')

# Plot for the second subplot
axs[1].plot(mean_HCHO_subset_SH[:, 0], mean_HCHO_subset_SH[:, 1], 'r-', label='HCHO')
axs[1].plot(mean_HCHO_subset_SH[:, 0], mean_iso_subset_SH[:, 1], 'b-', label='Model')
axs[1].set_title('Average anomaly in the Southern temperate zone', fontweight='bold')

# Plot for the third subplot
axs[2].plot(mean_HCHO_subset_equator[:, 0], mean_HCHO_subset_equator[:, 1], 'r-', label='HCHO')
axs[2].plot(mean_HCHO_subset_equator[:, 0], mean_iso_subset_equator[:, 1], 'b-', label='Model')
axs[2].set_title('Average anomaly in the equator', fontweight='bold')

# Set the same legend for all subplots
axs[0].legend()

# Set the same x-axis limits for all subplots
axs[0].set_xlim(0, 11)


# Set the x-axis ticks to be the month names and rotate them
for ax in axs:
    ax.set_xticks(range(12))
    ax.set_xticklabels(calendar.month_name[1:], rotation=45, ha='right')

# Add a dashed gray line at 0 for each subplot
for ax in axs:
    ax.axhline(0, color='gray', linestyle='--')

# Adjust the spacing between subplots
plt.subplots_adjust(hspace=0.4)

# Display the figure
plt.show()


