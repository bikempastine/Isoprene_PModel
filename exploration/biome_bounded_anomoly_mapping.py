"""
Normlise and average data and them mask by biome.
"""
# Name path to the data
ISO_PATH = "../results/iso_result_OMI_filterd_mon_mean.nc"
HCHO_PATH = "../data/OMI_iso_estimate/mon_average_OMI.nc"

# Import required packages
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
import calendar
import cartopy.crs as ccrs
import matplotlib.colors as colors

# Define the functions
def filter_with_nan(array, value, land_cover_int):
    """
    Filter the array by replacing certain values with NaN.

    Args:
    - array: NumPy array to filter
    - value: Value to keep in the array, other values will be replaced with NaN
    - land_cover_int: Land cover array to identify values to be replaced with NaN

    Returns:
    - filtered_array: Filtered array with NaN values
    """
    filtered_array = np.copy(array)

    # Remove water, urban, and unvegetated regions
    filtered_array[np.where((land_cover_int == 0) | (land_cover_int == 9) | (land_cover_int == 10))] = np.nan

    # Mask the values not indicated
    filtered_array[np.where(land_cover_int != value)] = np.nan

    # Error checking
    nan_count_before = np.sum(np.isnan(array))
    nan_count_after = np.sum(np.isnan(filtered_array))
    print(value,"Number of NaN values before:", nan_count_before)
    print(value,"Number of NaN values after:", nan_count_after)

    return filtered_array


def gridcell_to_m2(latitude):
    """
    Calculate the number of square meters in each grid cell in a 0.5 by 0.5 degree grid.

    Args:
    - latitude: Latitude value in degrees

    Returns:
    - meter_square_gridcell: Number of square meters in each grid cell
    """
    # Latitude lengths stay at about 111.1 km per degree
    half_degree_lat = 111111.1 / 2

    # Equation to get the length of each half degree of longitude (changes with latitude)
    half_degree_lon = np.cos(np.deg2rad(latitude)) * (111111.1 / 2)

    meter_square_gridcell = half_degree_lat * half_degree_lon
    return meter_square_gridcell


# Import land cover data
land_cover_dataset = netCDF4.Dataset("../data/LandCover_half.nc")
land_cover_dataset.set_auto_mask(True)

# Get the land cover variable
land_cover = land_cover_dataset["Land Cover"][:]

# Fix land cover, to have North America on the left and China on the right
halfway = np.shape(land_cover)[1] // 2
first_half = land_cover[:, :halfway]
second_half = land_cover[:, halfway:]
land_cover = np.concatenate((second_half, first_half), axis=1)

# Convert land_cover to int
land_cover_int = land_cover.astype(int)

# Open the model data
iso_dataset = netCDF4.Dataset(ISO_PATH)
iso_dataset.set_auto_mask(True)

# Get the isoprene variable
model = iso_dataset["J_times_IspS"][:]

# Filter out the fill values
model[model > 1000] = np.nan

# Normalize the data between the maximum and 0
model_mean = np.nanmean(model, axis=0)
model_anomaly = np.subtract(model, model_mean)
model_norm = np.divide(model_anomaly, model_mean)

iso_dataset.close()

# Create an empty result array
model_result_array = np.empty((12, 8))

# Loop over i from 0 to 11 and j from 1 to 8
for i in range(12):
    for j in range(1, 9):
        filtered_array = filter_with_nan(model_norm[i, :, :], j, land_cover_int)
        model_result_array[i, j-1] = np.nanmean(filtered_array)

# Create a figure and a single subplot
fig, ax = plt.subplots(figsize=(10, 4))

# Plot the lines for different biomes
biome_labels = ['Grasses/Cereal', 'Shrubs', 'Broadleaf Crops', 'Savannah', 'Evergreen Broadleaf Forest',
                'Deciduous Broadleaf Forest', 'Evergreen Needleleaf Forest', 'Deciduous Needleleaf Forest']

for i in range(8):
    ax.plot(range(12), model_result_array[:, i], label=biome_labels[i])

# Set the title and legend
ax.set_title('Average anomaly in isoprene emissions by biome')
ax.legend()

# Display the plot
plt.show()


# Import HCHO data
HCHO_dataset = netCDF4.Dataset(HCHO_PATH)
HCHO_dataset.set_auto_mask(False)

# Convert HCHO from isoprene/grid cell to per m2
latitudes = [float(lat) for lat in HCHO_dataset['lat'][:]]
no_m2_in_grid = [gridcell_to_m2(lat) for lat in latitudes]
tiled = np.tile(no_m2_in_grid, (len(HCHO_dataset['lon'][:]), 1)).T
isoprene_per_m2 = HCHO_dataset["EMworldC5H8"][:] / tiled
h_lat = HCHO_dataset["lat"][:]
HCHO_dataset.close()

# Calculate the threshold value for outliers
threshold = np.percentile(isoprene_per_m2, 99.98)

# Filter out ocean and places with 0 emission
isoprene_per_m2[isoprene_per_m2 == 0] = np.nan

# Set values above the threshold to NaN
isoprene_per_m2[isoprene_per_m2 > threshold] = np.nan

# Normalize the data
HCHO_mean = np.nanmean(isoprene_per_m2, axis=0)
HCHO_anomaly = np.subtract(isoprene_per_m2, HCHO_mean)
HCHO_norm = np.divide(HCHO_anomaly, HCHO_mean)

# Create an empty result array
result_array_HCHO = np.empty((12, 8))

# Loop over i from 0 to 11 and j from 1 to 8
for i in range(12):
    for j in range(1, 9):
        filtered_array = filter_with_nan(np.flipud(HCHO_norm)[i, :, :], j, land_cover_int)
        result_array_HCHO[i, j-1] = np.nanmean(filtered_array)

# Create a figure and a single subplot
fig, ax = plt.subplots(figsize=(10, 4))

# Plot the lines for different biomes
for i in range(8):
    ax.plot(range(12), result_array_HCHO[:, i], label=biome_labels[i])

# Set the title and legend
ax.set_title('Average anomaly in HCHO emissions by biome')
ax.legend()

# Display the plot
plt.show()

# Create a figure with subplots for each biome
fig, axs = plt.subplots(4, 2, figsize=(10, 12))

# Loop over the subplots and plot the data
for i, ax in enumerate(axs.flat):
    # Plot the model data with dashed line
    ax.plot(range(12), model_result_array[:, i], 'k-', label='Model')

    # Plot the satellite observation data with solid line
    ax.plot(range(12), result_array_HCHO[:, i], 'k--', label='SatelliteObservation')

    # Set the title for each subplot
    ax.set_title(f'Average anomaly in the {biome_labels[i]} biome', fontweight='bold')

    # Set the legend for each subplot
    ax.legend()

# Adjust the spacing between subplots
plt.tight_layout()

# Display the plot
plt.show()