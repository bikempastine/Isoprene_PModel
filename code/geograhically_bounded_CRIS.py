
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


# Import the Model data
CRIS = netCDF4.Dataset("../data/2012to2020_CrIS_Isoprene.nc")
CRIS.set_auto_mask(True)

CRIS.variables

CRIS_iso = CRIS["Isop"][:]
CRIS_lat = CRIS["lat"][:] # get the latitude data
CRIS.close()

CRIS_iso[CRIS_iso < 0] = np.nan 

CRIS_max = np.nanmax(CRIS_iso)/1e+16
CRIS_mean = np.nanmean(CRIS_iso)/1e+16

CRIS_max/CRIS_mean

import matplotlib.pyplot as plt


# Filter out NaN values
CRIS_iso[CRIS_iso < 1.5e+16] = np.nan 
filtered_data = CRIS_iso[~np.isnan(CRIS_iso)]

# Create histogram
hist, bins = np.histogram(filtered_data, bins='auto')

# Plot the histogram
plt.hist(filtered_data, bins='auto')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram without NaNs')
plt.show()

print(bins)
print(hist)


# Get the means per month 
mean_CRIS_iso_subset_NH= calculate_mean_subset_m(23.45, 90, CRIS_lat, CRIS_iso) #tropic of cancer up
mean_CRIS_iso_subset_equator = calculate_mean_subset_m(-23.45, 23.45, CRIS_lat, CRIS_iso) # between tropics of cancer and capricorn
mean_CRIS_iso_subset_SH = calculate_mean_subset_m(-90, -23.45, CRIS_lat, CRIS_iso) #tropic of capricorn down


# Import the Model data
model = netCDF4.Dataset("../results/iso_result_OMI_filterd_mon_mean.nc")
model.set_auto_mask(True)

# Get the isoprene estimate data
m_iso = model["J_times_IspS"][:]
m_lat = model["lat"][:] # get the latitude data
model.close()

m_iso[m_iso > 1000] = np.nan # set filled values to nan

# Get the means per month 
mean_iso_subset_NH= calculate_mean_subset_m(23.45, 90, m_lat, m_iso) #tropic of cancer up
mean_iso_subset_equator = calculate_mean_subset_m(-23.45, 23.45, m_lat, m_iso) # between tropics of cancer and capricorn
mean_iso_subset_SH = calculate_mean_subset_m(-90, -23.45, m_lat, m_iso) #tropic of capricorn down


## Plot the figures
# Create a figure and three subplots arranged vertically
fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# Plot for the first subplot
axs[0].plot(mean_CRIS_iso_subset_NH[:, 0], mean_CRIS_iso_subset_NH[:, 1], 'r-', label='HCHO')
axs[0].plot(mean_iso_subset_NH[:, 0], mean_iso_subset_NH[:, 1], 'b-', label='Model')
axs[0].set_title('Average anomaly in the Northern temperate zone', fontweight='bold')

# Plot for the second subplot
axs[1].plot(mean_CRIS_iso_subset_SH[:, 0], mean_CRIS_iso_subset_SH[:, 1], 'r-', label='HCHO')
axs[1].plot(mean_iso_subset_SH[:, 0], mean_iso_subset_SH[:, 1], 'b-', label='Model')
axs[1].set_title('Average anomaly in the Southern temperate zone', fontweight='bold')

# Plot for the third subplot
axs[2].plot(mean_CRIS_iso_subset_equator[:, 0], mean_CRIS_iso_subset_equator[:, 1], 'r-', label='HCHO')
axs[2].plot(mean_iso_subset_equator[:, 0], mean_iso_subset_equator[:, 1], 'b-', label='Model')
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



