"""
1. Loads and averages the isoprene data for each month
2. Loads the cleaned biome map
3. Filteres the isoprene data by biome and findes the average median scaled value for eah month
4. Plots and saves the results
"""

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import calendar
from scipy.stats import pearsonr
import csv

# Setting matplotlib parameters
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 11

def load_and_average_dataset(path, variable_name):
    """Load a dataset from a NetCDF file and take the average over each year."""
    dataset = netCDF4.Dataset(path)
    dataset.set_auto_mask(False)
    data = dataset[variable_name][:]
    dataset.close()

    # get the monthly means
    data_reshape = data.reshape(3, 12, data.shape[1], data.shape[2])
    data_monthly_averages = np.nanmean(data_reshape, axis=0)

    return data_monthly_averages

def median_scale_data(data):
    """Scale the data by its median value."""
    return data / np.nanmedian(data)

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

    # # Error checking
    # nan_count_before = np.sum(np.isnan(array))
    # nan_count_after = np.sum(np.isnan(filtered_array))
    # print(value,"Number of NaN values before:", nan_count_before)
    # print(value,"Number of NaN values after:", nan_count_after)

    return filtered_array


def apply_filter_and_mean(data, biome):
    median_data = median_scale_data(data) # Is this how I should do it? need to think it through
    result_array = np.empty((12, 8))

    for i in range(12):
        for j in range(1, 9):
            filtered_array = filter_with_nan(median_data[i, :, :], j, biome)
            result_array[i, j-1] = np.nanmean(filtered_array)

    return result_array


# Load the datasets
Model = load_and_average_dataset("clean_data/Model_comp_clean.nc", 'model_data')
OMI = load_and_average_dataset("clean_data/OMI_comp_clean.nc", 'omi_iso')
CrIS = load_and_average_dataset("clean_data/CrIS_comp_clean.nc", 'cris_iso')

# Load the biome map
biome_cdf = netCDF4.Dataset("../isoprene_factor/data/clean_Type3_biome.nc")
biome_cdf.set_auto_mask(False)
biome = biome_cdf['type 3'][:].astype(int)
biome_cdf.close()

# Filter by biome and find the mean per monthly average
model_biome_avg = apply_filter_and_mean(Model, biome)
omi_biome_avg = apply_filter_and_mean(OMI, biome)
cris_biome_avg = apply_filter_and_mean(CrIS, biome)

# Label the biomes
biome_labels = ['Grasses/Cereal', 'Shrubs', 'Broadleaf Crops', 'Savannah', 
                'Evergreen BL', 'Deciduous BL', 'Evergreen NL', 'Deciduous NL']


# Plot the montly isoprene emission patterns per biome
def plot_biome_averages(filename):
    """Plot biome averages."""
    # Create a figure with subplots for each biome
    fig, axs = plt.subplots(2, 4, figsize=(8, 6))
    
    # For x-axis labels
    months = [calendar.month_abbr[i] for i in range(1, 13)]
    # Include December
    x_ticks = list(range(12)[::4]) + [11]
    x_labels = months[::4] + [months[-1]]

    # Loop over the subplots and plot the data
    for i, ax in enumerate(axs.flat):
        # Plot the model data with dashed line
        ax.plot(range(12), model_biome_avg[:, i], 'k-', label=r'{\textit{f}}(T)×J')
        ax.plot(range(12), omi_biome_avg[:, i], 'k--', label=r'OMI')
        ax.plot(range(12), cris_biome_avg[:, i], 'k:', label=r'CrIS')

        # Set the title for each subplot
        ax.set_title(f'{biome_labels[i]}', fontweight='bold',  loc='left')
        
        # Remove top and right spines
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        
        # Label the x axis by months (abbr), every 4th month and December
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_labels)

        # Set y-axis to have 4 ticks and limits
        ax.set_ylim(0, 11)
        
        # Remove y-axis labels for all but the first column and add y-axis label for the first plot of each row
        if i % 4 != 0:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel(r'Median Normalised Isoprene')

    # Create a shared legend below all of the plots
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3)
    
    # Adjust the spacing between subplots
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.4)

    # Save the plot
    plt.savefig(filename)

    # Display the plot
    plt.show()


plot_biome_averages('../final_figures/biome_monthly_averages.png')


# corrolation

def calculate_correlation_by_biome(data1, data2, biome):
    # Create an array to store Pearson correlation coefficients for each biome
    pearson_r_by_biome = np.full((8, 2), np.nan)

    for biome_id in range(1, 9):  # Iterate through each biome (biome IDs from 1 to 8)
        data1_biome_filtered = np.array([])
        data2_biome_filtered = np.array([])

        for month in range(12):  # Iterate through each month (months from 0 to 11)
            data1_month_copy = np.copy(data1[month, :, :])
            data2_month_copy = np.copy(data2[month, :, :])

            # Apply biome filter to set non-biome elements to NaN
            data1_month_copy[np.where(biome != biome_id)] = np.nan

            # Create a mask to keep only valid (non-NaN) elements in both data arrays
            mask = np.logical_and(~np.isnan(data1_month_copy), ~np.isnan(data2_month_copy))

            # Concatenate the valid elements from each month
            data1_biome_filtered = np.concatenate((data1_biome_filtered, data1_month_copy[mask]))
            data2_biome_filtered = np.concatenate((data2_biome_filtered, data2_month_copy[mask]))

        # Calculate Pearson correlation coefficient for the biome and store in the array
        pearson_r_by_biome[biome_id - 1, 0] = pearsonr(data1_biome_filtered, data2_biome_filtered)[0]
        pearson_r_by_biome[biome_id - 1, 1] = pearsonr(data1_biome_filtered, data2_biome_filtered)[1]

    return pearson_r_by_biome

# Calcultate corrolations
biome_corr_model_omi = calculate_correlation_by_biome(Model, OMI, biome)
biome_corr_model_cris = calculate_correlation_by_biome(Model, CrIS, biome)

print("Pearson Correlation Coefficients (Model vs OMI) for each biome:")
print(biome_corr_model_omi)

print("\nPearson Correlation Coefficients (Model vs CrIS) for each biome:")
print(biome_corr_model_cris)


# Create a list of lists for the data
data = [[biome_labels[i], biome_corr_model_omi[i,0], biome_corr_model_cris[i,0]] for i in range(8)]

# Add header row at the top
header = ['biome', 'Model vs OMI', 'Model vs CrIS']

# Write the corr by iome data to a CSV file
with open('../results/biome_correlation_results.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(header)  # Write the header row
    writer.writerows(data)  # Write the data rows

# Made into a table in word and joined with other figure in word (see emails with zelim, subject tBLE)
