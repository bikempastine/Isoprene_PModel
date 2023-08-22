"""
1. Finds the pearson corrolation coefficient between each months for each set of data
2. plots it in a beautiful line graph
"""

import calendar
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs

# Global parameters
DATA_PATH = "clean_data/"
MODEL_FILE = "Model_comp_clean.nc"
OMI_FILE = "OMI_comp_clean.nc"
CRIS_FILE = "CrIS_comp_clean.nc"
MEGAN_FILE = "MEGAN_2005_2014.nc"
MODEL_VAR = 'model_data'
OMI_VAR = 'omi_iso'
CRIS_VAR = 'cris_iso'
MEGAN_VAR = 'megan_iso'
LATITUDE_VAR = 'lat'

# Setting matplotlib parameters
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 13

def open_and_read_ncdf(file_path, data_name):
    """Open a netCDF file, retrieve and return the required data and latitudes."""
    with netCDF4.Dataset(file_path) as cdf:
        cdf.set_auto_mask(False)
        data = cdf[data_name][:]
        lat = cdf[LATITUDE_VAR][:]
    return data, lat

def find_correlation(data1, data2, month):
    """
    Find the correlation coefficient between two datasets for a given month.
    
    Args:
        data (ndarray): Dataset array.
        month (int): Month index.
        
    Returns:
        float: Correlation coefficient.
    """
    valid_mask = np.logical_and(~np.isnan(data1[month,:,:]), ~np.isnan(data2[month,:,:]))
    x_valid = data1[month,:,:][valid_mask]
    y_valid = data2[month,:,:][valid_mask]
    corr_matrix = np.corrcoef(x_valid.flatten(), y_valid.flatten())

    return corr_matrix[0, 1]


# Open and clean the Model, OMI, MEGAN, and CrIS datasets
model_data, model_lat = open_and_read_ncdf(DATA_PATH + MODEL_FILE, MODEL_VAR)
omi_data, omi_lat = open_and_read_ncdf(DATA_PATH + OMI_FILE, OMI_VAR)
cris_data, cris_lat = open_and_read_ncdf(DATA_PATH + CRIS_FILE, CRIS_VAR)
megan_data, megan_lat = open_and_read_ncdf(DATA_PATH + MEGAN_FILE, MEGAN_VAR)
megan_data = megan_data[-36:, :, :]

# Calculate the correlations
correlations = {label: np.zeros(36) for label in ['omi_model', 'cris_model', 'cris_megan', 'omi_megan', 'cris_omi']}
for i in range(1, 36):
    correlations['omi_model'][i] = np.square(find_correlation(omi_data, model_data, i))
    correlations['cris_model'][i] = np.square(find_correlation(cris_data, model_data, i))
    correlations['cris_megan'][i] = np.square(find_correlation(cris_data, megan_data, i))
    correlations['omi_megan'][i] = np.square(find_correlation(omi_data, megan_data, i))
    correlations['cris_omi'][i] = np.square(find_correlation(cris_data, omi_data, i))



def plot_correlations(filename):
    """Plot correlation comparisons"""
    plt.figure(figsize=(8.5, 5.5))

    x = np.arange(0, 37)
    months = [calendar.month_abbr[(i % 12) + 1] + ' ' + str(2012 + (i // 12)) for i in range(37)]

    # Plot correlations
    plt.plot(x[1:36], correlations['omi_model'][1:], label= r'{\textit{f}}(T)×J vs OMI', color='firebrick', linestyle='--', linewidth=1.5)
    plt.plot(x[1:36], correlations['omi_megan'][1:], label= r'MEGAN vs OMI', color='black', linestyle='--', linewidth=1.5)
    plt.plot(x[1:36], correlations['cris_model'][1:], label= r'{\textit{f}}(T)×J vs CrIS', color='firebrick', linewidth=1.5)
    plt.plot(x[1:36], correlations['cris_megan'][1:], label= r'MEGAN vs CrIS', color='black', linewidth=1.5)
    
    # Set labels
    plt.ylabel(r'Coefficient of determination ($R^2$)')
    
    # Set y axis limit
    plt.ylim([0, 1])
    
    # Set xticks
    ax = plt.gca()
    ax.set_xticks(x[:])
    ax.set_xticklabels(['' if i % 6 != 0 else month for i, month in enumerate(months[:])], rotation=0, ha='center')

    for i, tick in enumerate(ax.xaxis.get_major_ticks()):
        if i % 6 == 0:
            tick.tick1line.set_markeredgewidth(1.6)
        else:
            tick.tick1line.set_markeredgewidth(0.8)

    # Remove top and right spines
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    # Add legend below the plot
    plt.legend(ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.1))  # Adjust the 'bbox_to_anchor' as needed to position the legend
    
    # Save the plot
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(filename, dpi=300, bbox_inches='tight')

    # Display the plot
    plt.show()




# Call function to plot correlations
plot_correlations('../final_figures/corrolation_line_plot.png')


def find_corr(data1, data2):
    x_valid = []
    y_valid = []
    for month in range(36):
        valid_mask = np.logical_and(~np.isnan(data1[month,:,:]), ~np.isnan(data2[month,:,:]))
        x_valid = np.concatenate([x_valid, data1[month,:,:][valid_mask]])
        y_valid =  np.concatenate([y_valid, data2[month,:,:][valid_mask]])
    
    corr_matrix = np.corrcoef(x_valid, y_valid)
    return corr_matrix[0,1]

find_corr(omi_data, model_data)
find_corr(cris_data, model_data)

find_corr(omi_data, megan_data)
find_corr(cris_data, megan_data)

find_corr(cris_data, omi_data)