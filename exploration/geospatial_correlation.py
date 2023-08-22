"""
Per-gridpoint time correlation of OMI, CRIS, MEGAN, and my MODEL

This script is used to:
1. load Model, OMI, CrIS, and MEGAN data from NetCDF files
2. calculate correlations between different pairs of datasets through time (ie for each gridcell)
3. plot and save those correlation maps.
"""

import netCDF4
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import cartopy.crs as ccrs
import matplotlib.colors as colors

def load_dataset(path, variable_name):
    """Load a dataset from a NetCDF file."""
    dataset = netCDF4.Dataset(path)
    dataset.set_auto_mask(False)
    data = dataset[variable_name][:]
    dataset.close()
    return data

# Load the datasets
Model = load_dataset("clean_data/Model_comp_clean.nc", 'model_data')
OMI = load_dataset("clean_data/OMI_comp_clean.nc", 'omi_iso')
CrIS = load_dataset("clean_data/CrIS_comp_clean.nc", 'cris_iso')
MEGAN= load_dataset("clean_data/MEGAN_2005_2014.nc", 'megan_iso')
MEGAN = MEGAN[-36:,:,:]
MEGAN[MEGAN < 1.7e-05] =np.nan

def calculate_corr_time(data1, data2, threshold=30):
    """
    Calculate a correlation array given two datasets.

    Parameters:
        data1 (numpy.ndarray): The first dataset.
        data2 (numpy.ndarray): The second dataset.
        threshold (int): The minimum number of non-NaN values required to compute the correlation.

    Returns:
        numpy.ndarray: A 2D array of correlation coefficients.
    """
    corr_array = np.zeros(data1.shape[1:])

    for lat in range(data1.shape[1]):
        for lon in range(data1.shape[2]):
            # get all of the time points in that gridcell
            data_A = np.array(data1[:, lat, lon])
            data_B = np.array(data2[:, lat, lon])

            # note the NANs as if there is one nan, that time point is discarded
            mask = np.logical_and(~np.isnan(data_A), ~np.isnan(data_B))

            # If there are fewer than threshold non-NaN values, set the correlation coefficient to NaN
            if np.count_nonzero(mask) < threshold:
                corr_array[lat, lon] = np.nan
            else:
                # Otherwise, compute the correlation coefficient
                masked_A = data_A[mask]
                masked_B = data_B[mask]
                corr_array[lat, lon] = pearsonr(masked_A, masked_B)[0]
    
    return corr_array

# Calculate correlations
correlations = {
    'model_cris': calculate_corr_time(Model, CrIS),
    'model_OMI': calculate_corr_time(Model, OMI),
    'megan_cris': calculate_corr_time(MEGAN, CrIS),
    'megan_OMI': calculate_corr_time(MEGAN, OMI),
    'cris_omi': calculate_corr_time(OMI, CrIS),
}

# enable LaTeX style for the graphs
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 13

def plot_correlations(correlations, title1, title2, filename):
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7.5,8),
                            subplot_kw={'projection': ccrs.PlateCarree()})

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.12)

    # Define the custom color map
    cmap = colors.LinearSegmentedColormap.from_list("", ["royalblue","navajowhite", "firebrick"])
    
    # Create a custom color map with discrete colors
    cmap = plt.cm.get_cmap(cmap, 9) 

    # Define the boundaries for your colors
    bounds = np.linspace(-1, 1, 9)

    # Create a normalization object using these boundaries
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # Plot the first correlation
    img1 = axs[0].imshow(correlations[0], origin='lower', cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    axs[0].set_title(r'{\textit{f}}(T)×J vs OMI', loc='left')
    axs[0].coastlines()

    # Plot the second correlation
    img2 = axs[1].imshow(correlations[1], origin='lower', cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    axs[1].set_title(r'{\textit{f}}(T)×J vs CrIS', loc='left')
    axs[1].coastlines()

    # Add a colorbar that is shared between the two subplots
    cbar = fig.colorbar(img1, ticks=bounds, ax=axs.ravel().tolist(), orientation='horizontal', fraction=0.04, pad=0.04)
    cbar.set_label('Pearson correlation coefficient (r)')

    # Save the figure
    plt.savefig(filename, dpi=300, bbox_inches='tight')

    # Show the plot
    plt.show()

# Plotting correlations
main_title = "Evaluating Modelled Isoprene: Pearson's R with OMI and CrIS Satellite Products (2012-2014)"
plot_correlations([correlations['model_OMI'], correlations['model_cris']],'f(T)×J vs OMI' , 'f(T)×J vs CrIS',  '../final_figures/spatial_model_correlations.png')










####################
##P value #######

def calculate_corr_time_pvalue(data1, data2, threshold=30):
    """
    Calculate a correlation array given two datasets.

    Parameters:
        data1 (numpy.ndarray): The first dataset.
        data2 (numpy.ndarray): The second dataset.
        threshold (int): The minimum number of non-NaN values required to compute the correlation.

    Returns:
        numpy.ndarray: A 2D array of correlation coefficients.
    """
    corr_array = np.zeros(data1.shape[1:])

    for lat in range(data1.shape[1]):
        for lon in range(data1.shape[2]):
            # get all of the time points in that gridcell
            data_A = np.array(data1[:, lat, lon])
            data_B = np.array(data2[:, lat, lon])

            # note the NANs as if there is one nan, that time point is discarded
            mask = np.logical_and(~np.isnan(data_A), ~np.isnan(data_B))

            # If there are fewer than threshold non-NaN values, set the correlation coefficient to NaN
            if np.count_nonzero(mask) < threshold:
                corr_array[lat, lon] = np.nan
            else:
                # Otherwise, compute the correlation coefficient
                masked_A = data_A[mask]
                masked_B = data_B[mask]
                corr_array[lat, lon] = pearsonr(masked_A, masked_B)[1]
    
    return corr_array



# Calculate correlations
correlations_p_value = {
    'model_cris': calculate_corr_time_pvalue(Model, CrIS),
    'model_OMI': calculate_corr_time_pvalue(Model, OMI),
    'megan_cris': calculate_corr_time_pvalue(MEGAN, CrIS),
    'megan_OMI': calculate_corr_time_pvalue(MEGAN, OMI),
    'cris_omi': calculate_corr_time_pvalue(OMI, CrIS),
}

# enable LaTeX style for the graphs
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 11

def plot_correlations(correlations, title1, title2, filename):
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7.5, 8),
                            subplot_kw={'projection': ccrs.PlateCarree()})

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.12)

    # Define the custom colors
    color_list = ['darkgreen', 'green', 'limegreen', 
                  'lightsalmon', '#FC9272', '#FB6A4A',
                  '#EF3B2C', '#CB181D', 'firebrick']
        
        #     FFA07A (lightsalmon)
        # FC9272
        # FB6A4A
        # EF3B2C
        # CB181D
        # 99000D
    
    # Define the boundaries for your colors
    bounds = [0, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1]

    # Create a colormap from your list of colors
    cmap = colors.ListedColormap(color_list)

    # Create a normalization object using these boundaries
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # Plot the first correlation
    img1 = axs[0].imshow(correlations[0], origin='lower', cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    axs[0].set_title(title1, loc='left')
    axs[0].coastlines()

    # Plot the second correlation
    img2 = axs[1].imshow(correlations[1], origin='lower', cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    axs[1].set_title(title2, loc='left')
    axs[1].coastlines()

    # Add a colorbar that is shared between the two subplots
    cbar = fig.colorbar(img1, ticks=bounds, ax=axs.ravel().tolist(), orientation='horizontal', fraction=0.04, pad=0.04)
    cbar.set_label('Pearson correlation p-value')

    # Save the figure (optional)
    plt.savefig(filename)

    # Show the plot
    plt.show()

# Plotting correlations
main_title = "Evaluating Modelled Isoprene: Pearson's R with OMI and CrIS Satellite Products (2012-2014)"
plot_correlations([correlations_p_value['model_OMI'], correlations_p_value['model_cris']],'f(T)×J vs OMI' , 'f(T)×J vs CrIS',  '../final_figures/appendix/spatial_model_correlations_p_val.png')







