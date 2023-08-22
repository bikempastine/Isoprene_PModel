"""
This script compares the OMI, CrIS, and the Isoprene Model datasets from 2012 to 2014.
The steps include:
1. Opening the data and cleaning it.
2. Applying a median feature scaling function to make data comparable.
3. Visualising the maps of the isoprene estimates for averaged NH and SH summer.
4. Visualising averaged monthly isoprene model estimates (for the appendix)
"""

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import colors
import calendar
import cartopy.crs as ccrs

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 13

def median_scale_data(data):
    """Scale the data by its median value."""
    return data / np.nanmedian(data)


def open_ncdf(file_path, data_name):
    """Open a netCDF file, retrieve and return the required data and latitudes."""
    cdf = netCDF4.Dataset(file_path)
    cdf.set_auto_mask(False)
    data = cdf[data_name][:]
    lat = cdf['lat'][:]
    cdf.close()
    return data, lat

def calculate_seasonal_averages(data):
    """
    Assumes data is a 3D array with shape (months, lat, lon) 
    and months is a multiple of 3 (as each season consists of 3 months).
    """
    # Reshape the data so that each season is grouped together
    reshaped_1 = data.reshape(3, 12, data.shape[1], data.shape[2])

    # Calculate the average over the seasonal axis
    yearly_averages = np.nanmean(reshaped_1, axis=0)

    reshaped_2 = yearly_averages.reshape((4, 3, 360, 720))

    seasonal_averages = np.nanmean(reshaped_2, axis=1)

    return seasonal_averages


# Open and scale the Model, OMI, and CrIS datasets
MODEL, Model_lat = open_ncdf("clean_data/Model_comp_clean.nc", 'model_data')
OMI, OMI_lat = open_ncdf("clean_data/OMI_comp_clean.nc", 'omi_iso')
CRIS, CrIS_lat = open_ncdf("clean_data/CrIS_comp_clean.nc", 'cris_iso')
MODEL_scaled, OMI_scaled, CRIS_scaled = map(median_scale_data, [MODEL, OMI, CRIS])

# Calculate seasonal averages
OMI_seasonal = calculate_seasonal_averages(OMI_scaled)
CRIS_seasonal = calculate_seasonal_averages(CRIS_scaled)
MODEL_seasonal = calculate_seasonal_averages(MODEL_scaled)

def plot_winter_summer(filename):
    data_sets = [MODEL_seasonal, OMI_seasonal, CRIS_seasonal]
    data_sets_names = ['Model', 'OMI', 'CRIS']
    seasons = ['December to February', 'June - August']

    # Create a custom color map with discrete colors
    cmap = colors.LinearSegmentedColormap.from_list("", ["paleturquoise","royalblue","rebeccapurple", "indigo"])

    # Define the boundaries for your colors
    bounds = np.linspace(0, 25, 9)

    # Create a normalization object using these boundaries
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(11, 8),
                            subplot_kw={'projection': ccrs.PlateCarree()})

    # Setting the data set names as titles for each column    
    axs[0,0].set_title(r'SH Summer (December - February)', loc='center')
    axs[0,1].set_title(r'NH Summer (June - August)', loc='center')

    # Setting the season names as labels for each row
    # for ax, row in zip(axs[:,0], data_sets_names):
    #     fig.text(0.13, 0.33, row, ha='center', va='center', rotation='vertical')

    fig.text(0.14, 0.77, r'{\textit{f}}(T)×J', ha='center', va='center', rotation='vertical')
    fig.text(0.14, 0.52, r'OMI', ha='center', va='center', rotation='vertical')
    fig.text(0.14, 0.29, r'CrIS', ha='center', va='center', rotation='vertical')

    for i, data_set in enumerate(data_sets):
        for j in range(2):
            img = axs[i, j].imshow(data_set[j], origin='lower', extent=[-180, 180, -90, 90], cmap=cmap, norm=norm)
            axs[i, j].coastlines()

    fig.subplots_adjust(wspace=0.01, hspace=0.04)

    cbar = fig.colorbar(img, ticks= [0, 5, 10, 15, 20, 25], ax=axs.ravel().tolist(), orientation='horizontal', extend = 'max', fraction=0.04, pad=0.04)
    cbar.set_label('median scaled isoprene')

    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

plot_winter_summer('../final_figures/normalised_isoprene.png')


# Average for each month
model_reshape = MODEL.reshape(3, 12, MODEL.shape[1], MODEL.shape[2])
model_averages = np.nanmean(model_reshape, axis=0)

def plot_all_model(filename):
    # Create a custom color map with discrete colors
    cmap = colors.LinearSegmentedColormap.from_list("", ["paleturquoise","royalblue","rebeccapurple", "indigo"])

    # Define the boundaries for your colors
    bounds = np.linspace(0, 37, 9)

    # Create a normalization object using these boundaries
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(11, 9),
                            subplot_kw={'projection': ccrs.PlateCarree()})

    months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

    for i in range(12):
        row = i // 4  # Calculate row index
        col = i % 4   # Calculate column index
        img = axs[row, col].imshow(model_averages[i], origin='lower', extent=[-180, 180, -90, 90], cmap=cmap, norm=norm)
        axs[row, col].coastlines()
        axs[row, col].set_title(f'{months[i]}',loc='left')  # Setting the month names as titles for each plot

    fig.subplots_adjust(wspace=0.01, hspace=0.01, bottom=0.15)
    
    cbar = fig.colorbar(img, ticks= bounds, ax=axs.ravel().tolist(), orientation='horizontal',fraction=0.03, pad =0.01)
    cbar.set_label(r'{\textit{f}}(T)×J')

    plt.savefig(filename)
    plt.show()

plot_all_model('../final_figures/model_averaged_plot.png')





def plot_winter_summer(filename):
    data_sets = [MODEL_seasonal, OMI_seasonal, CRIS_seasonal]
    data_sets_names = ['Model', 'OMI', 'CRIS']

    cmap = colors.LinearSegmentedColormap.from_list("", ["paleturquoise","royalblue","rebeccapurple", "indigo"])
    bounds = np.linspace(0, 25, 9)
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(11, 8),
                            subplot_kw={'projection': ccrs.PlateCarree()})
    
    axs[0,0].set_title(r'December - February', loc='left')
    axs[0,1].set_title(r'June - August', loc='left')

    fig.text(0.14, 0.77, r'{\textit{f}}(T)×J', ha='center', va='center', rotation='vertical')
    fig.text(0.14, 0.52, r'OMI', ha='center', va='center', rotation='vertical')
    fig.text(0.14, 0.29, r'CrIS', ha='center', va='center', rotation='vertical')

    # Define extent for South America and the specified region in the USA
    south_america_extent = [-90, -30, -60, 15]
    usa_extent = [-107.5, -65, 17, 42]  # Puerto Rico to New Mexico, and up to New York

    for i, data_set in enumerate(data_sets):
        # First column - South America
        axs[i, 0].set_extent(south_america_extent)
        axs[i, 0].imshow(data_set[0], origin='lower', extent=[-180, 180, -90, 90], cmap=cmap, norm=norm)
        axs[i, 0].coastlines()

        # Second column - Specified USA region
        axs[i, 1].set_extent(usa_extent)
        axs[i, 1].imshow(data_set[1], origin='lower', extent=[-180, 180, -90, 90], cmap=cmap, norm=norm)
        axs[i, 1].coastlines()

    fig.subplots_adjust(wspace=0.01, hspace=0.04)

    # cbar = fig.colorbar(fig, ticks= [0, 5, 10, 15, 20, 25], ax=axs.ravel().tolist(), orientation='horizontal', extend = 'max', fraction=0.04, pad=0.04)
    # cbar.set_label('median normalised isoprene')

    #plt.savefig(filename)
    plt.show()

plot_winter_summer('../final_figures/normalised_isoprene.png')
