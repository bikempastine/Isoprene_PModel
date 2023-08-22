import numpy as np
import netCDF4
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap

# Open the Model data
Model_cdf = netCDF4.Dataset("../data/model_clean.nc")
Model_cdf.set_auto_mask(False)
model = Model_cdf['model_data'][:]
Model_cdf.close()

# Open the OMI data
OMI_cdf = netCDF4.Dataset("../data/OMI_clean.nc")
OMI_cdf.set_auto_mask(False)
omi = OMI_cdf['omi_iso'][:]
OMI_cdf.close()

# filter out when Model is above 1 (I think here the error is driving the high F values)
model[model < 1] = np.nan
omi[model < 1] = np.nan

# calculate the isoprene factor
isoprene_factor = omi/model
isoprene_factor_mean = np.nanmean(isoprene_factor, axis = 0)
isoprene_factor_std = np.nanstd(isoprene_factor, axis = 0)

# enable LaTeX style for the graphs
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 13


def plot_isoprene_factor(isoprene_factors, filename):
    # Extract mean and standard deviation data
    mean_factor, std_factor = isoprene_factors

    # Define custom colormap for mean isoprene factor
    cmap_mean = colors.LinearSegmentedColormap.from_list(
        "", ["lightsteelblue", "navajowhite", "lightsalmon", "coral",
             "salmon", "lightcoral", "indianred", "firebrick"]
    )
    cmap_mean = plt.cm.get_cmap(cmap_mean, 26)

    # Define boundaries for colormap
    bounds_mean = np.linspace(-30, 190, 26)
    norm_mean = colors.BoundaryNorm(bounds_mean, cmap_mean.N)

    # Define custom colormap for standard deviation of isoprene factor
    cmap_std = colors.LinearSegmentedColormap.from_list(
        "", ["paleturquoise", "royalblue", "rebeccapurple", "indigo"]
    )
    cmap_std = plt.cm.get_cmap(cmap_std, 31)

    # Define boundaries for colormap
    bounds_std = np.linspace(0, 310, 31)
    norm_std = colors.BoundaryNorm(bounds_std, cmap_std.N)

    # Create figure with two subplots
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10, 10),
                            subplot_kw={'projection': ccrs.PlateCarree()}, 
                            gridspec_kw={'hspace': 0.20})
    
    # Plot mean isoprene factor
    img_mean = axs[0].imshow(
        mean_factor * 1e6 - np.nanmean(mean_factor * 1e6),
        origin='lower', cmap=cmap_mean, norm=norm_mean,
        extent=[-180, 180, -90, 90]
    )
    axs[0].set_title('Average Isoprene Factor per Gridcell', loc='left')
    axs[0].coastlines()
    axs[0].set_xlabel('Longitude')
    axs[0].set_ylabel('Latitude')

    # Add colorbar for mean isoprene factor
    cbar_mean = fig.colorbar(
        img_mean, ax=axs[0], label='Deviation from Mean Isoprene Factor', 
        ticks=[-30, 0, 30, 60, 90, 120, 150, 180], extend = 'max'
    )
    cbar_mean.ax.tick_params(size=0)
    cbar_mean.ax.set_title(r'×$10^{-6}$', pad=12, loc='left', fontsize = 10)

    # Plot standard deviation of isoprene factor
    img_std = axs[1].imshow(
        std_factor * 1e6, origin='lower', cmap=cmap_std, norm=norm_std, 
        extent=[-180, 180, -90, 90]
    )
    axs[1].set_title(r'Standard Deviation of Isoprene Factor over Time', loc='left')
    axs[1].coastlines()
    axs[1].set_xlabel('Longitude')
    axs[1].set_ylabel('Latitude')

    # Add colorbar for standard deviation of isoprene factor
    cbar_std = fig.colorbar(
        img_std, ax=axs[1], label=r'Standard Deviation of Isoprene Factor', 
        ticks=[0, 50, 100, 150, 200, 250, 300], extend = 'max'
    )
    cbar_std.ax.tick_params(size=0)
    cbar_std.ax.set_title(r'×$10^{-6}$', pad=12, loc='left', fontsize = 10)
    
    plt.subplots_adjust(wspace=0.1)

    # Save the plot
    plt.savefig(filename, dpi=300, bbox_inches='tight')

    plt.show()


# Call function with your data
plot_isoprene_factor([isoprene_factor_mean, isoprene_factor_std], '../../final_figures/isoprene_F_map.png')


# Spatial Descriptive Statistics for each time slice
spatial_variances = np.nanvar(isoprene_factor* 1e6, axis=(1, 2))
spatial_variances.shape

# Temporal Descriptive Statistics for each spatial location
temporal_variances = np.nanvar(isoprene_factor* 1e6, axis=0)
temporal_variances.shape

# Average spatial and temporal variances
average_spatial_variance = np.nanmean(spatial_variances)
average_temporal_variance = np.nanmean(temporal_variances)

print(f"Average Spatial Variance: {average_spatial_variance}")
print(f"Average Temporal Variance: {average_temporal_variance}")

# Determine which variance is greater
if average_spatial_variance > average_temporal_variance:
    print("Data varies more with space.")
elif average_spatial_variance < average_temporal_variance:
    print("Data varies more through time.")
else:
    print("Data varies equally with space and time.")