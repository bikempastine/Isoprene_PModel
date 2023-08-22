
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
import calendar
from scipy.stats import pearsonr
import csv


# Setting matplotlib parameters
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 12

def calculate_mean_subset_m(latitude_min, latitude_max, latitude, data):
    """
    Calculate the mean anomaly per month for the provided latitude range
    """
    lat_indices = np.where((latitude >= latitude_min) & (latitude <= latitude_max))[0]
    data_subset = data[:, lat_indices, :]
    data_mean_subset = np.nanmean(data_subset, axis=0)
    anomaly_subset = np.subtract(data_subset, data_mean_subset)
    norm_subset = np.divide(anomaly_subset, data_mean_subset)

    norm_subset_reshape = norm_subset.reshape(3, 12, norm_subset.shape[1], norm_subset.shape[2])
    norm_subset_monthly_averages = np.nanmean(norm_subset_reshape, axis=0)
   
    mean_data_subset = np.column_stack((np.arange(12), [np.nanmean(norm_subset_monthly_averages[i, :, :]) for i in range(12)]))
    return mean_data_subset


def open_ncdf(file_path, data_name):
    """
    Open a netCDF file, retrieve and return the required data and latitudes.
    """
    with netCDF4.Dataset(file_path) as cdf:
        cdf.set_auto_mask(False)
        data = cdf[data_name][:]
        lat = cdf['lat'][:]
    return data, lat


def get_mean_lat_subsets(filepath, dataname):
    """
    Open netCDF file and calculate mean subset data.
    """
    data, lat = open_ncdf(filepath, dataname)
    data_nh = calculate_mean_subset_m(23.45, 90, lat, data)
    data_equator = calculate_mean_subset_m(-23.45, 23.45, lat, data)
    data_sh = calculate_mean_subset_m(-90, -23.45, lat, data)

    return data_nh, data_equator, data_sh

def plot_subplots(data_list, label_list, title_list, filename):
    """
    This function plots multiple subplots, based on input data.
    """

    # Create a figure and subplot structure
    # This creates 3 subplots in vertical arrangement (3 rows, 1 column)
    fig, axs = plt.subplots(1,3, figsize=(10, 4), sharex=True)

    # Loop over each subplot and plot data
    for idx, ax in enumerate(axs):

        # Extract data for this subplot
        data1 = data_list[0][idx]
        data2 = data_list[1][idx]
        data3 = data_list[2][idx]

        # Add a horizontal line at y=0
        ax.axhline(0, color='whitesmoke')

        # Plot each data set on the same subplot
        ax.plot(data1[:, 0], data1[:, 1], 'k-', label=r'{\textit{f}}(T)×J')
        ax.plot(data2[:, 0], data2[:, 1], 'k--', label = f'{label_list[1]}')
        ax.plot(data3[:, 0], data3[:, 1], 'k:', label = f'{label_list[2]}')

        # Set title for this subplot
        ax.set_title(title_list[idx])

        # Set x-axis ticks as months and rotate labels for readability
        ax.set_xticks(range(12))
        ax.set_xticklabels(calendar.month_abbr[1:], rotation=90, ha='center')

        # set y axis limit so its the same for each one
        ax.set_ylim(-1, 2.5)

        # only show y tick labels and main label on first subplot
        if idx != 0:
            ax.set_yticklabels([])
       
    axs[0].set_ylabel('Average Anomaly')     

    # Adjust subplot layout 
    plt.subplots_adjust(hspace=0.4, bottom=0.2)  # Make space for the legend at the bottom

    # Add a ledgend
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3)

    # Save the figure
    plt.savefig(filename, dpi=300, bbox_inches='tight')

    # Display the figure with subplots
    plt.show()


model_data = get_mean_lat_subsets("clean_data/Model_comp_clean.nc", 'model_data')
omi_data = get_mean_lat_subsets("clean_data/OMI_comp_clean.nc", 'omi_iso')
cris_data = get_mean_lat_subsets("clean_data/CrIS_comp_clean.nc", 'cris_iso')

plot_subplots([model_data, omi_data, cris_data],
              ['f(T)×J', 'OMI', 'CrIS'],
              ['Northern Temperate Zone',
               'Tropical Zone ',
               'Southern Temperate Zone'],
               '../final_figures/lat_bound_anomoly_line_graph.png')



# corrolations
def calcualte_corr_by_lat(data1_list, data2_list):
    
    data1, data1_lat = open_ncdf(data1_list[0], data1_list[1])
    data2, data2_lat= open_ncdf(data2_list[0], data2_list[1])

    def subset_by_lat(latitude_min, latitude_max, data1, latitude1, data2, latitude2):
        lat_indices1 = np.where((latitude1 >= latitude_min) & (latitude1 <= latitude_max))[0]
        data1_subset = data1[:, lat_indices1, :]

        lat_indices2 = np.where((latitude2 >= latitude_min) & (latitude2 <= latitude_max))[0]
        data2_subset = data2[:, lat_indices2, :]

        return [data1_subset, data2_subset]
    
    data_nh = subset_by_lat(23.45, 90, data1, data1_lat,data2, data2_lat)
    data_equator = subset_by_lat(-23.45, 23.45, data1, data1_lat,data2, data2_lat)
    data_sh = subset_by_lat(-90, -23.45, data1, data1_lat,data2, data2_lat)

    def find_corr(data_list):
        A_filtered = np.array([])
        B_filtered = np.array([])
        
        for month in range(12):
            data_A = np.copy(data_list[0][month,:,:])
            data_B = np.copy(data_list[1][month,:,:])

            mask = np.logical_and(~np.isnan(data_A), ~np.isnan(data_B))

            # Concatenate the valid elements from each month
            A_filtered = np.concatenate((A_filtered, data_A[mask]))
            B_filtered = np.concatenate((B_filtered, data_B[mask]))

        corr, p = pearsonr(A_filtered, B_filtered)
        return corr
    
    corr_nh = find_corr(data_nh)
    corr_eq = find_corr(data_equator)
    corr_sh = find_corr(data_sh)

    zone_names = ['Northern Temperate Zone','Tropical Zone ','Southern Temperate Zone']
    correlations = [corr_nh, corr_eq, corr_sh]

    zone_corr = list(zip(zone_names, correlations))

    return zone_corr



Model_list = ["clean_data/Model_comp_clean.nc", 'model_data']
OMI_list = ["clean_data/OMI_comp_clean.nc", 'omi_iso']
CrIS_list = ["clean_data/CrIS_comp_clean.nc", 'cris_iso']

# Calculate correlations for Model vs OMI and Model vs CrIS
corr_model_omi = calcualte_corr_by_lat(Model_list, OMI_list)
corr_model_cris = calcualte_corr_by_lat(Model_list, CrIS_list)

# Save the results in a CSV file
output_file = "../results/lat_bound_correlation_results.csv"

# Prepare data to write into the CSV file
data_to_write = [('zone', 'Model vs OMI', 'Model vs CrIS')]
zone_names = ['Northern Temperate Zone','Tropical Zone ','Southern Temperate Zone']
for zone_name, corr_omi, corr_cris in zip(zone_names, corr_model_omi, corr_model_cris):
    data_to_write.append((zone_name, corr_omi[1], corr_cris[1]))

# Write data to CSV
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(data_to_write)

