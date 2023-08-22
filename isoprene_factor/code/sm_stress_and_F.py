import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# enable LaTeX style for the graphs
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 11

# open the data
isoprene_factor_data = np.load('../data/isoprene_factor_2005_2014.npy')
soil_moisture_data = np.load('../data/sm_stress_2005_2014.npy')


def mask_values(sm_data, array2, stress_condition):
    """
    Create masked arrays where NaN and specific stress condition values are masked.

    Parameters:
        sm_data (np.ndarray): soil moisture.
        array2 (np.ndarray): Second input array.
        stress_condition (str): 'y' for stressed conditions, 'n' for non-stressed conditions.

    Returns:
        np.ndarray, np.ndarray: Masked arrays for array1 and array2 based on the stress condition.
    """
    array2_valid = []

    for month in range(120):
        
        if stress_condition == 'y' :
            sm_mask = (sm_data[month,:,:] < 1)
        elif stress_condition == 'n' :
            sm_mask = (sm_data[month,:,:] == 1)

        valid_mask = np.logical_and(sm_mask, ~np.isnan(array2[month,:,:]))
        array2_valid = np.concatenate([array2_valid, array2[month,:,:][valid_mask]])

    return array2_valid

# Mask the data for stressed and non-stressed conditions
isoprene_factor_stressed = mask_values(soil_moisture_data, isoprene_factor_data, 'y')
isoprene_factor_non_stressed = mask_values(soil_moisture_data, isoprene_factor_data, 'n')

# Perform Mann-Whitney U test
U, p = stats.mannwhitneyu(isoprene_factor_stressed, isoprene_factor_non_stressed)

# Calculate median and mean values
np.median(isoprene_factor_stressed)
np.median(isoprene_factor_non_stressed)

np.mean(isoprene_factor_stressed)
np.mean(isoprene_factor_non_stressed)

# Create a boxplot

# Mask the data for isoprene factor, OMI, and model for stressed and non-stressed conditions
isoprene_factor_stressed = mask_values(soil_moisture_data, isoprene_factor_data * 1e6, 'y')
isoprene_factor_non_stressed = mask_values(soil_moisture_data, isoprene_factor_data * 1e6, 'n')

# Create subplots to show the three boxplots side by side
def boxplot_omi_model_F_sm_stress(filename):

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    # Plot isoprene factor boxplot
    axs[0].boxplot([isoprene_factor_stressed, isoprene_factor_non_stressed], showfliers=False)
    axs[0].set_xticklabels(['Soil Moisture Stressed', 'Soil Moisture Unstressed'])
    axs[0].set_ylabel(r'Isoprene Factor ×$10^{-6}$')
    axs[0].set_title("Boxplot of Isoprene Factor Grouped by Soil Moisture Stress")

    # Plot OMI boxplot
    axs[1].boxplot([omi_stressed, omi_non_stressed], showfliers=False)
    axs[1].set_xticklabels(['Soil Moisture Stressed', 'Soil Moisture Unstressed'])
    axs[1].set_ylabel(r'OMI ×$10^{-5}$')
    axs[1].set_title("Boxplot of OMI Grouped by Soil Moisture Stress")

    # Plot model boxplot
    axs[2].boxplot([model_stressed, model_non_stressed], showfliers=False)
    axs[2].set_xticklabels(['Soil Moisture Stressed', 'Soil Moisture Unstressed'])
    axs[2].set_ylabel('Model')
    axs[2].set_title("Boxplot of Model Grouped by Soil Moisture Stress")

    plt.tight_layout()

    plt.savefig(filename)
    plt.show()


boxplot_omi_model_F_sm_stress('../figures/sm_boxplots.png')





#####################

def part_values(sm_data, array2, stress_condition):
    """
    Create masked arrays where NaN and specific stress condition values are masked.

    Parameters:
        sm_data (np.ndarray): soil moisture.
        array2 (np.ndarray): Second input array.
        stress_condition (str): 'y' for stressed conditions, 'n' for non-stressed conditions.

    Returns:
        np.ndarray, np.ndarray: Masked arrays for array1 and array2 based on the stress condition.
    """
    array2_valid = []

    for month in range(120):
        
        if stress_condition == 'n' :
            sm_mask = (sm_data[month,:,:] == 1)
        if stress_condition == 'y' :
            sm_mask = (sm_data[month,:,:] < 1)
        if stress_condition == 0.8 :
            sm_mask = (np.logical_and(sm_data[month, :, :] >= 0.8 , sm_data[month, :, :] < 1))
        if stress_condition == 0.6 :
            sm_mask = (np.logical_and(sm_data[month, :, :] >= 0.6 , sm_data[month, :, :] < 0.8))
        if stress_condition == 0.4 :
            sm_mask = (np.logical_and(sm_data[month, :, :] >= 0.4 , sm_data[month, :, :] < 0.6))
        if stress_condition == 0.2 :
            sm_mask = (sm_data[month,:,:] < 0.4)

        valid_mask = np.logical_and(sm_mask, ~np.isnan(array2[month,:,:]))
        array2_valid = np.concatenate([array2_valid, array2[month,:,:][valid_mask]])

    return array2_valid

# Mask the data for stressed and non-stressed conditions
stress_values = ['n','y', 0.8, 0.6, 0.4, 0.2]
isoprene_factor_parted = list(range(len(stress_values))) 

for i in range(len(stress_values)):
    isoprene_factor_parted[i] = part_values(soil_moisture_data, isoprene_factor_data* 1e6, stress_condition=stress_values[i])


# Create the box plot
def box_plot_sm_parted_F():

    # Define line properties
    medianprops = {'color': 'black', 'linewidth': 1.5}
    boxprops = {'color': 'black', 'linewidth': 1.5}
    whiskerprops = {'color': 'black', 'linewidth': 1}
    capprops = {'color': 'black', 'linewidth': 1}
    meanpointprops = dict(marker='.', markeredgecolor='black',
                      markerfacecolor='firebrick')

    # Initialise the figure
    plt.figure(figsize=(5.5, 4))
    bx = plt.boxplot(isoprene_factor_parted, vert=False, showfliers=False, showmeans=True, patch_artist=True,
                meanprops=meanpointprops,
                medianprops=medianprops, boxprops=boxprops, 
                whiskerprops=whiskerprops, capprops=capprops)  
    
    for box in bx['boxes'][0:2]:
        box.set(facecolor = 'lightgrey' )

    for box in bx['boxes'][2:7]:
        box.set(facecolor = 'white' )

    plt.xlabel(r'Isoprene Factor ×$10^{-6}$',  loc='right')
    plt.yticks(range(1,7), [r'$\beta$ = 1', r'$\beta$ $<$ 1', r'1 $>$ $\beta$ $\geq$ 0.8', 
                            r'0.8 $>$ $\beta$ $\geq$ 0.6',r'0.6 $>$ $\beta$ $\geq$ 0.4', r'$\beta$ $<$ 0.4'] )
    
    
    #plt.ylabel('Soil Moisture Stress')
    plt.title(r'Isoprene Factor by Soil Moisture Stress ($\beta$)', fontweight='bold',  loc='left' )
    plt.subplots_adjust(left=0.25)
    
    plt.show()
    

box_plot_sm_parted_F()



H, p = kruskal(*isoprene_factor_parted) #star used to unpack
print('H statistic:', H)
print('p value:', p) 


np.nanmedian()




plt.hist(soil_moisture_data.flatten(), bins=500, edgecolor='black')
# Add labels and title
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram')

# Show the histogram
plt.show()