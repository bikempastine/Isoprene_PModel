
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kruskal
import scikit_posthocs as sp

# enable LaTeX style for the graphs
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 11

# Load the biome map
biome_cdf = netCDF4.Dataset("../data/clean_Type3_biome.nc")
biome_cdf.set_auto_mask(False)
biome = biome_cdf['type 3'][:]
biome_cdf.close()

isoprene_factor = np.load('../data/isoprene_factor_2005_2014.npy')

# exploration

# Get the unique values and their counts

biome_labels = [
    "Water Bodies",
    "Grasses/Cereal",
    "Shrublands",
    "Broadleaf Crops",
    "Savannas",
    "Evergreen BL",
    "Deciduous BL",
    "Evergreen NL",
    "Deciduous NL",
    "Unvegetated",
    "Urban and Built-up Lands"
]

# List of lists to hold the isoprene_factor values for each biome
data = [[] for _ in range(11)]  # assuming 11 unique biomes

# Loop through each time point
for t in range(isoprene_factor.shape[0]):
    # Loop through each unique biome
    for i in range(11):  # assuming 11 unique biomes
        # Create a copy of the isoprene_factor for the current time point
        isoprene_factor_t = isoprene_factor[t].copy() * 1e6
        # Set to NaN where the biome isn't equal to the current biome
        isoprene_factor_t[biome != i] = np.nan
        # Flatten the array and filter out the NaNs
        isoprene_factor_t = isoprene_factor_t[~np.isnan(isoprene_factor_t)]
        # Append the remaining values to the corresponding list
        data[i].extend(isoprene_factor_t.tolist())


def plot_biome_bound_f(filename):
    # Slice the data and biome_labels lists to include only biomes 1 to 9
    data_to_plot = data[1:9]
    biome_labels_to_plot = biome_labels[1:9]

    # Define line properties
    medianprops = {'color': 'black', 'linewidth': 1.5}
    boxprops = {'color': 'black', 'linewidth': 1.5}
    whiskerprops = {'color': 'black', 'linewidth': 1}
    capprops = {'color': 'black', 'linewidth': 1}
    meanpointprops = dict(marker='.', markeredgecolor='black',
                      markerfacecolor='firebrick')

    # Initialise the figure
    plt.figure(figsize=(5.5, 4))
    plt.boxplot(data_to_plot, vert=False, showfliers=False, showmeans=True,
                meanprops=meanpointprops,
                medianprops=medianprops, boxprops=boxprops, 
                whiskerprops=whiskerprops, capprops=capprops)  
    
    plt.yticks(range(1, 9), biome_labels_to_plot)  # range should match the number of biomes
    plt.title(f' Isoprene Factor by Biome', fontweight='bold',  loc='left' )
    plt.xlabel(r'Isoprene Factor ×$10^{-6}$',  loc='right')

    # Adjust the subplots to make sure the y labels are not cut off
    plt.subplots_adjust(left=0.25)

    # Save the figure
    plt.savefig(filename)
    
    plt.show()

plot_biome_bound_f('../../final_figures/biome_bound_f_box_and_whickers.png')


# Prepare your data
data_to_test = [np.array(data[i]).ravel() for i in range(len(data))]

# Perform the Kruskal-Wallis H-test (one-way ANOVA on ranks)
H, p = kruskal(*data_to_test) #star used to unpack
print('H statistic:', H)
print('p value:', p) 

# p <0.05 so this means that at least one population mean for a biome is different from the population mean of another
# now we do a posthoc tst to find out which are different 

# Perform Dunn's test for pairwise comparisons
# 'p_adjust' parameter is for multiple test correction
posthoc = sp.posthoc_dunn(data_to_test, p_adjust='bonferroni')

print(posthoc) #all significantly different from each other. 

# I performed a Kruskal-Wallis non-parametric test on the monthly isoprene factor results from 2005 to 2014 for each biome. 
# The results indicated that there is a significant difference in the median isoprene factor across different biomes (H = 278324, p < 0.001).

# To further investigate which specific pairs of biomes exhibited significant differences in isoprene factors, 
# I conducted a post hoc analysis using Dunn's test with Bonferroni correction for multiple comparisons. 
# This analysis revealed that all pairs of biomes were significantly different from each other at a 5% confidence level.

# import seaborn as sns

# # Set the labels
# labels = biome_labels

# # Create the heatmap
# plt.figure(figsize=(8, 6))
# sns.heatmap(posthoc, xticklabels=labels, yticklabels=labels, annot=True, cmap='coolwarm', cbar=True)

# # Show the plot
# plt.title('Pairwise p-value matrix')
# plt.show()
