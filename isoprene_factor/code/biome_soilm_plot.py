# need to first run 
# sm_stress_and_F.py
# biome_and_F.py

import matplotlib.pyplot as plt

# Create the box plot for the Biome data
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

    # Initialise the figure with two subplots side by side
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7, 8))

    # First subplot for Biome data
    ax1 = axs[0]
    ax1.boxplot(data_to_plot, vert=False, showfliers=False, showmeans=True,
                meanprops=meanpointprops,
                medianprops=medianprops, boxprops=boxprops, 
                whiskerprops=whiskerprops, capprops=capprops)  

    ax1.yaxis.tick_left()
    ax1.set_yticks(range(1, 9))
    ax1.set_yticklabels(biome_labels_to_plot)
    ax1.set_title(f'Isoprene Factor by Biome', fontweight='bold', loc = 'left')
    #ax1.set_xlabel(r'Isoprene Factor ×$10^{-6}$')

    # Second subplot for Soil Moisture Stress data
    ax2 = axs[1]
    bx = ax2.boxplot(isoprene_factor_parted, vert=False, showfliers=False, showmeans=True, patch_artist=True,
                meanprops=meanpointprops,
                medianprops=medianprops, boxprops=boxprops, 
                whiskerprops=whiskerprops, capprops=capprops)  

    for box in bx['boxes'][0:2]:
        box.set(facecolor='lightgrey')

    for box in bx['boxes'][2:7]:
        box.set(facecolor='white')

    ax2.set_yticks(range(1, 7))
    ax2.set_yticklabels([r'$\beta$ = 1', r'$\beta$ $<$ 1', r'1 $>$ $\beta$ $\geq$ 0.8', 
                         r'0.8 $>$ $\beta$ $\geq$ 0.6', r'0.6 $>$ $\beta$ $\geq$ 0.4', r'$\beta$ $<$ 0.4'])
    ax2.set_title(r'Isoprene Factor by Soil Moisture Stress ($\beta$)', fontweight='bold', loc= 'left')
    ax2.set_xlabel(r'Isoprene Factor ×$10^{-6}$', loc = 'right')

    # Adjust the subplots layout
    plt.subplots_adjust(wspace=0.5,left=0.25)

    # Set x-axis range to 0 to 150 for both plots
    for ax in axs:
        ax.set_xlim(0, 160)

    # Save the figure
    plt.savefig(filename, dpi=300, bbox_inches='tight')

    plt.show()


# Call the function to create the plots side by side
plot_biome_bound_f('../../final_figures/box_plots.png')

