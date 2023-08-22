import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd
from scipy.stats import pearsonr


# Read the CSV file into a DataFrame
data = pd.read_csv('../data/cleaned_regress_for_R.csv')

data = pd.read_csv('../data/lasso_results_2.csv')

# Display the first few rows of the DataFrame
print(data.head())


df = data.drop(columns=['time', 'lat', 'lon','F','sm_categorical','biome'])
df.head()
# Rename columns
df = df.rename(columns={'lue': 'LUE', 'jmax': 'Jmax', 'gpp': 'GPP', 'vcmax': 'Vcmax', 
                        'sm_continuous' :'SM stress', 'temp':'Temp', 'swdown' : 'SW down', 'fpar': 'FPAR'})

# Calculate correlation matrix
corr = df.corr()
p_value_matrix = df.apply(lambda x: df.apply(lambda y: pearsonr(x, y)[1]))

# Create a heatmap
plt.figure(figsize=(10, 8))  # Set figure size
sns.heatmap(corr, annot=True, cmap='coolwarm', vmin=-1, vmax=1)

# Display the heatmap
plt.show()


unique_vars = np.unique(data['Var1'])
unique_lambda = np.unique(data['lambda'])


# Setting matplotlib parameters
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 15


plt.figure(figsize=(8, 6))

# Specifying hex codes for colors; add or replace as needed
color_map = {
    'train_data.J': '#D32F2F',  
    'Biome': '#CBD6E2',
    'train_data.fpar': '#EBB8DD',  
    'train_data.gpp': '#388E3C',  # rich green
    'train_data.jmax': '#FBC02D',  # bright yellow
    'train_data.lue': '#53BDA5',  # vibrant purple
    'train_data.sm_continuous': '#E64A19',  # deep orange
    'train_data.swdown': '#7B1FA2',  
    'train_data.temp': '#0288D1',  
    'train_data.vcmax':   '#253342'
}

# Rename variables for the legend
legend_names = {
    'Biome': "Biomes",
    'train_data.J': "J",
    'train_data.fpar': 'FPAR', 
    'train_data.gpp': 'GPP', 
    'train_data.jmax': 'Jmax',
    'train_data.lue': 'LUE', 
    'train_data.sm_continuous': 'Soil Moisture Stress', 
    'train_data.swdown': 'Shortwave Radiation',
    'train_data.temp': 'Ambient Temperature', 
    'train_data.vcmax': 'Vcmax'
}

biome_plotted = False  # Track if we've already plotted a biome line for the legend

for idx, var in enumerate(unique_vars):
    mask = [v == var for v in data['Var1']]
    
    if 'train_data.biome' in var:
        # If we haven't added a biome line to the legend, add it and set the flag to True
        if not biome_plotted:
            plt.plot(np.log(data['lambda'][mask]), data['value'][mask], label=legend_names['Biome'], color=color_map['Biome'], linewidth=2)
            biome_plotted = True
        # Plot the biome lines, but don't add to the legend
        plt.plot(np.log(data['lambda'][mask]), data['value'][mask], color=color_map['Biome'], linewidth=2)
    else:
        plt.plot(np.log(data['lambda'][mask]), data['value'][mask], label=legend_names[var], color=color_map[var], linewidth=2)

plt.xlabel(r'$\log(\lambda)$')
plt.ylabel('Normalised Coefficient Value')
plt.ylim(-1, 1)
plt.title('Lasso Regression Coeficient Trace Plot')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=5)
#plt.tight_layout()
plt.savefig('trace_plot.png', dpi=300, bbox_inches='tight')
plt.show()
