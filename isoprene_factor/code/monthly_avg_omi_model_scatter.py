import numpy as np
import netCDF4
import cartopy.crs as ccrs
from matplotlib import pyplot as plt

# Open the Model data
Model_cdf = netCDF4.Dataset("../data/model_monthly_average_clean.nc")
Model_cdf.set_auto_mask(False)
model = Model_cdf['model_data'][:]
Model_cdf.close()

# Open the OMI data
OMI_cdf = netCDF4.Dataset("../data/OMI_monthly_average_clean.nc")
OMI_cdf.set_auto_mask(False)
omi = OMI_cdf['omi_iso'][:]
OMI_cdf.close()

# filter out when Model below 1 (I think here the error is driving the high F values)
# model[model < 1] = np.nan
# omi[model < 1] = np.nan

# model[omi < 3e-05] = np.nan
# omi[omi < 3e-05] = np.nan

# scatter plot
plt.scatter(model, omi, color='black', alpha=0.02, s=1)
plt.xlabel('model')
plt.ylabel('omi')
plt.title('Scatter Plot of omi vs model')

avg_F = np.nanmean(omi)/np.nanmean(model)
x_line = np.linspace(np.nanmin(model), np.nanmax(model), 100)

plt.plot(x_line, avg_F * x_line, color='red', linewidth=2, label = 'Isoprene F')
#plt.plot(x_line, (avg_F/2) * x_line, color='blue', linewidth=2, label = 'Isoprene F/2')
plt.plot(x_line, (avg_F/1.5) * x_line, color='green', linewidth=2, label = 'Isoprene F/1.5')
plt.legend(loc = 'right')
plt.show()


def flatten_and_rm_nans(data_list):
    A_filtered = np.array([])
    B_filtered = np.array([])
        
    for month in range(12):
        data_A = np.copy(data_list[0][month,:,:])
        data_B = np.copy(data_list[1][month,:,:])

        mask = np.logical_and(~np.isnan(data_A), ~np.isnan(data_B))

        # Concatenate the valid elements from each month
        A_filtered = np.concatenate((A_filtered, data_A[mask]))
        B_filtered = np.concatenate((B_filtered, data_B[mask]))

    return [A_filtered, B_filtered]


model_flat, omi_flat = flatten_and_rm_nans([model, omi])

coefficients = np.polyfit(model_flat, omi_flat, 1)
slope = coefficients[0]
intercept = coefficients[1]

# Generate points for the regression line
regression_line = slope * model_flat + intercept


plt.scatter(model_flat, omi_flat, color='black', alpha=0.02, s=1)
plt.xlabel('model')
plt.ylabel('omi')
plt.title('Scatter Plot of omi vs model')
avg_F = np.nanmean(omi)/np.nanmean(model)
x_line = np.linspace(np.nanmin(model), np.nanmax(model), 100)
y_line = avg_F * x_line
plt.plot(x_line, y_line, color='red', linewidth=2, label = 'F')
plt.plot(model_flat, regression_line, color='green', label='Regression Line')
# plt.text(x_line[80]-3, y_line[80]- 0.0003, 'Isoprene Factor', color='red', fontsize=10, ha='left', va='bottom', backgroundcolor='white')
plt.legend()
plt.show()




import matplotlib.pyplot as plt
import numpy as np

# Assuming 'model' and 'omi' are your two variables you want to plot, and 'biome' is your color coding variable
# Also assuming that model and omi are of the shape (12, 360, 720) and biome is of shape (360, 720)

# Load the biome map
biome_cdf = netCDF4.Dataset("../data/clean_Type3_biome.nc")
biome_cdf.set_auto_mask(False)
biome = biome_cdf['type 3'][:].astype(int)
biome_cdf.close()


# Flatten the arrays
model_flat = model[0,:,:].reshape(-1)
omi_flat = omi[0,:,:].reshape(-1)
# biome_flat = np.repeat(biome.reshape(-1), 12)  # Repeat biome values for each month
biome_flat = biome.reshape(-1)

# Create a mask for biome ids from 1 through 8
mask = np.isin(biome_flat, range(1, 9))

# Apply the mask
model_flat = model_flat[mask]
omi_flat = omi_flat[mask]
biome_flat = biome_flat[mask]

# Create a discrete colormap (replace with any colormap suitable for your data)
cmap = plt.get_cmap('tab10', 8)  # 'tab10' is a colormap with 10 distinct colors. We want 8 colors for 8 biomes.

# Create scatter plot
plt.figure(figsize=(10,8))
scatter = plt.scatter(model_flat, omi_flat, c=biome_flat, cmap=cmap, s=1)

# Create colorbar
cbar = plt.colorbar(scatter, ticks=range(1, 9))  # Adjust the ticks according to your biome range
cbar.set_label('Biome')

plt.show()



