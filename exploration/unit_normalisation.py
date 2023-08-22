import numpy as np
import netCDF4
import matplotlib.pyplot as plt

def gridcell_to_m2(latitude):
    """Calculate the number of square meters in each grid cell in a 0.5 by 0.5 degree grid."""
    # Latitude lengths stay at about 111.1 km per degree
    half_degree_lat = 111111.1 / 2

    # Equation to get the length of each half degree of longitude (changes with latitude)
    half_degree_lon = np.cos(np.deg2rad(latitude)) * (111111.1 / 2)

    meter_square_gridcell = half_degree_lat * half_degree_lon
    return meter_square_gridcell

# Read netCDF file and access data
HCHO = netCDF4.Dataset("../data/OMI_iso_estimate/mon_average_OMI.nc")
HCHO.set_auto_mask(False)

# Convert HCHO from isoprene/grid call to per m2
latitudes = [float(lat) for lat in HCHO['lat'][:]]  # Convert latitude values to float

# Calculate the number of square meters in each grid cell
no_m2_in_grid = [gridcell_to_m2(lat) for lat in latitudes]

# Repeat the number of square meters across each longitude
tiled = np.tile(no_m2_in_grid, (len(HCHO['lon'][:]), 1)).T

# Calculate the isoprene emitted per square meter
isoprene_per_m2 = HCHO["EMworldC5H8"][:] / tiled

# Filter out zeros and set them to NaN
isoprene_per_m2[isoprene_per_m2 == 0] = np.nan

# Filter out the outlying HCHO data
isoprene_per_m2[isoprene_per_m2 > 0.003] = np.nan

# Calculate statistics
isoprene_per_m2_max = np.nanmax(isoprene_per_m2) * 1e+3
isoprene_per_m2_mean = np.nanmean(isoprene_per_m2) * 1e+3
isoprene_per_m2_max/isoprene_per_m2_mean

# Filter out NaN values for plotting
filtered_data = isoprene_per_m2[~np.isnan(isoprene_per_m2)]

# Create histogram
hist, bins = np.histogram(filtered_data, bins='auto')

# Plot histogram
plt.hist(filtered_data, bins='auto')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram without NaNs')
plt.show()

# Import the Model data
model = netCDF4.Dataset("../results/iso_result_OMI_filterd_mon_mean.nc")
model.set_auto_mask(True)

# Get the isoprene estimate data
m_iso = model["J_times_IspS"][:]

# Set filled values and zeros to NaN
m_iso[m_iso > 8.3e+35] = np.nan
m_iso[m_iso == 0] = np.nan

# Calculate statistics
m_iso_max = np.nanmax(m_iso)
m_iso_mean = np.nanmean(m_iso)
m_iso_max/m_iso_mean

# Filter out NaN values for plotting
filtered_data = m_iso[~np.isnan(m_iso)]

# Create histogram
hist, bins = np.histogram(filtered_data, bins='auto')

# Plot histogram
plt.hist(filtered_data, bins='auto')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram without NaNs')
plt.show()


# Import the CRis data
CRIS = netCDF4.Dataset("../data/2012to2020_CrIS_Isoprene.nc")
CRIS.set_auto_mask(True)

# Get the isoprene estimate data
CRIS_iso = CRIS["Isop"][:]

CRIS_iso.close()

# Set zeros to NaN
CRIS_iso[CRIS_iso <= 0] = np.nan

# Calculate statistics
CRIS_max = np.nanmax(CRIS_iso)/1e+16
CRIS_mean = np.nanmean(CRIS_iso)/1e+16
CRIS_max/CRIS_mean

# Filter out NaN values for plotting
filtered_data = CRIS_iso[~np.isnan(CRIS_iso)]

# Create histogram
hist, bins = np.histogram(filtered_data, bins='auto')

# Plot the histogram
plt.hist(filtered_data, bins='auto')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram without NaNs')
plt.show()

def scale_data(data):
    min_val = 0
    max_val = np.nanmax(data)
    scaled_data = (data - min_val) / (max_val - min_val)
    return scaled_data

CRIS_scaled = scale_data(CRIS_iso)
np.nanmean(CRIS_scaled)

OMI_scaled = scale_data(isoprene_per_m2)
np.nanmean(OMI_scaled)

MODEL_scaled = scale_data(m_iso)
np.nanmean(MODEL_scaled)


import matplotlib.pyplot as plt
import cartopy.crs as ccrs

i =1
# Create the subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Plot the first image
image1 = axes[0].imshow(np.flipud(OMI_scaled[i,:,:]), origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
axes[0].coastlines()
axes[0].set_global()
axes[0].set_title(f'OMI {i}', fontweight='bold')

# Plot the second image
image2 = axes[1].imshow(CRIS_scaled[i,:,:], origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
axes[1].coastlines()
axes[1].set_global()
axes[1].set_title(f'CRIS {i}', fontweight='bold')


# Plot the third image
image3 = axes[2].imshow(MODEL_scaled[i,:,:], origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
axes[2].coastlines()
axes[2].set_global()
axes[2].set_title(f'Model {i}', fontweight='bold')

# Add colorbar
cbar = fig.colorbar(image3, ax=axes.ravel().tolist(), fraction=0.022, pad=0.03, location='bottom')

# Show the plot
plt.show()



import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Define the desired values of i
i_values = range(12)

for i in i_values:
    # Create the subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Plot the first image
    image1 = axes[0].imshow(np.flipud(OMI_scaled[i,:,:]), origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
    axes[0].coastlines()
    axes[0].set_global()
    axes[0].set_title(f'OMI {i}', fontweight='bold')

    # Plot the second image
    image2 = axes[1].imshow(CRIS_scaled[i,:,:], origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
    axes[1].coastlines()
    axes[1].set_global()
    axes[1].set_title(f'CRIS {i}', fontweight='bold')

    # Plot the third image
    image3 = axes[2].imshow(MODEL_scaled[i,:,:], origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
    axes[2].coastlines()
    axes[2].set_global()
    axes[2].set_title(f'Model {i}', fontweight='bold')

    # Add colorbar
    cbar = fig.colorbar(image3, ax=axes.ravel().tolist(), fraction=0.022, pad=0.03, location='bottom')

    # Save the figure
    filename = f"Figures/global_{i}.png"
    plt.savefig(filename)

    # Close the figure
    plt.close(fig)