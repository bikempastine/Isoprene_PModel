"""
Example of plotting the HCHO-Isoprene data, including going from /grid cell to /m2
"""

# Load in the required packages
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors

## Load the HCHO data
HCHO = netCDF4.Dataset("../data/Top-down_C5H8_OMI_2010.nc")
HCHO.set_auto_mask(False)

## Conveert isoprene values from the units kg/gridcell/month to kg/m2/month
def gridcell_to_m2(latitude):
    "Find the number of m2 in each grid cell in a 0.5 by 0.5 degree grid"
    half_degree_lat = 111111.1/2 # latitude lengths stay at about 111.1km per degree
    half_degree_lon = np.cos(np.deg2rad(latitude)) * (111111.1/2) # equation to get the length of each half degree of longitude (changes with latitude)
    meter_square_gridcell = half_degree_lat * half_degree_lon
    return meter_square_gridcell

latitudes = [float(lat) for lat in HCHO['lat'][:]] #get each latitude degree as a float 
no_m2_in_grid = [gridcell_to_m2(lat) for lat in latitudes] #gets the number of m2 blocks in each grid cell
tiled = np.tile(no_m2_in_grid, (len(HCHO['lon'][:]), 1)).T #repeat across each longitude as distance remains the same across each latitude degree

# Get the isoprene emmited per m2
isoprene_per_m2 = (HCHO["EMworldC5H8"][0, :, :])/(tiled) #January 2010

# Flip the map: comes in upside down
flipped_isoprene = np.flipud(isoprene_per_m2) #January, 2010


## Apply a land mask


## Save the data
# Create a new netCDF file for storing the variables
output_file = netCDF4.Dataset("../results/HCHO_Jan_2010.nc", "w")

# Create dimensions in the output netCDF file
output_file.createDimension('lat', HCHO.dimensions["lat"].size)  
output_file.createDimension('lon', HCHO.dimensions["lon"].size)  
output_file.createDimension('time', None)  # unlimited axis (can be appended to).

# Copy dimention variables to the output file from ds
output_file.createVariable("time", HCHO["time"].dtype, ("time",))[:] = 1
output_file.createVariable("lat", HCHO["lat"].dtype, ("lat",))[:] = HCHO.variables['lat'][:]
output_file.createVariable("lon", HCHO["lon"].dtype, ("lon",))[:] = HCHO.variables['lon'][:]

# Create and assign data to variables in the output file
output_file.createVariable("HCHO_Iso", isoprene_per_m2.dtype, ("time", "lat", "lon"))[:] = isoprene_per_m2

# Close the netCDF to prevent corruption
output_file.close()


## Plotting
# Create a new figure and axis
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Define colormap and normalize the data
cmap = plt.cm.get_cmap('PuRd')
norm = colors.Normalize(vmin = flipped_isoprene.min(), vmax=np.percentile(flipped_isoprene[flipped_isoprene != 0], 75)) # want to scle the colourbar better

# Plot the data on the map
image = ax.imshow(flipped_isoprene , origin="lower", extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
ax.coastlines()

# Add a colorbar
cbar = plt.colorbar(image, ax=ax)
cbar.set_label('kg/ m2 / month')

# Set the map extent
ax.set_global()

# Add title and labels
plt.title('HCHO-Isoprene January 2010 ')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

plt.savefig('../results/HCHO_Jan_2010')
plt.show()

# Close the netCDF to avoid corruption
HCHO.close()
