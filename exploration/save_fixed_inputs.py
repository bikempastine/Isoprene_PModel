import numpy as np
import netCDF4
from scipy.stats import skew

## Clean and Save Model Data
# Open the original Model data
Model_cdf = netCDF4.Dataset("../data/model_2012_to_2014.nc")
Model_cdf.set_auto_mask(False)

# Extract the Model data
Model = Model_cdf['J_times_IspS'][:]


# Clean the Model data by setting invalid values to NaN
Model[(Model > 8.3e+35) | (Model == 0)] = np.nan


# Create a new netCDF file to write the cleaned Model data
new_model_cdf = netCDF4.Dataset("clean_data/Model_comp_clean.nc", 'w')

# Create dimensions in the new netCDF file
lat_dim = new_model_cdf.createDimension('lat', Model.shape[1])
lon_dim = new_model_cdf.createDimension('lon', Model.shape[2])
time_dim = new_model_cdf.createDimension('time', Model.shape[0])

# Create variables in the new netCDF file
lat_var = new_model_cdf.createVariable('lat', 'f4', ('lat',))
lon_var = new_model_cdf.createVariable('lon', 'f4', ('lon',))
time_var = new_model_cdf.createVariable('time', 'i4', ('time',))
model_var = new_model_cdf.createVariable('model_data', 'f4', ('time', 'lat', 'lon',))

# Assign values to variables in the new netCDF file
lat_var[:] = Model_cdf.variables['lat'][:]
lon_var[:] = Model_cdf.variables['lon'][:]
time_var[:] = Model_cdf.variables['time'][:]
model_var[:] = Model

# Close the new netCDF file and the original Model data file
new_model_cdf.close()
Model_cdf.close()


## Clean and Save CrIS Data
# Open the original CrIS data
CrIS_cdf = netCDF4.Dataset("../data/CrIS_2012_to_2014.nc")
CrIS_cdf.set_auto_mask(False)

# Extract the CrIS data
CRIS_flip = CrIS_cdf['Isop'][:]
CRIS_lon_flip = CrIS_cdf['lon'][:]


# Rearrange CrIS data to account for the longitude split
halfway = np.shape(CRIS_flip)[2] // 2
first_half = CRIS_flip[:, :, :halfway]
second_half = CRIS_flip[:, :, halfway:]
CrIS = np.concatenate((second_half, first_half), axis=2)
CrIS[CrIS <= 0] = np.nan


first_half_lon = CRIS_lon_flip[:halfway]
second_half_lon = CRIS_lon_flip[halfway:]
CRIS_lon = np.concatenate((second_half_lon, first_half_lon))


# Create a new netCDF file to write the cleaned CrIS data
new_cris_cdf = netCDF4.Dataset("clean_data/CrIS_comp_clean.nc", 'w')

# Create dimensions in the new netCDF file
lat_dim = new_cris_cdf.createDimension('lat', CrIS.shape[1])
lon_dim = new_cris_cdf.createDimension('lon', CrIS.shape[2])
time_dim = new_cris_cdf.createDimension('time', CrIS.shape[0])

# Create variables in the new netCDF file
lat_var = new_cris_cdf.createVariable('lat', 'f4', ('lat',))
lon_var = new_cris_cdf.createVariable('lon', 'f4', ('lon',))
time_var = new_cris_cdf.createVariable('time', 'i4', ('time',))
cris_iso = new_cris_cdf.createVariable('cris_iso', 'f4', ('time', 'lat', 'lon',))

# Assign values to variables in the new netCDF file
lat_var[:] = CrIS_cdf.variables['lat'][:]
lon_var[:] = CRIS_lon
time_var[:] = CrIS_cdf.variables['time'][:]
cris_iso[:] = CrIS

# Close the new netCDF file and the original CrIS data file
new_cris_cdf.close()
CrIS_cdf.close()


## Clean and Save OMI Data
# Open the original OMI data
OMI_cdf = netCDF4.Dataset("../data/OMI_iso_estimate/OMI_2012_to_2014.nc")
OMI_cdf.set_auto_mask(False)

# Extract the latitude values from OMI data
latitudes = [float(lat) for lat in OMI_cdf['lat'][:]]

def gridcell_to_m2(latitude):
    """
    Calculate the number of square meters in each grid cell in a 0.5 by 0.5 degree grid.
    
    Args:
        latitude (float): Latitude value.
        
    Returns:
        float: Number of square meters in each grid cell.
    """
    half_degree_lat = 111111.1 / 2
    half_degree_lon = np.cos(np.deg2rad(latitude)) * (111111.1 / 2)
    meter_square_gridcell = half_degree_lat * half_degree_lon
    return meter_square_gridcell

# Calculate the number of square meters in each grid cell
no_m2_in_grid = [gridcell_to_m2(lat) for lat in latitudes]
tiled = np.tile(no_m2_in_grid, (len(OMI_cdf['lon'][:]), 1)).T

# Calculate the isoprene emitted per square meter
OMI = OMI_cdf["EMworldC5H8"][:] / tiled

# Filter out the outlying HCHO data
OMI[(OMI > 0.003) | (OMI < 1.16729704e-07) | (OMI == 0)] = np.nan


# Create a new netCDF file to write the cleaned OMI data
new_omi_cdf = netCDF4.Dataset("clean_data/OMI_comp_clean.nc", 'w')

# Create dimensions in the new netCDF file
lat_dim = new_omi_cdf.createDimension('lat', OMI.shape[1])
lon_dim = new_omi_cdf.createDimension('lon', OMI.shape[2])
time_dim = new_omi_cdf.createDimension('time', OMI.shape[0])

# Create variables in the new netCDF file
lat_var = new_omi_cdf.createVariable('lat', 'f4', ('lat',))
lon_var = new_omi_cdf.createVariable('lon', 'f4', ('lon',))
time_var = new_omi_cdf.createVariable('time', 'i4', ('time',))
omi_iso = new_omi_cdf.createVariable('omi_iso', 'f4', ('time', 'lat', 'lon',))

# Assign values to variables in the new netCDF file
lat_var[:] = OMI_cdf.variables['lat'][:] *-1
lon_var[:] = OMI_cdf.variables['lon'][:]
time_var[:] = OMI_cdf.variables['time'][:]
omi_iso[:] = np.flip(OMI, axis=1)

# Close the new netCDF file and the original OMI data file
new_omi_cdf.close()
OMI_cdf.close()

new_omi_cdf = netCDF4.Dataset("clean_data/OMI_comp_clean.nc")
np.flip(new_omi_cdf.variables['lat'][:]) *-1