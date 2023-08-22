import netCDF4
import numpy as np

def clean_and_resave_netcdf(original_data_path, new_file_path, var_name):
    """Clean data and save to a new NetCDF file"""
    
    # Open the original data
    original_cdf = netCDF4.Dataset(original_data_path)
    original_cdf.set_auto_mask(False)

    # Extract the data
    data_flip = original_cdf[var_name][:]
    lon_flip = original_cdf['lon'][:]

    # Rearrange data to account for the longitude split
    halfway = np.shape(data_flip)[1] // 2
    first_half = data_flip[ :, :halfway]
    second_half = data_flip[ :, halfway:]
    rearranged_data = np.concatenate((second_half, first_half), axis=1)

    first_half_lon = lon_flip[:halfway]
    second_half_lon = lon_flip[halfway:]
    rearranged_lon = np.concatenate((second_half_lon, first_half_lon))

    # Create a new netCDF file to write the cleaned data
    new_cdf = netCDF4.Dataset(new_file_path, 'w')

    # Create dimensions in the new netCDF file
    lat_dim = new_cdf.createDimension('lat', rearranged_data.shape[0])
    lon_dim = new_cdf.createDimension('lon', rearranged_data.shape[1])

    # Create variables in the new netCDF file
    lat_var = new_cdf.createVariable('lat', 'f4', ('lat',))
    lon_var = new_cdf.createVariable('lon', 'f4', ('lon',))
    new_var = new_cdf.createVariable(var_name, 'f4', ( 'lat', 'lon',))

    # Assign values to variables in the new netCDF file
    lat_var[:] = original_cdf.variables['lat'][:]
    lon_var[:] = rearranged_lon
    new_var[:] = rearranged_data

    # Close the new netCDF file and the original data file
    new_cdf.close()
    original_cdf.close()


clean_and_resave_netcdf(original_data_path = "../isoprene_factor/data/Type3_LC0.5.nc", 
                        new_file_path = "../isoprene_factor/data/clean_Type3_biome.nc", 
                        var_name = 'type 3')
