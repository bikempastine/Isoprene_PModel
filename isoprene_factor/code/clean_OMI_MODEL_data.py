import numpy as np
import netCDF4

def clean_save_model_data(input_file, output_file):
    # Open the original Model data
    Model_cdf = netCDF4.Dataset(input_file)
    Model_cdf.set_auto_mask(False)

    # Extract the Model data
    Model = Model_cdf['J_times_IspS'][:]

    # Clean the Model data by setting invalid values to NaN
    Model[(Model > 8.3e+35) | (Model == 0)] = np.nan

    # Create a new netCDF file to write the cleaned Model data
    new_model_cdf = netCDF4.Dataset(output_file, 'w')

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


def clean_save_omi_data(input_file, output_file):
    # Open the original OMI data
    OMI_cdf = netCDF4.Dataset(input_file)
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
    new_omi_cdf = netCDF4.Dataset(output_file, 'w')

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


clean_save_model_data('../data/model.nc', '../data/model_clean.nc')
clean_save_model_data('../data/model_monthly_average.nc', '../data/model_monthly_average_clean.nc')
clean_save_omi_data('../data/OMI.nc', '../data/OMI_clean.nc')
clean_save_omi_data('../data/OMI_monthly_average.nc', '../data/OMI_monthly_average_clean.nc')