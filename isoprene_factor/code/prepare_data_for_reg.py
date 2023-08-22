import numpy as np
import netCDF4

model_cdf = netCDF4.Dataset("../../results/iso_result_OMI_filtered.nc")
model_cdf.set_auto_mask(False)
biome = biome_cdf['type 3'][:]
biome_cdf.close()


soil_moisture_data = np.load('../data/sm_stress_2005_2014.npy')
isoprene_factor = np.load('../data/isoprene_factor_2005_2014.npy')

def elementwise_bool_check(array1, array2, array3):
    result = np.empty((array1.shape[0], array1.shape[1]),dtype ='bool')
    for lat in range(array1.shape[0]):
        for lon in range(array1.shape[1]):
            result[lat,lon] = all([array1[lat,lon], array2[lat,lon], array3[lat,lon]])
    return result


def filter_data(data1, sm_stress):
    x_valid = []
    for month in range(120):
        
        if sm_stress == 'y' :
            sm_mask = (soil_moisture_data[month,:,:] < 1)
        elif sm_stress == 'n' :
            sm_mask = (soil_moisture_data[month,:,:] == 1)

        valid_mask = elementwise_bool_check(sm_mask, ~np.isnan(data1[month,:,:]), ~np.isnan(soil_moisture_data[month,:,:]))
        x_valid = np.concatenate([x_valid, data1[month,:,:][valid_mask]])
    
    return x_valid


n_sm_F = filter_data(isoprene_factor* 1e6, sm_stress = 'n').tolist()
y_sm_F = filter_data(isoprene_factor* 1e6, sm_stress = 'y').tolist()
F_data = n_sm_F + y_sm_F

n_rep = ['n'] * len(n_sm_F)
y_rep = ['y'] * len(y_sm_F)
catagory_sm = n_rep + y_rep


def filter_data(data1, soil_moisture_data):
    x_valid = []
    sm_valid = []
    for month in range(120):
        sm_mask = (soil_moisture_data[month,:,:] < 1)

        valid_mask = elementwise_bool_check(sm_mask, ~np.isnan(data1[month,:,:]), ~np.isnan(soil_moisture_data[month,:,:]))
        x_valid = np.concatenate([x_valid, data1[month,:,:][valid_mask]])
        sm_valid = np.concatenate([sm_valid, soil_moisture_data[month,:,:][valid_mask]])
    
    return x_valid.tolist(), sm_valid.tolist()

F_w_sm_stress , sm_valid = filter_data(isoprene_factor* 1e6,soil_moisture_data )
len(F_w_sm_stress)
len(sm_valid)


# by biome
def F_by_biome(data1, biome):
    # Create an array to store Pearson correlation coefficients for each biome
    biome_list = np.array([])
    F_by_biome_list = np.array([])

    for biome_id in range(1, 9):  # Iterate through each biome (biome IDs from 1 to 8)
        data1_biome_filtered = np.array([])

        for month in range(120):  # Iterate through each month (months from 0 to 11)
            data1_month_copy = np.copy(data1[month, :, :])

            # Apply biome filter
            biome_mask = (biome == biome_id)

            # Create a mask to keep only valid (non-NaN) elements in both data arrays
            mask = elementwise_bool_check(biome_mask, ~np.isnan(data1_month_copy), ~np.isnan(biome))

            # Concatenate the valid elements from each month
            data1_biome_filtered = np.concatenate((data1_biome_filtered, data1_month_copy[mask]))


        F_by_biome_list = np.concatenate((F_by_biome_list, data1_biome_filtered))
        biome_list = np.concatenate((biome_list, np.full(data1_biome_filtered.size, biome_id)))

    return F_by_biome_list.tolist(), biome_list.tolist()

# Load the biome map
biome_cdf = netCDF4.Dataset("../data/clean_Type3_biome.nc")
biome_cdf.set_auto_mask(False)
biome = biome_cdf['type 3'][:]
biome_cdf.close()



F_by_biome , biome_list = F_by_biome(isoprene_factor*1e6, biome)

np.any(np.isnan(F_by_biome))