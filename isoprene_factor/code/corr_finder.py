import numpy as np
import netCDF4
import matplotlib.pyplot as plt

# Open the Model data
Model_cdf = netCDF4.Dataset("../data/model_clean.nc")
Model_cdf.set_auto_mask(False)
model = Model_cdf['model_data'][:]
Model_cdf.close()

# Open the OMI data
OMI_cdf = netCDF4.Dataset("../data/OMI_clean.nc")
OMI_cdf.set_auto_mask(False)
omi = OMI_cdf['omi_iso'][:]
OMI_cdf.close()


def find_corr(data1, data2):
    x_valid = []
    y_valid = []
    for month in range(120):
        valid_mask = np.logical_and(~np.isnan(data1[month,:,:]), ~np.isnan(data2[month,:,:]))
        x_valid = np.concatenate([x_valid, data1[month,:,:][valid_mask]])
        y_valid =  np.concatenate([y_valid, data2[month,:,:][valid_mask]])
    
    corr_matrix = np.corrcoef(x_valid, y_valid)
    return corr_matrix[0,1]


find_corr(model,omi)


#################################
####### By SM_stress ############
#################################

soil_moisture_data = np.load('../data/sm_stress_2005_2014.npy')


def elementwise_bool_check(array1, array2, array3):
    result = np.empty((array1.shape[0], array1.shape[1]),dtype ='bool')
    for lat in range(array1.shape[0]):
        for lon in range(array1.shape[1]):
            result[lat,lon] = all([array1[lat,lon], array2[lat,lon], array3[lat,lon]])
    return result

def find_corr_sm(data1, data2, sm_stress):
    x_valid = []
    y_valid = []
    for month in range(120):
        
        if sm_stress == 'y' :
            sm_mask = (soil_moisture_data[month,:,:] < 1)
        elif sm_stress == 'n' :
            sm_mask = (soil_moisture_data[month,:,:] == 1)

        valid_mask = elementwise_bool_check(sm_mask, ~np.isnan(data1[month,:,:]), ~np.isnan(data2[month,:,:]))
        x_valid = np.concatenate([x_valid, data1[month,:,:][valid_mask]])
        y_valid =  np.concatenate([y_valid, data2[month,:,:][valid_mask]])
    
    corr_matrix = np.corrcoef(x_valid, y_valid)
    return corr_matrix[0,1]


find_corr_sm(omi, model, sm_stress = 'n')
find_corr_sm(omi, model, sm_stress = 'y')
