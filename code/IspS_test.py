# import the required packages
import numpy as np
import netCDF4
from matplotlib import pyplot as plt
import IspS_functions as IspS #import the functions from the IspS_functions.py script

REFERANCE_TEMP = 30 

## Test: SS in the range given in Chapter 2 of Catherine's PhD Figure 2.5

# Generate a range of numbers from 15 - 50 degrees C
x_values = range(15, 51)

# Apply the SS function to each value in the range and normalize
SS_values = np.array([IspS.SS(x) for x in x_values])
SS_normalized = (SS_values - np.min(SS_values)) / (np.max(SS_values) - np.min(SS_values)) #max/min normalization

# Create and view a plot, see that it is the same as the blue SS line
plt.plot(x_values, SS_normalized)
plt.xlabel('Temperature (C)')
plt.ylabel('Normalised Response')
plt.title('SS varying Tk')
plt.show()


# Apply the IspS function (SS(T)/SS(REFERANCE_TEMP)) to each value in the range 
IspS_values = np.array([IspS.IspS_function(x, REFERANCE_TEMP) for x in x_values])

# Create and view a plot, see that it is the same shape as the blue SS line
# This makes sense since we just devided by a constant value (SS 30), or whatever you chose
plt.plot(x_values, IspS_values)
plt.xlabel('Temperature in Celsius')
plt.ylabel('f(T)')
plt.title('f(T) evaluated at the reference temp')
plt.show()

## Test: test on satalite temperature data

# Load an example dataset containing the main variables (downloaded from David Orme's repository)
ds = netCDF4.Dataset("../data/pmodel_inputs.nc")
ds.set_auto_mask(False)
temp = ds["temp"][:]
ds.close()

# Apply the function to each value in the range
IspS_temp_results = np.vectorize(IspS.IspS_function)(temp, REFERANCE_TEMP)

# Display the results
ax = plt.imshow(IspS_temp_results[0, :, :], origin="lower", extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()
