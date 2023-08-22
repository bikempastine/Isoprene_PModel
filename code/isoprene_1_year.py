
import math
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
from pyrealm import pmodel
import isoprene_functions as Iso_fun


def run_isoprene_model(input_file, output_file, reference_temp):
    """
    Run the isoprene model and store the results in a netCDF file.

    Inputs:
    - input_file: Path to the input netCDF file.
    - output_file: Path to save the output netCDF file.
    - reference_temp: Reference temperature to run the IspS_function at.

    Outputs:
    - Saves the results of the isoprene model in the output netCDF file.
    """

    # Load the input dataset
    ds = netCDF4.Dataset(input_file)
    ds.set_auto_mask(False)

    # Extract the required variables
    temp = ds["temp"][:]
    co2 = ds["CO2"][:]
    elev = ds["elevation"][:]
    vpd = ds["VPD"][:]
    fapar = ds["fAPAR"][:]
    ppfd = ds["ppfd"][:]

    # Manipulate the inputs
    patm = pmodel.calc_patm(elev)
    temp[temp < -25] = np.nan
    vpd[vpd < 0] = 0

    # Run the PModel
    env = pmodel.PModelEnvironment(tc=temp, co2=co2, patm=patm, vpd=vpd)
    model = pmodel.PModel(env, method_jmaxlim='smith19')
    model.estimate_productivity(fapar, ppfd)

    # Get the necessary variables for J and Jv
    gammastar = env.gammastar
    ci = model.optchi.ci
    vcmax = model.vcmax
    kmm = env.kmm
    a_j = (model.vcmax / model.optchi.mjoc) * model.optchi.mj

    # Get J and Jv
    J = Iso_fun.find_J(a_j, ci, gammastar)
    Jv = Iso_fun.find_Jv(vcmax, ci, gammastar, kmm)

    # Run the IspS function
    IspS = np.vectorize(Iso_fun.IspS_function)(temp, reference_temp)

    # Produce objects that will be stored
    GPP = model.gpp
    Jmax = model.jmax
    LUE = model.lue
    Vcmax = model.vcmax

    J_times_IspS_umol = J * IspS

    # Create a new netCDF file for storing the variables
    output = netCDF4.Dataset(output_file, "w")

    # Create dimensions in the output netCDF file
    output.createDimension('lat', ds.dimensions["latitude"].size)
    output.createDimension('lon', ds.dimensions["longitude"].size)
    output.createDimension('time', None)  # unlimited axis (can be appended to).

    # Copy dimension variables to the output file from ds
    output.createVariable("time", ds["Time"].dtype, ("time",))[:] = ds.variables['Time'][:]
    output.createVariable("lat", ds["latitude"].dtype, ("lat",))[:] = ds.variables['latitude'][:]
    output.createVariable("lon", ds["longitude"].dtype, ("lon",))[:] = ds.variables['longitude'][:]

    # Create and assign data to variables in the output file
    output.createVariable("gpp", GPP.dtype, ("time", "lat", "lon"))[:] = GPP
    output.createVariable("jmax", Jmax.dtype, ("time", "lat", "lon"))[:] = Jmax
    output.createVariable("lue", LUE.dtype, ("time", "lat", "lon"))[:] = LUE
    output.createVariable("vcmax", Jv.dtype, ("time", "lat", "lon"))[:] = Vcmax
    output.createVariable("J_times_IspS", J_times_IspS_umol.dtype, ("time", "lat", "lon"))[:] = J_times_IspS_umol
    output.createVariable("J", J.dtype, ("time", "lat", "lon"))[:] = J
    output.createVariable("Jv", Jv.dtype, ("time", "lat", "lon"))[:] = Jv

    # Close the netCDF to prevent corruption
    ds.close()
    output.close()

    # Print to the console that the run is complete
    return(print("Done and saved in:", output_file))


# Example usage
input_file = "../data/pmodel_inputs.nc"
output_file = "../results/Iso_test_results.nc"
reference_temp = 30
run_isoprene_model(input_file, output_file, reference_temp)

ds = netCDF4.Dataset(output_file)
ds.variables
ds.close()


import netCDF4
ds = netCDF4.Dataset("../data/SWdown_2005/SWdown_WFDE5_CRU_200501_v2.0.nc")
ds.variables
ds.close()