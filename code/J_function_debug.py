####################################################

from pyrealm import pmodel
import numpy as np
import matplotlib.pyplot as plt

# Make an np array of varying ppfd values
ppfd_array = np.arange(0, 300500, 500)
rep = len(ppfd_array) #lengths of each var needs to match so we need to get out length of ppfd_array

# Fill an np.array of the size of timesteps in rep with ta reasonable value
filled_mean_temp = np.full(rep, 4)
filled_mean_co2 = np.full(rep, 390)
filled_mean_patm = np.full(rep, 97367)
filled_mean_vpd = np.full(rep, 1119)
filled_mean_fapar = np.full(rep, 0.32)


# Calculate photosynthetic environment with constant environmental values
env = pmodel.PModelEnvironment(tc = filled_mean_temp, co2 = filled_mean_co2, patm = filled_mean_patm, vpd = filled_mean_vpd)
env.summarize()

# Run the P model
model = pmodel.PModel(env, method_jmaxlim = 'smith19') #run with smith19 to be consistent with J calculation
model.estimate_productivity(filled_mean_fapar, ppfd_array)
model.summarize()

# Get necessary variables for running the models
gammastar = env.gammastar
ci = model.optchi.ci
vcmax = model.vcmax
kmm = env.kmm
a_j = (model.vcmax / model.optchi.mjoc) * model.optchi.mj #no easy way to get this directly, backcalculating is the best way


# Define to J and Jv functions
def find_J(a_j, ci, gammastar):
    """Coming from Smith et al., 2019 (Eq. 4)
    Units for ci ang gammastar are Pa, which is consistent with the Pmodel
    Units for Aj is umol m-2 s-1, could not find units in pmodel, assumed same as in farquhar so consistent"""
    J = 4 * a_j * ((ci + 2 * gammastar) / (ci - gammastar))
    return J

def find_Jv(vcmax, ci, gammastar, kmm):
    """Comming from Morfopoulos et al, 2014 (New Physiologist)
    Units for ci, gammastar, and kmm is ubar, so convert to Pa for consistency
    Units for vcmax is umol m-2 s-1which is what is reported as the output by the pmodel code doccumentation (but not in the paper)"""
    Jv = 4 * vcmax * ((ci + 2 * gammastar) / (ci + kmm))
    return Jv


# Run functions to find J and Jv
J = find_J(a_j, ci, gammastar)
Jv = find_Jv(vcmax, ci*10, gammastar*10, kmm*10) #convert the relevant vars from Pa to ubar

# Calculate average J and Jv results
np.nanmean(J)
np.nanmean(Jv)
#theyre the same

# Create line graph for J and Jv
plt.plot(ppfd_array, J, color='blue', label='J')
plt.plot(ppfd_array, Jv, color='red', label='Jv')
plt.title('Line Graph')
plt.xlabel('ppfd')
plt.ylabel('Value')
plt.legend()
plt.show()
# Theyre the same