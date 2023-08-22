
from matplotlib import pyplot as plt
import numpy as np
import netCDF4
import cartopy.crs as ccrs
from numpy import ndarray as NDArray


# import data
aet_nc = netCDF4.Dataset("clean_data/aet_2005_2014.nc")
aet_nc.set_auto_mask(True)
aet = aet_nc['aet'][:]
aet_nc.close()

pet_nc = netCDF4.Dataset("clean_data/pet_2005_2014.nc")
pet_nc.set_auto_mask(True)
pet = pet_nc['pet'][:]
pet_nc.close()

sm_lim_nc = netCDF4.Dataset("clean_data/sm_lim_2005_2014.nc")
sm_lim_nc.set_auto_mask(True)
sm_lim = sm_lim_nc['sm_lim'][:]
sm_lim_nc.close()


def calc_soilmstress_stocker(
    soilm: NDArray,
    meanalpha: NDArray = np.array(1.0),
    const: PModelConst = PModelConst(),
    ) -> NDArray:
    r"""Calculate Stocker's empirical soil moisture stress factor.

    This function calculates a penalty factor :math:`\beta(\theta)` for well-watered GPP
    estimates as an empirically derived stress factor :cite:p:`Stocker:2020dh`. The
    factor is calculated as a function of relative soil moisture (:math:`m_s`, fraction
    of field capacity) and average aridity, quantified by the local annual mean ratio of
    actual over potential evapotranspiration (:math:`\bar{\alpha}`).

    The value of :math:`\beta` is defined relative to two soil moisture thresholds
    (:math:`\theta_0, \theta^{*}`) as:

      .. math::
        :nowrap:

        \[
            \beta =
                \begin{cases}
                    q(m_s - \theta^{*})^2 + 1,  & \theta_0 < m_s <= \theta^{*} \\
                    1, &  \theta^{*} < m_s,
                \end{cases}
        \]

    where :math:`q` is an aridity sensitivity parameter setting the stress factor at
    :math:`\theta_0`:

    .. math:: q=(1 - (a + b \bar{\alpha}))/(\theta^{*} - \theta_{0})^2

    Default parameters of :math:`a=0` and :math:`b=0.7330` are as described in Table 1
    of :cite:t:`Stocker:2020dh` specifically for the 'FULL' use case, with
    ``method_jmaxlim="wang17"``, ``do_ftemp_kphio=TRUE``.

    Note that it is possible to use the empirical soil moisture stress factor effect on
    GPP to back calculate realistic Jmax and Vcmax values within the calculations of the
    P Model. This is applied, for example, in the `rpmodel` implementation.
    The :mod:`pyrealm.pmodel` module treats this factor purely as a penalty that can be
    applied after the estimation of GPP.

    Args:
        soilm: Relative soil moisture as a fraction of field capacity
            (unitless). Defaults to 1.0 (no soil moisture stress).
        meanalpha: Local annual mean ratio of actual over potential
            evapotranspiration, measure for average aridity. Defaults to 1.0.
        const: Instance of :class:`~pyrealm.constants.pmodel_const.PModelConst`.

    PModel Parameters:
        theta0: lower bound of soil moisture
            (:math:`\theta_0`, ``soilmstress_theta0``).
        thetastar: upper bound of soil moisture
            (:math:`\theta^{*}`, ``soilmstress_thetastar``).
        a: aridity parameter (:math:`a`, ``soilmstress_a``).
        b: aridity parameter (:math:`b`, ``soilmstress_b``).

    Returns:
        A numeric value or values for :math:`\beta`

    Examples:
        >>> # Relative reduction (%) in GPP due to soil moisture stress at
        >>> # relative soil water content ('soilm') of 0.2:
        >>> round((calc_soilmstress(0.2) - 1) * 100, 5)
        -11.86667
    """

    # TODO - move soilm params into standalone param class for this function -
    #        keep the PModelConst cleaner?

    # Check inputs, return shape not used
    #_ = check_input_shapes(soilm, meanalpha)

    # Calculate outstress
    y0 = const.soilmstress_a + const.soilmstress_b * meanalpha
    beta = (1.0 - y0) / (const.soilmstress_theta0 - const.soilmstress_thetastar) ** 2
    outstress = 1.0 - beta * (soilm - const.soilmstress_thetastar) ** 2

    # Filter wrt to thetastar
    outstress = np.where(soilm <= const.soilmstress_thetastar, outstress, 1.0)

    # Clip
    outstress = np.clip(outstress, 0.0, 1.0)

    return outstress



alpha = aet/pet

mean_alpha_yearly = np.full((1,360, 720), fill_value = np.nan)

for i in range(10):
    mean_alpha = np.mean(alpha[0+(i*12):12+(i*12),:,:], axis =0)
    mean_alpha_rep = np.repeat(mean_alpha[np.newaxis, ...], 12, axis=0)
    mean_alpha_yearly = np.concatenate((mean_alpha_yearly, mean_alpha_rep), axis=0)

mean_alpha_yearly = mean_alpha_yearly[1:,:,:] 


sm_stress = calc_soilmstress_stocker(soilm= sm_lim, meanalpha=mean_alpha_yearly)

sm_stress_masked = np.ma.masked_array(sm_stress, mask=sm_lim.mask)
sm_stress_flip = np.flip(sm_stress_masked, axis=1)
sm_stress_flip_nan = sm_stress_flip.filled(np.nan)

np.save('clean_data/sm_stress_2005_2014.npy', sm_stress_flip_nan)



## reload the data
#loaded_array = np.load('clean_data/sm_stress_2005_2014.npy')

# plot to check
# fig = plt.figure(figsize=(10,5))
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
# img = ax.imshow(loaded_array[18,:,:], origin='lower', cmap='RdBu_r', extent=[-180, 180, -90, 90])
# ax.coastlines()
# ax.gridlines()
# plt.colorbar(img, ax=ax, orientation='horizontal', label='sm')
# #plt.savefig('Figures/cris_omi_correlation.png')
# plt.show()

np.save('my_array.npy', my_array)