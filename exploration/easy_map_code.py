import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4

plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
plt.imshow(x, extent=[-180, 180, -90, 90],
           origin='lower', transform=ccrs.PlateCarree(), cmap='viridis')
ax.coastlines()
plt.colorbar(label='Data Values')
plt.show()


biome_cdf = netCDF4.Dataset("../data/clean_Type3_biome.nc")
biome_cdf.set_auto_mask(False)
biome = biome_cdf['type 3'][:]
biome_cdf.close()
x= biome
x[x != 7] = np.nan