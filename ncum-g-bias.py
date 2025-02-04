import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors

#load data
ncum_g_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc')
ncum_g_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day3rf_jjas2020-24.nc')
ncum_g_5 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day5rf_jjas2020-24.nc')
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc' )

#change lat and lon to match
ncum_g_1 = ncum_g_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_3 = ncum_g_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_5 = ncum_g_5.rename({'latitude': 'lat', 'longitude': 'lon'})
#regrid
ncum_g_1_regridded = ncum_g_1.interp_like(obs, method='nearest')
ncum_g_3_regridded = ncum_g_3.interp_like(obs, method='nearest')
ncum_g_5_regridded = ncum_g_5.interp_like(obs, method='nearest')
#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_1_regridded = ncum_g_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_3_regridded = ncum_g_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_5_regridded = ncum_g_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

# mean rainfall
mean_obs =obs['rf'].mean(dim='time') 
mean_ncum_g_1 = ncum_g_1_regridded['APCP_surface'].mean(dim='time')
mean_ncum_g_3 = ncum_g_3_regridded['APCP_surface'].mean(dim='time')
mean_ncum_g_5 = ncum_g_5_regridded['APCP_surface'].mean(dim='time')

# biases
bias_ncum_g_1 = mean_ncum_g_1 - mean_obs
bias_ncum_g_3 = mean_ncum_g_3 - mean_obs
bias_ncum_g_5 = mean_ncum_g_5 - mean_obs

# Plot Bias NCUM-G
levels = np.arange(-10, 12, 1)
cmap = plt.get_cmap('RdBu', len(levels) - 1) 
norm = mcolors.BoundaryNorm(levels, cmap.N) 


plt.figure(figsize=(14, 6))
ax1 = plt.subplot(1, 3, 1, projection=ccrs.PlateCarree())
bias_ncum_g_1.plot(ax=ax1, cmap=cmap, norm=norm)
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax1.set_title('Mean Rainfall Bias NCUM-G_day_1')

ax2 = plt.subplot(1, 3, 2, projection=ccrs.PlateCarree())
bias_ncum_g_3.plot(ax=ax2, cmap=cmap, norm=norm)
ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax2.set_title('Mean Rainfall Bias NCUM-G_day_3')

ax3 = plt.subplot(1, 3, 3, projection=ccrs.PlateCarree())
bias_ncum_g_5.plot(ax=ax3, cmap=cmap, norm=norm, cbar_kwargs={'label': 'Bias (mm)'})
ax3.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax3.set_title('Mean Rainfall Bias NCUM-G_day_5')

plt.tight_layout()
plt.tight_layout()
plt.show()
