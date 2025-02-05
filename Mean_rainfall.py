import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import dask.array as da
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors

#load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc', )
ncum_g_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc', chunks={'time': 10})
ncum_g_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day3rf_jjas2020-24.nc', chunks={'time': 10})
ncum_g_5 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day5rf_jjas2020-24.nc', chunks={'time': 10})
#ncum_r_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2023.nc')
#ncum_r_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day3rf_jjas2023.nc')
#ncum_r_5 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day5rf_jjas2023.nc')

#change lat and lon to match
ncum_g_1 = ncum_g_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_3 = ncum_g_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_5 = ncum_g_5.rename({'latitude': 'lat', 'longitude': 'lon'})
#ncum_r_1 = ncum_r_1.rename({'latitude': 'lat', 'longitude': 'lon'})
#ncum_r_3 = ncum_r_3.rename({'latitude': 'lat', 'longitude': 'lon'})
#ncum_r_5 = ncum_r_5.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_g_1_regridded = ncum_g_1.interp_like(obs, method='nearest')
ncum_g_3_regridded = ncum_g_3.interp_like(obs, method='nearest')
ncum_g_5_regridded = ncum_g_5.interp_like(obs, method='nearest')
#ncum_r_regridded_1 = ncum_r_1.interp_like(obs, method='nearest')
#ncum_r_regridded_3 = ncum_r_3.interp_like(obs, method='nearest')
#ncum_r_5_regridded = ncum_r_5.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))#62
ncum_g_1_regridded = ncum_g_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_3_regridded = ncum_g_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_5_regridded = ncum_g_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

# mean rainfall
obs_mean = obs.mean(dim='time')
ncum_g_1_mean = ncum_g_1_regridded.mean(dim='time')
ncum_g_3_mean = ncum_g_3_regridded.mean(dim='time')
ncum_g_5_mean = ncum_g_5_regridded.mean(dim='time')
#ncum_r_1_mean = ncum_r_1_regridded.mean(dim='time')
#ncum_r_3_mean = ncum_r_3_regridded.mean(dim='time')
#ncum_r_5_mean = ncum_r_5_regridded.mean(dim='time')
'''
#plot
fig, axes = plt.subplots(2, 2, figsize=(15, 10), subplot_kw={'projection': '3d'})

obs_mean['rf'].plot(ax=axes[0, 0], cmap='viridis')
axes[0, 0].set_title('Observed Mean Rainfall')

ncum_g_1_mean['APCP_Surface'].plot(ax=axes[0, 1], cmap='viridis')
axes[0, 1].set_title('NCUM G Day 1 Mean Rainfall')

ncum_g_3_mean['APCP_Surface'].plot(ax=axes[1, 0], cmap='viridis')
axes[1, 0].set_title('NCUM G Day 3 Mean Rainfall')

ncum_g_5_mean['APCP_Surface'].plot(ax=axes[1, 1], cmap='viridis')
axes[1, 1].set_title('NCUM G Day 5 Mean Rainfall')

plt.tight_layout()
plt.show()
'''
plt.figure(figsize=(14, 12))

# Plot for Observed Mean Rainfall
ax1 = plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
obs_mean['rf'].plot(ax=ax1, cmap='Blues', transform=ccrs.PlateCarree())
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax1.set_title('Mean Rainfall Observation')

# Plot for NCUM-G Day 1 Mean Rainfall
ax2 = plt.subplot(2, 2, 2, projection=ccrs.PlateCarree())
ncum_g_1_mean['APCP_surface'].plot(ax=ax2, cmap='Blues', transform=ccrs.PlateCarree())
ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax2.set_title('Mean Rainfall NCUM-G Day 1')

# Plot for NCUM-G Day 3 Mean Rainfall
ax3 = plt.subplot(2, 2, 3, projection=ccrs.PlateCarree())
ncum_g_3_mean['APCP_surface'].plot(ax=ax3, cmap='Blues', transform=ccrs.PlateCarree())
ax3.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax3.set_title('Mean Rainfall NCUM-G Day 3')

# Plot for NCUM-G Day 5 Mean Rainfall
ax4 = plt.subplot(2, 2, 4, projection=ccrs.PlateCarree())
ncum_g_5_mean['APCP_surface'].plot(ax=ax4, cmap='Blues', transform=ccrs.PlateCarree())
ax4.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax4.set_title('Mean Rainfall NCUM-G Day 5')

plt.tight_layout()
plt.show()