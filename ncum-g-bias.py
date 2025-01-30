import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates

#load data
ncum_g = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc')
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc' )

#change lat and lon to match
ncum_g = ncum_g.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_g_regridded = ncum_g.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(62, 106))
ncum_g_regridded = ncum_g_regridded.sel(lat=slice(6, 41), lon=slice(62, 106))

# mean rainfall
mean_obs =obs['rf'].mean(dim='time') 
mean_ncum_g = ncum_g_regridded['APCP_surface'].mean(dim='time')

# biases
bias_ncum_g = mean_ncum_g - mean_obs

# Plot Bias NCUM-G
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.figure(figsize=(14, 6))
ax1 = plt.subplot(1, 2, 1, projection=ccrs.PlateCarree())
bias_ncum_g.plot(ax=ax1, cmap='RdBu',
                 vmin=-10, vmax=10, 
                 cbar_kwargs={'label': 'Bias (mm)'}
                 )
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax1.set_title('Mean Rainfall Bias NCUM-G')
plt.show()
