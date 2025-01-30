import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates

#load data
ncum_r = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc')
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc' )

#change lat and lon to match
#ncum_r = ncum_r.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_r_regridded = ncum_r.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(62, 106))
ncum_r_regridded = ncum_r_regridded.sel(lat=slice(6, 41), lon=slice(62, 106))

# mean rainfall
mean_obs =obs['rf'].mean(dim='time') 
mean_ncum_r = ncum_r_regridded['APCP_20_24'].mean(dim='time')

# biases
bias_ncum_r = mean_ncum_r - mean_obs

# Plot Bias NCUM-G
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.figure(figsize=(14, 6))
ax1 = plt.subplot(1, 2, 1, projection=ccrs.PlateCarree())
bias_ncum_r.plot(ax=ax1, cmap='RdBu',
                 vmin=-10, vmax=10, 
                 cbar_kwargs={'label': 'Bias (mm)'}
                 )
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax1.set_title('Mean Rainfall Bias NCUM-R')
plt.show()
