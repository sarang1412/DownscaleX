import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#load data
ncum_r = xr.open_dataset('e:\\Dissertation\\ncumr_day1rf_jjas2023.nc')
ncum_g = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc')
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc' )

#change lat and lon to match
ncum_g = ncum_g.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_g_regridded = ncum_g.interp_like(obs, method='nearest')
ncum_r_regridded = ncum_r.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(62, 106))
ncum_g_regridded = ncum_g_regridded.sel(lat=slice(6, 41), lon=slice(62, 106))
ncum_r_regridded = ncum_r_regridded.sel(lat=slice(6, 41), lon=slice(62, 106))

#proper date format
obs['time'] = pd.to_datetime(obs['time'].values)
ncum_g_regridded['time'] = pd.to_datetime(ncum_g_regridded['time'].values)
ncum_r_regridded['time'] = pd.to_datetime(ncum_r_regridded['time'].values)

mean_obs = obs['rf'].mean(dim=('lat', 'lon')) 
mean_ncum_g = ncum_g_regridded['APCP_surface'].mean(dim=('lat', 'lon'))
mean_ncum_r = ncum_r_regridded['APCP_24'].mean(dim=('lat', 'lon'))

# biases
bias_ncum_g = mean_ncum_g - mean_obs
bias_ncum_r = mean_ncum_r - mean_obs

#yearly average
bias_ncum_g = bias_ncum_g.groupby('time.year').mean()
bias_ncum_r = bias_ncum_r.groupby('time.year').mean()

#plot
plt.figure(figsize=(12,6))
plt.plot(bias_ncum_g['year'],bias_ncum_g, label='NCUM-G bias')
plt.plot(bias_ncum_r['year'],bias_ncum_r, label='NCUM-R bias')
plt.title('Comparison of Bias')
plt.xlabel('Time')
plt.ylabel('Mean Rainfall Bias(mm)')
plt.legend()
plt.grid(True)
plt.show()