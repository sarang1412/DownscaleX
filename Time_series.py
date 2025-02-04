import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc' )
ncum_g_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc')
ncum_g_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day3rf_jjas2020-24.nc')
ncum_g_5 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day5rf_jjas2020-24.nc')
#ncum_r_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2023.nc')
#ncum_r_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2023.nc')
#ncum_r_5 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2023.nc')

#change lat and lon to match
ncum_g_1 = ncum_g_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_3 = ncum_g_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_5 = ncum_g_5.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_g_1_regridded = ncum_g_1.interp_like(obs, method='nearest')
ncum_g_3_regridded = ncum_g_3.interp_like(obs, method='nearest')
ncum_g_5_regridded = ncum_g_5.interp_like(obs, method='nearest')
#ncum_r_1_regridded = ncum_r_1.interp_like(obs, method='nearest')
#ncum_r_3_regridded = ncum_r_3.interp_like(obs, method='nearest')
#ncum_r_5_regridded = ncum_r_5.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_1_regridded = ncum_g_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_3_regridded = ncum_g_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_5_regridded = ncum_g_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_3_regridded = ncum_r_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_5_regridded = ncum_r_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

#proper date format
obs['time'] = pd.to_datetime(obs['time'].values)
ncum_g_1_regridded['time'] = pd.to_datetime(ncum_g_1_regridded['time'].values)
ncum_g_3_regridded['time'] = pd.to_datetime(ncum_g_3_regridded['time'].values)
ncum_g_5_regridded['time'] = pd.to_datetime(ncum_g_5_regridded['time'].values)
#ncum_r_1_regridded['time'] = pd.to_datetime(ncum_r_1_regridded['time'].values)
#ncum_r_3_regridded['time'] = pd.to_datetime(ncum_r_3_regridded['time'].values)
#ncum_r_5_regridded['time'] = pd.to_datetime(ncum_r_5_regridded['time'].values)
mean_obs = obs['rf'].mean(dim=('lat', 'lon')) 
mean_ncum_g_1 = ncum_g_1_regridded['APCP_surface'].mean(dim=('lat', 'lon'))
mean_ncum_g_3 = ncum_g_3_regridded['APCP_surface'].mean(dim=('lat', 'lon'))
mean_ncum_g_5 = ncum_g_5_regridded['APCP_surface'].mean(dim=('lat', 'lon'))
#mean_ncum_r_1 = ncum_r_1_regridded['APCP_24'].mean(dim=('lat', 'lon'))
#mean_ncum_r_1 = ncum_r_1_regridded['APCP_24'].mean(dim=('lat', 'lon'))
#mean_ncum_r_1 = ncum_r_1_regridded['APCP_24'].mean(dim=('lat', 'lon'))

# biases
bias_ncum_g_1 = mean_ncum_g_1 - mean_obs
bias_ncum_g_3 = mean_ncum_g_3 - mean_obs
bias_ncum_g_5 = mean_ncum_g_5 - mean_obs
#bias_ncum_r_1 = mean_ncum_r_1 - mean_obs
#bias_ncum_r_3 = mean_ncum_r_3 - mean_obs
#bias_ncum_r_5 = mean_ncum_r_5 - mean_obs

#yearly average
bias_ncum_g_1 = bias_ncum_g_1.groupby('time.year').mean()
bias_ncum_g_3 = bias_ncum_g_3.groupby('time.year').mean()
bias_ncum_g_5 = bias_ncum_g_5.groupby('time.year').mean()
#bias_ncum_r_1 = bias_ncum_r_1.groupby('time.year').mean()
#bias_ncum_r_3 = bias_ncum_r_3.groupby('time.year').mean()
#bias_ncum_r_5 = bias_ncum_r_5.groupby('time.year').mean()

#plot
plt.figure(figsize=(12,6))
plt.plot(bias_ncum_g_1['year'],bias_ncum_g_1, label='NCUM-G bias day 1',color='red')
plt.plot(bias_ncum_g_3['year'],bias_ncum_g_3, label='NCUM-G bias day 3',color='orange')
plt.plot(bias_ncum_g_5['year'],bias_ncum_g_5, label='NCUM-G bias day 5',color='yellow')
#plt.plot(bias_ncum_r_1['year'],bias_ncum_r_1, label='NCUM-R bias day 1',color='purple')
#plt.plot(bias_ncum_r_3['year'],bias_ncum_r_3, label='NCUM-R bias day 3',color='blue')
#plt.plot(bias_ncum_r_5['year'],bias_ncum_r_5, label='NCUM-R bias day 5',color='green')
plt.title('Comparison of Bias')
plt.xlabel('Time')
plt.xticks(bias_ncum_g_1['year'].values)
plt.ylabel('Mean Rainfall Bias(mm)')
plt.legend()
plt.grid(True)
plt.show()