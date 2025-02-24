import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import rioxarray
from shapely.geometry import mapping

#load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc', )
ncum_g_1 = xr.open_dataset('e:\\Dissertation\\data\\ncum_day1rf_jjas2021-24.nc')
ncum_g_2 = xr.open_dataset('e:\\Dissertation\\data\\ncum_day2rf_jjas2021-24.nc')
ncum_g_3 = xr.open_dataset('e:\\Dissertation\\data\\ncum_day3rf_jjas2021-24.nc')
#ncum_r_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2023.nc')
#ncum_r_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day3rf_jjas2023.nc')
#ncum_r_5 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day5rf_jjas2023.nc')

# Trim dataset
obs = obs.sel(time=slice('2021-06-01T03:00:00.000000000', '2024-09-30T03:00:00.000000000'))


#change lat and lon to match
ncum_g_1 = ncum_g_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_3 = ncum_g_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_2 = ncum_g_2.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_g_1_regridded = ncum_g_1.interp_like(obs, method='nearest')
ncum_g_3_regridded = ncum_g_3.interp_like(obs, method='nearest')
ncum_g_2_regridded = ncum_g_2.interp_like(obs, method='nearest')
#ncum_r_1_regridded = ncum_r_1.interp_like(obs, method='nearest')
#ncum_r_3_regridded = ncum_r_3.interp_like(obs, method='nearest')
#ncum_r_5_regridded = ncum_r_5.interp_like(obs, method='nearest')

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs = obs.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_1_regridded = ncum_g_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_3_regridded = ncum_g_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_2_regridded = ncum_g_2_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_1_regridded = ncum_g_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_3_regridded = ncum_g_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_2_regridded = ncum_g_2_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_3_regridded = ncum_r_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_5_regridded = ncum_r_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

#proper date format
obs['time'] = pd.to_datetime(obs['time'].values)
ncum_g_1_regridded['time'] = pd.to_datetime(ncum_g_1_regridded['time'].values)
ncum_g_3_regridded['time'] = pd.to_datetime(ncum_g_3_regridded['time'].values)
ncum_g_2_regridded['time'] = pd.to_datetime(ncum_g_2_regridded['time'].values)
#ncum_r_1_regridded['time'] = pd.to_datetime(ncum_r_1_regridded['time'].values)
#ncum_r_3_regridded['time'] = pd.to_datetime(ncum_r_3_regridded['time'].values)
#ncum_r_5_regridded['time'] = pd.to_datetime(ncum_r_5_regridded['time'].values)
mean_obs = obs['rf'].mean(dim=('lat', 'lon')) 
mean_ncum_g_1 = ncum_g_1_regridded['APCP_surface'].mean(dim=('lat', 'lon'))
mean_ncum_g_3 = ncum_g_3_regridded['APCP_surface'].mean(dim=('lat', 'lon'))
mean_ncum_g_2 = ncum_g_2_regridded['APCP_surface'].mean(dim=('lat', 'lon'))
#mean_ncum_r_1 = ncum_r_1_regridded['APCP_24'].mean(dim=('lat', 'lon'))
#mean_ncum_r_1 = ncum_r_1_regridded['APCP_24'].mean(dim=('lat', 'lon'))
#mean_ncum_r_1 = ncum_r_1_regridded['APCP_24'].mean(dim=('lat', 'lon'))

# biases
bias_ncum_g_1 = mean_ncum_g_1 - mean_obs
bias_ncum_g_3 = mean_ncum_g_3 - mean_obs
bias_ncum_g_2 = mean_ncum_g_2 - mean_obs
#bias_ncum_r_1 = mean_ncum_r_1 - mean_obs
#bias_ncum_r_3 = mean_ncum_r_3 - mean_obs
#bias_ncum_r_5 = mean_ncum_r_5 - mean_obs

#yearly average
bias_ncum_g_1 = bias_ncum_g_1.groupby('time.year').mean()
bias_ncum_g_3 = bias_ncum_g_3.groupby('time.year').mean()
bias_ncum_g_2 = bias_ncum_g_2.groupby('time.year').mean()
#bias_ncum_r_1 = bias_ncum_r_1.groupby('time.year').mean()
#bias_ncum_r_3 = bias_ncum_r_3.groupby('time.year').mean()
#bias_ncum_r_5 = bias_ncum_r_5.groupby('time.year').mean()

#plot
plt.figure(figsize=(12,6))
plt.plot(bias_ncum_g_1['year'],bias_ncum_g_1, label='NCUM-G bias day 1',color='red')
plt.plot(bias_ncum_g_2['year'],bias_ncum_g_2, label='NCUM-G bias day 2',color='orange')
plt.plot(bias_ncum_g_3['year'],bias_ncum_g_3, label='NCUM-G bias day 3',color='yellow')
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