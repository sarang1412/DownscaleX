import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import rioxarray
from shapely.geometry import mapping

#load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc', )
ncum_r_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2021-24.nc')
ncum_r_2 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day2rf_jjas2021-24.nc')
ncum_r_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day3rf_jjas2021-24.nc')

# Trim dataset
obs = obs.sel(time=slice('2021-06-01T03:00:00.000000000', '2024-09-30T03:00:00.000000000'))


#change lat and lon to match
ncum_r_1 = ncum_r_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_3 = ncum_r_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_2 = ncum_r_2.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_r_1_regridded = ncum_r_1.interp_like(obs, method='nearest')
ncum_r_3_regridded = ncum_r_3.interp_like(obs, method='nearest')
ncum_r_2_regridded = ncum_r_2.interp_like(obs, method='nearest')

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs = obs.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r_1_regridded = ncum_r_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r_3_regridded = ncum_r_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r_2_regridded = ncum_r_2_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_3_regridded = ncum_r_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_2_regridded = ncum_r_2_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

#proper date format
obs['time'] = pd.to_datetime(obs['time'].values)
ncum_r_1_regridded['time'] = pd.to_datetime(ncum_r_1_regridded['time'].values)
ncum_r_3_regridded['time'] = pd.to_datetime(ncum_r_3_regridded['time'].values)
ncum_r_2_regridded['time'] = pd.to_datetime(ncum_r_2_regridded['time'].values)

# Filter data for threshold 0 - 30 mm
obs_filtered = obs.where((obs['rf'] >= 0) & (obs['rf'] <= 30), drop=True)
ncum_r_1_filtered = ncum_r_1_regridded.where((ncum_r_1_regridded['APCP_surface'] >= 0) & (ncum_r_1_regridded['APCP_surface'] <= 30), drop=True)
ncum_r_3_filtered = ncum_r_3_regridded.where((ncum_r_3_regridded['APCP_surface'] >= 0) & (ncum_r_3_regridded['APCP_surface'] <= 30), drop=True)
ncum_r_2_filtered = ncum_r_2_regridded.where((ncum_r_2_regridded['APCP_surface'] >= 0) & (ncum_r_2_regridded['APCP_surface'] <= 30), drop=True)

# Calculate mean
mean_obs = obs_filtered['rf'].mean(dim=('lat', 'lon'))
mean_ncum_r_1 = ncum_r_1_filtered['APCP_surface'].mean(dim=('lat', 'lon'))
mean_ncum_r_3 = ncum_r_3_filtered['APCP_surface'].mean(dim=('lat', 'lon'))
mean_ncum_r_2 = ncum_r_2_filtered['APCP_surface'].mean(dim=('lat', 'lon'))

# biases
bias_ncum_r_1 = mean_ncum_r_1 - mean_obs
bias_ncum_r_3 = mean_ncum_r_3 - mean_obs
bias_ncum_r_2 = mean_ncum_r_2 - mean_obs

#yearly average
bias_ncum_r_1 = bias_ncum_r_1.groupby('time.year').mean()
bias_ncum_r_3 = bias_ncum_r_3.groupby('time.year').mean()
bias_ncum_r_2 = bias_ncum_r_2.groupby('time.year').mean()

#plot
plt.figure(figsize=(12,6))

plt.xlabel('Time')
plt.xticks(bias_ncum_r_1['year'].values)
plt.ylabel('Mean Rainfall Bias(mm)')
plt.xlabel('Time')
plt.title('Comparison of Bias')
plt.plot(bias_ncum_r_1['year'],bias_ncum_r_1, label='NCUM-R bias day 1',color='purple')
plt.plot(bias_ncum_r_2['year'],bias_ncum_r_2, label='NCUM-R bias day 2',color='blue')
plt.plot(bias_ncum_r_3['year'],bias_ncum_r_3, label='NCUM-R bias day 3',color='green')
plt.ylabel('Mean Rainfall Bias(mm)')
plt.xticks(bias_ncum_g_1['year'].values)
plt.legend()
plt.grid(True)
plt.show()