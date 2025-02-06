import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import dask.array as da
import geopandas as gpd
import rioxarray
from shapely.geometry import mapping

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

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs = obs.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_1_regridded = ncum_g_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_3_regridded = ncum_g_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_5_regridded = ncum_g_5_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

#proper date format
obs['time'] = pd.to_datetime(obs['time'].values)
ncum_g_1_regridded['time'] = pd.to_datetime(ncum_g_1_regridded['time'].values)
ncum_g_3_regridded['time'] = pd.to_datetime(ncum_g_3_regridded['time'].values)
ncum_g_5_regridded['time'] = pd.to_datetime(ncum_g_5_regridded['time'].values)
#ncum_r_1_regridded['time'] = pd.to_datetime(ncum_r_1_regridded['time'].values)
#ncum_r_3_regridded['time'] = pd.to_datetime(ncum_r_3_regridded['time'].values)
#ncum_r_5_regridded['time'] = pd.to_datetime(ncum_r_5_regridded['time'].values)

# daily averages for each day of JJAS months over the years 2020-2024
daily_avg_obs = obs.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_g_1 = ncum_g_1_regridded.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_g_3 = ncum_g_3_regridded.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_g_5 = ncum_g_5_regridded.groupby('time.dayofyear').mean(dim='time')
#daily_avg_ncum_r_1 = ncum_r_1_regridded.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_g_1 = daily_avg_ncum_g_1.dropna(dim='dayofyear', how='all')
daily_avg_ncum_g_3 = daily_avg_ncum_g_3.dropna(dim='dayofyear', how='all')
daily_avg_ncum_g_5 = daily_avg_ncum_g_5.dropna(dim='dayofyear', how='all')

#variance
obs_variance = daily_avg_obs['rf'].var(dim=('lat', 'lon'))
ncum_g_1_variance = daily_avg_ncum_g_1['APCP_surface'].var(dim=('lat', 'lon'))
ncum_g_3_variance = daily_avg_ncum_g_3['APCP_surface'].var(dim=('lat', 'lon'))
ncum_g_5_variance = daily_avg_ncum_g_5['APCP_surface'].var(dim=('lat', 'lon'))
#ncum_r_1_variance = daily_avg_ncum_r_1['APCP_24'].var(dim=('lat', 'lon'))
#ncum_r_3_variance = daily_avg_ncum_r_3['APCP_24'].var(dim=('lat', 'lon'))
#ncum_r_5_variance = daily_avg_ncum_r_5['APCP_24'].var(dim=('lat', 'lon'))

# Plot
data_variances = [
    (obs_variance, 'Observed data 10Km'),
    (ncum_g_1_variance, 'NCUM-G 12Km to 10Km (day 1 forecast)'),
    (ncum_g_3_variance, 'NCUM-G 12Km to 10Km (day 3 forecast)'),
    (ncum_g_5_variance, 'NCUM-G 12Km to 10Km (day 5 forecast)')#,
    #(ncum_r_1_variance, 'NCUM-R 4km to 10Km (day 1 forecast)'),
    #(ncum_r_3_variance, 'NCUM-R 4km to 10Km (day 3 forecast)'),
    #(ncum_r_5_variance, 'NCUM-R 4km to 10Km (day 5 forecast)')
]

plt.plot(ncum_g_1_variance['dayofyear'], obs_variance, label="Observed data 10Km",color = "black")
plt.plot(ncum_g_1_variance['dayofyear'], ncum_g_1_variance, label="NCUM-G 12Km to 10Km (day 1 forecast)",color = "red")
plt.plot(ncum_g_1_variance['dayofyear'], ncum_g_3_variance, label="NCUM-G 12Km to 10Km (day 3 forecast)",color = "blue")
plt.plot(ncum_g_1_variance['dayofyear'], ncum_g_5_variance, label="NCUM-G 12Km to 10Km (day 5 forecast)",color = "green")
plt.xlabel('Time')
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%B'))
plt.ylabel('Variance')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()