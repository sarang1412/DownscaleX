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
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc')
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

# daily averages for each day of JJAS months over the years 2020-2024
daily_avg_obs = obs.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_r_1 = ncum_r_1_regridded.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_r_3 = ncum_r_3_regridded.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_r_2 = ncum_r_2_regridded.groupby('time.dayofyear').mean(dim='time')

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

daily_avg_obs = daily_avg_obs.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
daily_avg_ncum_r_1 = daily_avg_ncum_r_1.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
daily_avg_ncum_r_2 = daily_avg_ncum_r_2.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
daily_avg_ncum_r_3 = daily_avg_ncum_r_3.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

#drom NaN values
daily_avg_ncum_r_1 = daily_avg_ncum_r_1.dropna(dim='dayofyear', how='all')
daily_avg_ncum_r_3 = daily_avg_ncum_r_3.dropna(dim='dayofyear', how='all')
daily_avg_ncum_r_2 = daily_avg_ncum_r_2.dropna(dim='dayofyear', how='all')

#variance
obs_variance = daily_avg_obs['rf'].var(dim=('lat', 'lon'))
ncum_r_1_variance = daily_avg_ncum_r_1['APCP_surface'].var(dim=('lat', 'lon'))
ncum_r_3_variance = daily_avg_ncum_r_3['APCP_surface'].var(dim=('lat', 'lon'))
ncum_r_2_variance = daily_avg_ncum_r_2['APCP_surface'].var(dim=('lat', 'lon'))

# Plot
data_variances = [
    (obs_variance, 'Observed data 10Km'),
    (ncum_r_1_variance, 'NCUM-R 4km to 10Km (day 1 forecast)'),
    (ncum_r_2_variance, 'NCUM-R 4km to 10Km (day 2 forecast)'),
    (ncum_r_3_variance, 'NCUM-R 4km to 10Km (day 3 forecast)')
]

plt.plot(ncum_r_1_variance['dayofyear'], obs_variance, label="Observed data 10Km",color = "black")
plt.plot(ncum_r_1_variance['dayofyear'], ncum_r_1_variance, label="NCUM-R 12Km to 10Km (day 1 forecast)",color = "red")
plt.plot(ncum_r_1_variance['dayofyear'], ncum_r_2_variance, label="NCUM-R 12Km to 10Km (day 2 forecast)",color = "blue")
plt.plot(ncum_r_1_variance['dayofyear'], ncum_r_3_variance, label="NCUM-R 12Km to 10Km (day 3 forecast)",color = "green")
plt.xlabel('Time')
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%B'))
plt.ylabel('Variance')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()