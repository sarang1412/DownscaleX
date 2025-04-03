import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import mapping
import xarray as xr 
import numpy as np 
import matplotlib.dates as mdates

#load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc')
ncum_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc')
ncum_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day3rf_jjas2020-24.nc')
ncum_5 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day5rf_jjas2020-24.nc')
EQM_ncum_1 = xr.open_dataset('e:\\Dissertation\\data\\EQM_ncumg_day1rf_jjas2020-24-0p25.nc')
EQM_ncum_3 = xr.open_dataset('e:\\Dissertation\\data\\EQM_ncumg_day3rf_jjas2020-24-0p25.nc')
EQM_ncum_5 = xr.open_dataset('e:\\Dissertation\\data\\EQM_ncumg_day5rf_jjas2020-24-0p25.nc')
PQM_ncum_1 = xr.open_dataset('e:\\Dissertation\\data\\PQM_ncumg_day1rf_jjas2020-24-0p25.nc')
PQM_ncum_3 = xr.open_dataset('e:\\Dissertation\\data\\PQM_ncumg_day3rf_jjas2020-24-0p25.nc')
PQM_ncum_5 = xr.open_dataset('e:\\Dissertation\\data\\PQM_ncumg_day5rf_jjas2020-24-0p25.nc')
GPQM_ncum_1 = xr.open_dataset('e:\\Dissertation\\data\\GPQM_ncumg_day1rf_jjas2020-24-0p25.nc')
GPQM_ncum_3 = xr.open_dataset('e:\\Dissertation\\data\\GPQM_ncumg_day3rf_jjas2020-24-0p25.nc')
GPQM_ncum_5 = xr.open_dataset('e:\\Dissertation\\data\\GPQM_ncumg_day5rf_jjas2020-24-0p25.nc')

# Trim dataset
#obs = obs.sel(time=slice('2021-06-01T03:00:00.000000000', '2024-09-30T03:00:00.000000000'))

#change lat and lon to match
ncum_1 = ncum_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_3 = ncum_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_5 = ncum_5.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
obs = obs.interp_like(ncum_1, method='nearest')
EQM_ncum_1_regridded = EQM_ncum_1.interp_like(obs, method='nearest')
EQM_ncum_3_regridded = EQM_ncum_3.interp_like(obs, method='nearest')
EQM_ncum_5_regridded = EQM_ncum_5.interp_like(obs, method='nearest')
PQM_ncum_1_regridded = PQM_ncum_1.interp_like(obs, method='nearest')
PQM_ncum_3_regridded = PQM_ncum_3.interp_like(obs, method='nearest')
PQM_ncum_5_regridded = PQM_ncum_5.interp_like(obs, method='nearest')
GPQM_ncum_1_regridded = GPQM_ncum_1.interp_like(obs, method='nearest')
GPQM_ncum_3_regridded = GPQM_ncum_3.interp_like(obs, method='nearest')
GPQM_ncum_5_regridded = GPQM_ncum_5.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))#62
ncum_1 = ncum_1.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_5 = ncum_5.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_3 = ncum_3.sel(lat=slice(6, 41), lon=slice(65, 106))
EQM_ncum_1_regridded = EQM_ncum_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
EQM_ncum_5_regridded = EQM_ncum_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
EQM_ncum_3_regridded = EQM_ncum_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
PQM_ncum_1_regridded = PQM_ncum_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
PQM_ncum_5_regridded = PQM_ncum_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
PQM_ncum_3_regridded = PQM_ncum_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
GPQM_ncum_1_regridded = GPQM_ncum_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
GPQM_ncum_5_regridded = GPQM_ncum_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
GPQM_ncum_3_regridded = GPQM_ncum_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs_data = obs.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_1_data = ncum_1.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_5_data = ncum_5.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_3_data = ncum_3.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
EQM_ncum_1_data = EQM_ncum_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
EQM_ncum_5_data = EQM_ncum_5_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
EQM_ncum_3_data = EQM_ncum_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
PQM_ncum_1_data = PQM_ncum_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
PQM_ncum_5_data = PQM_ncum_5_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
PQM_ncum_3_data = PQM_ncum_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
GPQM_ncum_1_data = GPQM_ncum_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
GPQM_ncum_5_data = GPQM_ncum_5_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
GPQM_ncum_3_data = GPQM_ncum_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

# daily averages for each day of JJAS months over the years 2020-2024
daily_avg_obs = obs_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_1 = ncum_1_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_3 = ncum_3_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_5 = ncum_5_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_EQM_ncum_1 = EQM_ncum_1_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_EQM_ncum_3 = EQM_ncum_3_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_EQM_ncum_5 = EQM_ncum_5_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_PQM_ncum_1 = PQM_ncum_1_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_PQM_ncum_3 = PQM_ncum_3_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_PQM_ncum_5 = PQM_ncum_5_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_GPQM_ncum_1 = GPQM_ncum_1_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_GPQM_ncum_3 = GPQM_ncum_3_data.groupby('time.dayofyear').mean(dim='time')
daily_avg_GPQM_ncum_5 = GPQM_ncum_5_data.groupby('time.dayofyear').mean(dim='time')

#drom NaN values
daily_avg_ncum_1 = daily_avg_ncum_1.dropna(dim='dayofyear', how='all')
daily_avg_ncum_3 = daily_avg_ncum_3.dropna(dim='dayofyear', how='all')
daily_avg_ncum_5 = daily_avg_ncum_5.dropna(dim='dayofyear', how='all')
daily_avg_EQM_ncum_1 = daily_avg_EQM_ncum_1.dropna(dim='dayofyear', how='all')
daily_avg_EQM_ncum_3 = daily_avg_EQM_ncum_3.dropna(dim='dayofyear', how='all')
daily_avg_EQM_ncum_5 = daily_avg_EQM_ncum_5.dropna(dim='dayofyear', how='all') 
daily_avg_PQM_ncum_1 = daily_avg_PQM_ncum_1.dropna(dim='dayofyear', how='all') 
daily_avg_PQM_ncum_3 = daily_avg_PQM_ncum_3.dropna(dim='dayofyear', how='all') 
daily_avg_PQM_ncum_5 = daily_avg_PQM_ncum_5.dropna(dim='dayofyear', how='all') 
daily_avg_GPQM_ncum_1 = daily_avg_GPQM_ncum_1.dropna(dim='dayofyear', how='all')
daily_avg_GPQM_ncum_3 = daily_avg_GPQM_ncum_3.dropna(dim='dayofyear', how='all')
daily_avg_GPQM_ncum_5 = daily_avg_GPQM_ncum_5.dropna(dim='dayofyear', how='all')

#variance
obs_variance = daily_avg_obs['rf'].var(dim=('lat', 'lon'))
ncum_1_variance = daily_avg_ncum_1['APCP_surface'].var(dim=('lat', 'lon'))
ncum_3_variance = daily_avg_ncum_3['APCP_surface'].var(dim=('lat', 'lon'))
ncum_5_variance = daily_avg_ncum_5['APCP_surface'].var(dim=('lat', 'lon'))
ncum_1_EQM_variance = daily_avg_EQM_ncum_1['RAINFALL'].var(dim=('lat', 'lon'))
ncum_3_EQM_variance = daily_avg_EQM_ncum_3['RAINFALL'].var(dim=('lat', 'lon'))
ncum_5_EQM_variance = daily_avg_EQM_ncum_5['RAINFALL'].var(dim=('lat', 'lon'))
ncum_1_PQM_variance = daily_avg_PQM_ncum_1['RAINFALL'].var(dim=('lat', 'lon'))
ncum_3_PQM_variance = daily_avg_PQM_ncum_3['RAINFALL'].var(dim=('lat', 'lon'))
ncum_5_PQM_variance = daily_avg_PQM_ncum_5['RAINFALL'].var(dim=('lat', 'lon'))
ncum_1_GPQM_variance = daily_avg_GPQM_ncum_1['RAINFALL'].var(dim=('lat', 'lon'))
ncum_3_GPQM_variance = daily_avg_GPQM_ncum_3['RAINFALL'].var(dim=('lat', 'lon'))
ncum_5_GPQM_variance = daily_avg_GPQM_ncum_5['RAINFALL'].var(dim=('lat', 'lon'))

# Create a 2x2 subplot layout
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# No Bias Correction
axs[0, 0].plot(ncum_1_variance['dayofyear'], obs_variance, label="Observed data 10Km", color="black")
axs[0, 0].plot(ncum_1_variance['dayofyear'], ncum_1_variance, label="NCUM-G 12Km to 10Km (day 1 forecast)", color="red")
axs[0, 0].plot(ncum_1_variance['dayofyear'], ncum_3_variance, label="NCUM-G 12Km to 10Km (day 2 forecast)", color="blue")
axs[0, 0].plot(ncum_1_variance['dayofyear'], ncum_5_variance, label="NCUM-G 12Km to 10Km (day 3 forecast)", color="green")
axs[0, 0].set_title("No Bias Correction")
axs[0, 0].set_xlabel('Time')
axs[0, 0].set_ylabel('Variance')
axs[0, 0].xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
axs[0, 0].xaxis.set_major_formatter(mdates.DateFormatter('%B'))
axs[0, 0].legend()
axs[0, 0].grid(True)

# EQM
axs[0, 1].plot(ncum_1_variance['dayofyear'], obs_variance, label="Observed data 10Km", color="black")
axs[0, 1].plot(ncum_1_variance['dayofyear'], ncum_1_EQM_variance, label="NCUM-G EQM (day 1 forecast)", color="red")
axs[0, 1].plot(ncum_1_variance['dayofyear'], ncum_3_EQM_variance, label="NCUM-G EQM (day 2 forecast)", color="blue")
axs[0, 1].plot(ncum_1_variance['dayofyear'], ncum_5_EQM_variance, label="NCUM-G EQM (day 3 forecast)", color="green")
axs[0, 1].set_title("EQM")
axs[0, 1].set_xlabel('Time')
axs[0, 1].set_ylabel('Variance')
axs[0, 1].xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
axs[0, 1].xaxis.set_major_formatter(mdates.DateFormatter('%B'))
axs[0, 1].legend()
axs[0, 1].grid(True)

# PQM
axs[1, 0].plot(ncum_1_variance['dayofyear'], obs_variance, label="Observed data 10Km", color="black")
axs[1, 0].plot(ncum_1_variance['dayofyear'], ncum_1_PQM_variance, label="NCUM-G PQM (day 1 forecast)", color="red")
axs[1, 0].plot(ncum_1_variance['dayofyear'], ncum_3_PQM_variance, label="NCUM-G PQM (day 2 forecast)", color="blue")
axs[1, 0].plot(ncum_1_variance['dayofyear'], ncum_5_PQM_variance, label="NCUM-G PQM (day 3 forecast)", color="green")
axs[1, 0].set_title("PQM")
axs[1, 0].set_xlabel('Time')
axs[1, 0].set_ylabel('Variance')
axs[1, 0].xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
axs[1, 0].xaxis.set_major_formatter(mdates.DateFormatter('%B'))
axs[1, 0].legend()
axs[1, 0].grid(True)

# GPQM
axs[1, 1].plot(ncum_1_variance['dayofyear'], obs_variance, label="Observed data 10Km", color="black")
axs[1, 1].plot(ncum_1_variance['dayofyear'], ncum_1_GPQM_variance, label="NCUM-G GPQM (day 1 forecast)", color="red")
axs[1, 1].plot(ncum_1_variance['dayofyear'], ncum_3_GPQM_variance, label="NCUM-G GPQM (day 2 forecast)", color="blue")
axs[1, 1].plot(ncum_1_variance['dayofyear'], ncum_5_GPQM_variance, label="NCUM-G GPQM (day 3 forecast)", color="green")
axs[1, 1].set_title("GPQM")
axs[1, 1].set_xlabel('Time')
axs[1, 1].set_ylabel('Variance')
axs[1, 1].xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
axs[1, 1].xaxis.set_major_formatter(mdates.DateFormatter('%B'))
axs[1, 1].legend()
axs[1, 1].grid(True)

# Adjust layout and show the plot
plt.tight_layout()
plt.show()