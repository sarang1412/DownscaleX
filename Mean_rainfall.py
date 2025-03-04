import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import dask.array as da
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
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
ncum_g_2 = ncum_g_2.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_3 = ncum_g_3.rename({'latitude': 'lat', 'longitude': 'lon'})
#ncum_r_1 = ncum_r_1.rename({'latitude': 'lat', 'longitude': 'lon'})
#ncum_r_3 = ncum_r_3.rename({'latitude': 'lat', 'longitude': 'lon'})
#ncum_r_5 = ncum_r_5.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_g_1_regridded = ncum_g_1.interp_like(obs, method='nearest')
ncum_g_2_regridded = ncum_g_2.interp_like(obs, method='nearest')
ncum_g_3_regridded = ncum_g_3.interp_like(obs, method='nearest')
#ncum_r_regridded_1 = ncum_r_1.interp_like(obs, method='nearest')
#ncum_r_regridded_3 = ncum_r_3.interp_like(obs, method='nearest')
#ncum_r_5_regridded = ncum_r_5.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))#62
ncum_g_1_regridded = ncum_g_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_2_regridded = ncum_g_2_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_3_regridded = ncum_g_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
#ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

# mean rainfall
obs_mean = obs.mean(dim='time')
ncum_g_1_mean = ncum_g_1_regridded.mean(dim='time')
ncum_g_2_mean = ncum_g_2_regridded.mean(dim='time')
ncum_g_3_mean = ncum_g_3_regridded.mean(dim='time')
#ncum_r_1_mean = ncum_r_1_regridded.mean(dim='time')
#ncum_r_3_mean = ncum_r_3_regridded.mean(dim='time')
#ncum_r_5_mean = ncum_r_5_regridded.mean(dim='time')

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs_mean = obs_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_1_mean = ncum_g_1_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_2_mean = ncum_g_2_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_3_mean = ncum_g_3_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

plt.figure(figsize=(14, 12))

# Plot for Observed Mean Rainfall
ax1 = plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
obs_mean['rf'].plot(ax=ax1, cmap='Blues', transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.9}, vmin=0,vmax=64
)
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax1, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax1.set_title('Mean Rainfall Observation', fontsize=12)

# Plot for NCUM-G Day 1 Mean Rainfall
ax2 = plt.subplot(2, 2, 2, projection=ccrs.PlateCarree())
ncum_g_1_mean['APCP_surface'].plot(ax=ax2, cmap='Blues', transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.9}, vmin=0,vmax=64
)
ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax2, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax2.set_title('Mean Rainfall NCUM-G Day 1', fontsize=12)

# Plot for NCUM-G Day 2 Mean Rainfall
ax3 = plt.subplot(2, 2, 3, projection=ccrs.PlateCarree())
ncum_g_2_mean['APCP_surface'].plot(ax=ax3, cmap='Blues', transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.9}, vmin=0,vmax=64
)
ax3.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax3, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax3.set_title('Mean Rainfall NCUM-G Day 2', fontsize=12)

# Plot for NCUM-G Day 3 Mean Rainfall
ax4 = plt.subplot(2, 2, 4, projection=ccrs.PlateCarree())
ncum_g_3_mean['APCP_surface'].plot(ax=ax4, cmap='Blues', transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.9}, vmin=0,vmax=64
)
ax4.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax4, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax4.set_title('Mean Rainfall NCUM-G Day 3', fontsize=12)

plt.tight_layout()
plt.show()