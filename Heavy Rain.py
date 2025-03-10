import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from shapely.geometry import mapping
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

#load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc', )
ncum_r_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2021-24.nc')
ncum_r_2 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day2rf_jjas2021-24.nc')
ncum_r_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day3rf_jjas2021-24.nc')
ncum_r_2
# Trim dataset
obs = obs.sel(time=slice('2021-06-01T03:00:00.000000000', '2024-09-30T03:00:00.000000000'))

#change lat and lon to match
ncum_r_1 = ncum_r_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_2 = ncum_r_2.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_3 = ncum_r_3.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_r_1_regridded = ncum_r_1.interp_like(obs, method='nearest')
ncum_r_3_regridded = ncum_r_3.interp_like(obs, method='nearest')
ncum_r_2_regridded = ncum_r_2.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))#62
ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_2_regridded = ncum_r_2_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_3_regridded = ncum_r_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))


# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs = obs.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r1 = ncum_r_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r2 = ncum_r_2_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r3 = ncum_r_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

obs_VH, obs_EX = obs.where((obs >= 115.6) & (obs <= 204.4), drop=True).fillna(0), obs.where(obs >= 204.5, drop=True).fillna(0)
ncum_r1_VH, ncum_r1_EX = ncum_r1.where((ncum_r1 >= 115.6) & (ncum_r1 <= 204.4), drop=True).fillna(0), ncum_r1.where(ncum_r1 >= 204.5, drop=True).fillna(0)
ncum_r2_VH, ncum_r2_EX = ncum_r2.where((ncum_r2 >= 115.6) & (ncum_r2 <= 204.4), drop=True).fillna(0), ncum_r2.where(ncum_r2 >= 204.5, drop=True).fillna(0)
ncum_r3_VH, ncum_r3_EX = ncum_r3.where((ncum_r3 >= 115.6) & (ncum_r3 <= 204.4), drop=True).fillna(0), ncum_r3.where(ncum_r3 >= 204.5, drop=True).fillna(0)

obs_VH_2D, obs_EX_2D = xr.where(obs_VH.any(dim='time'), 50, 0), xr.where(obs_EX.any(dim='time'), 100, 0)
ncum_r1_VH_2D, ncum_r1_EX_2D = xr.where(ncum_r1_VH.any(dim='time'), 50, 0), xr.where(ncum_r1_EX.any(dim='time'), 100, 0)
ncum_r2_VH_2D, ncum_r2_EX_2D = xr.where(ncum_r2_VH.any(dim='time'), 50, 0), xr.where(ncum_r2_EX.any(dim='time'), 100, 0)
ncum_r3_VH_2D, ncum_r3_EX_2D = xr.where(ncum_r3_VH.any(dim='time'), 50, 0), xr.where(ncum_r3_EX.any(dim='time'), 100, 0)

obs_final = np.sign(obs_EX_2D - obs_VH_2D)
ncum_r1_final = np.sign(ncum_r1_EX_2D - ncum_r1_VH_2D)
ncum_r2_final = np.sign(ncum_r2_EX_2D - ncum_r2_VH_2D)
ncum_r3_final = np.sign(ncum_r3_EX_2D - ncum_r3_VH_2D)

hex_colors = ['#64fffc','#04c4ff00','#ff01fe']
hexa = LinearSegmentedColormap.from_list('custom_gradient', hex_colors)



ax1 = plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
obs_final['rf'].plot(ax=ax1, cmap=hexa, transform=ccrs.PlateCarree(), add_colorbar=False)
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax1, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax1.set_title('Rainfall Intensity - Observation', fontsize=12)

# Plot for NCUM-G Day 1 Mean Rainfall
ax2 = plt.subplot(2, 2, 2, projection=ccrs.PlateCarree())
ncum_r3_final['APCP_surface'].plot(ax=ax2, cmap=hexa, transform=ccrs.PlateCarree(), add_colorbar=False)
ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax2, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax2.set_title('Rainfall Intensity - NCUM- Day 1', fontsize=12)

# Plot for NCUM-G Day 2 Mean Rainfall
ax3 = plt.subplot(2, 2, 3, projection=ccrs.PlateCarree())
ncum_r3_final['APCP_surface'].plot(ax=ax3, cmap=hexa, transform=ccrs.PlateCarree(), add_colorbar=False)
ax3.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax3, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax3.set_title('Rainfall Intensity - NCUM-R Day 2', fontsize=12)

# Plot for NCUM-G Day 3 Mean Rainfall
ax4 = plt.subplot(2, 2, 4, projection=ccrs.PlateCarree())
ncum_r3_final['APCP_surface'].plot(ax=ax4, cmap=hexa, transform=ccrs.PlateCarree(), add_colorbar=False)
ax4.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax4, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax4.set_title('Rainfall Intensity - NCUM-R Day 3', fontsize=12)

plt.tight_layout()
plt.show()