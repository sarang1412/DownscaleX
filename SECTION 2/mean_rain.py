import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from shapely.geometry import mapping
from matplotlib.colors import LinearSegmentedColormap

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
EQM_ncum_1_regridded = EQM_ncum_1.interp_like(ncum_1, method='nearest')
EQM_ncum_3_regridded = EQM_ncum_3.interp_like(ncum_3, method='nearest')
EQM_ncum_5_regridded = EQM_ncum_5.interp_like(ncum_5, method='nearest')
PQM_ncum_1_regridded = PQM_ncum_1.interp_like(ncum_1, method='nearest')
PQM_ncum_3_regridded = PQM_ncum_3.interp_like(ncum_3, method='nearest')
PQM_ncum_5_regridded = PQM_ncum_5.interp_like(ncum_5, method='nearest')
GPQM_ncum_1_regridded = GPQM_ncum_1.interp_like(ncum_1, method='nearest')
GPQM_ncum_3_regridded = GPQM_ncum_3.interp_like(ncum_3, method='nearest')
GPQM_ncum_5_regridded = GPQM_ncum_5.interp_like(ncum_5, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106)) # capture full monsoon extend in bay of bengal
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

# mean rainfall
obs_mean = obs.mean(dim='time')
ncum_1_mean = ncum_1.mean(dim='time')
ncum_5_mean = ncum_5.mean(dim='time')
ncum_3_mean = ncum_3.mean(dim='time')
EQM_ncum_1_mean = EQM_ncum_1_regridded.mean(dim='time')
EQM_ncum_5_mean = EQM_ncum_5_regridded.mean(dim='time')
EQM_ncum_3_mean = EQM_ncum_3_regridded.mean(dim='time')
PQM_ncum_1_mean = PQM_ncum_1_regridded.mean(dim='time')
PQM_ncum_5_mean = PQM_ncum_5_regridded.mean(dim='time')
PQM_ncum_3_mean = PQM_ncum_3_regridded.mean(dim='time')
GPQM_ncum_1_mean = GPQM_ncum_1_regridded.mean(dim='time')
GPQM_ncum_5_mean = GPQM_ncum_5_regridded.mean(dim='time')
GPQM_ncum_3_mean = GPQM_ncum_3_regridded.mean(dim='time')

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs_mean = obs_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_1_mean = ncum_1_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_5_mean = ncum_5_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_3_mean = ncum_3_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
EQM_ncum_1_mean = EQM_ncum_1_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
EQM_ncum_5_mean = EQM_ncum_5_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
EQM_ncum_3_mean = EQM_ncum_3_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
PQM_ncum_1_mean = PQM_ncum_1_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
PQM_ncum_5_mean = PQM_ncum_5_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
PQM_ncum_3_mean = PQM_ncum_3_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
GPQM_ncum_1_mean = GPQM_ncum_1_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
GPQM_ncum_5_mean = GPQM_ncum_5_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
GPQM_ncum_3_mean = GPQM_ncum_3_mean.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

plt.figure(figsize=(20, 20))

# Create a custom colormap
hex_colors = ['#feecbe','#dcfecb','#95ff98','#64fffc','#04c4ff','#0066ff','#9364ff','#dc64ff','#ff01fe']
hexa = LinearSegmentedColormap.from_list('custom_gradient', hex_colors)

# Plot for Observed Mean Rainfall
ax1 = plt.subplot(4, 4, 1, projection=ccrs.PlateCarree())
obs_mean['rf'].plot(ax=ax1, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax1, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax1.set_title('Mean Rainfall Observation', fontsize=18)
#ax1.set_extent([65, 106, 6, 41], crs=ccrs.PlateCarree())

# Plot for NCUM-G Day 1 Mean Rainfall
ax2 = plt.subplot(4, 4, 2, projection=ccrs.PlateCarree())
ncum_1_mean['APCP_surface'].plot(ax=ax2, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax2, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax2.set_title('Mean Rainfall NCUM-G Day 1', fontsize=18)
#ax2.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

# Plot for NCUM-G Day 2 Mean Rainfall
ax3 = plt.subplot(4, 4, 3, projection=ccrs.PlateCarree())
ncum_3_mean['APCP_surface'].plot(ax=ax3, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax3.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax3, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax3.set_title('Mean Rainfall NCUM-G Day 3', fontsize=18)
#ax3.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

# Plot for NCUM-G Day 3 Mean Rainfall
ax4 = plt.subplot(4, 4, 4, projection=ccrs.PlateCarree())
ncum_5_mean['APCP_surface'].plot(ax=ax4, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax4.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax4, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax4.set_title('Mean Rainfall NCUM-G Day 5', fontsize=18)
#ax4.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

ax5 = plt.subplot(4, 4, 6, projection=ccrs.PlateCarree())
EQM_ncum_1_mean['RAINFALL'].plot(ax=ax5, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax5.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax5, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax5.set_title('Mean Rainfall NCUM-G Day 1 (EQM)', fontsize=18)
#ax5.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

# Plot for NCUM-G Day 2 Mean Rainfall
ax6 = plt.subplot(4, 4, 7, projection=ccrs.PlateCarree())
EQM_ncum_3_mean['RAINFALL'].plot(ax=ax6, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax6.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax6, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax6.set_title('Mean Rainfall NCUM-G Day 3 (EQM)', fontsize=18)
#ax6.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

# Plot for NCUM-G Day 3 Mean Rainfall
ax7 = plt.subplot(4, 4, 8, projection=ccrs.PlateCarree())
EQM_ncum_5_mean['RAINFALL'].plot(ax=ax7, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax7.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax7, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax7.set_title('Mean Rainfall NCUM-G Day 5 (EQM)', fontsize=18)
#ax7.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

ax8 = plt.subplot(4, 4, 10, projection=ccrs.PlateCarree())
EQM_ncum_1_mean['RAINFALL'].plot(ax=ax8, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax8.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax8, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax8.set_title('Mean Rainfall NCUM-G Day 1 (PQM)', fontsize=18)
#ax8.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

# Plot for NCUM-G Day 2 Mean Rainfall
ax9 = plt.subplot(4, 4, 11, projection=ccrs.PlateCarree())
EQM_ncum_3_mean['RAINFALL'].plot(ax=ax9, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax6.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax9, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax9.set_title('Mean Rainfall NCUM-G Day 3 (PQM)', fontsize=18)
#ax9.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

# Plot for NCUM-G Day 3 Mean Rainfall
ax9 = plt.subplot(4, 4, 12, projection=ccrs.PlateCarree())
EQM_ncum_5_mean['RAINFALL'].plot(ax=ax9, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax9.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax9, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax9.set_title('Mean Rainfall NCUM-G Day 5 (GPQM)', fontsize=18)
#ax9.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

ax10 = plt.subplot(4, 4, 14, projection=ccrs.PlateCarree())
EQM_ncum_1_mean['RAINFALL'].plot(ax=ax10, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax10.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax10, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax10.set_title('Mean Rainfall NCUM-G Day 1 (GPQM)', fontsize=18)
#ax10.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

# Plot for NCUM-G Day 2 Mean Rainfall
ax11 = plt.subplot(4, 4, 15, projection=ccrs.PlateCarree())
EQM_ncum_3_mean['RAINFALL'].plot(ax=ax11, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax11.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax11, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax11.set_title('Mean Rainfall NCUM-G Day 3 (GPQM)', fontsize=18)
#ax11.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

# Plot for NCUM-G Day 3 Mean Rainfall
ax12 = plt.subplot(4, 4, 16, projection=ccrs.PlateCarree())
EQM_ncum_5_mean['RAINFALL'].plot(ax=ax12, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'mean rainfall(mm/day)','shrink': 0.8}, vmin=0,vmax=52
)
ax12.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax12, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax12.set_title('Mean Rainfall NCUM-G Day 5 (GPQM)', fontsize=18)
#ax12.set_extent([88, 98, 21, 30], crs=ccrs.PlateCarree())

plt.tight_layout()
plt.show()