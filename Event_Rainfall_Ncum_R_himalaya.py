import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from shapely.geometry import mapping
from matplotlib.colors import LinearSegmentedColormap

#load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc', )
ncum_r_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2021-24.nc')
ncum_r_2 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day2rf_jjas2021-24.nc')
ncum_r_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day3rf_jjas2021-24.nc')

#change lat and lon to match
ncum_r_1 = ncum_r_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_2 = ncum_r_2.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_3 = ncum_r_3.rename({'latitude': 'lat', 'longitude': 'lon'})

# Trim dataset
obs = obs.sel(time=slice('2021-06-01T03:00:00.000000000', '2024-09-30T03:00:00.000000000'))

#regrid
ncum_r_1_regridded = ncum_r_1#.interp_like(obs, method='nearest')
ncum_r_3_regridded = ncum_r_3#.interp_like(obs, method='nearest')
ncum_r_2_regridded = ncum_r_2#.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))#62
ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_2_regridded = ncum_r_2_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_3_regridded = ncum_r_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

# Select event
t = "2023-07-10T0:00:00.000000000"
obs_event = obs.sel(time= t, method='nearest')
ncum_r_day_1 = ncum_r_1_regridded.sel(time=t, method='nearest')
ncum_r_day_3 = ncum_r_3_regridded.sel(time=t, method='nearest')
ncum_r_day_2 = ncum_r_2_regridded.sel(time=t, method='nearest')

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
shape = shape.set_crs(epsg=4326, allow_override=True)
shape = shape.to_crs(epsg=4326)

obs_event = obs_event.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r_1_event = ncum_r_day_1.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r_2_event = ncum_r_day_2.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r_3_event = ncum_r_day_3.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

plt.figure(figsize=(14, 12))
hex_colors = ['#feecbe','#dcfecb','#95ff98','#64fffc','#04c4ff','#0066ff','#9364ff','#dc64ff','#ff01fe']
hexa = LinearSegmentedColormap.from_list('custom_gradient', hex_colors)
# Plot for Observed Rainfall
city = "Dehradun"
clon = 78.03
clat = 30.32
extent = [74, 79, 28, 33] 
ax1 = plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
obs_event['rf'].plot(ax=ax1, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'rainfall(mm/day)','shrink': 0.9}, vmin=0,vmax=120 
)
ax1.set_extent(extent, crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax1.plot(clon, clat, "ko", transform=ccrs.PlateCarree())
ax1.text(clon + 0.03, clat, city, fontsize=12, transform=ccrs.PlateCarree())
shape.boundary.plot(ax=ax1, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax1.set_title('Event Rainfall - Observation', fontsize=12)

# Plot for NCUM-G Day 1
ax2 = plt.subplot(2, 2, 2, projection=ccrs.PlateCarree())
ncum_r_1_event['APCP_surface'].plot(ax=ax2, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'rainfall(mm/day)','shrink': 0.9}, vmin=0,vmax=120
)
ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax2.set_extent(extent, crs=ccrs.PlateCarree())
ax2.plot(clon, clat, "ko", transform=ccrs.PlateCarree())
ax2.text(clon + 0.03, clat, city, fontsize=12, transform=ccrs.PlateCarree())
shape.boundary.plot(ax=ax2, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax2.set_title('Event Rainfall - NCUM-R Day 1', fontsize=12)

# Plot for NCUM-G Day 2
ax3 = plt.subplot(2, 2, 3, projection=ccrs.PlateCarree())
ncum_r_2_event['APCP_surface'].plot(ax=ax3, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'rainfall(mm/day)','shrink': 0.9}, vmin=0,vmax=120
)
ax3.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax3.set_extent(extent, crs=ccrs.PlateCarree())
ax3.plot(clon, clat, "ko", transform=ccrs.PlateCarree())
ax3.text(clon + 0.03, clat, city, fontsize=12, transform=ccrs.PlateCarree())
shape.boundary.plot(ax=ax3, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax3.set_title('Event Rainfall - NCUM-R Day 2', fontsize=12)

# Plot for NCUM-G Day 3
ax4 = plt.subplot(2, 2, 4, projection=ccrs.PlateCarree())
ncum_r_3_event['APCP_surface'].plot(ax=ax4, cmap=hexa, transform=ccrs.PlateCarree(),cbar_kwargs={'label': 'rainfall(mm/day)','shrink': 0.9}, vmin=0,vmax=120
)
ax4.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax4.set_extent(extent, crs=ccrs.PlateCarree())
ax4.plot(clon, clat, "ko", transform=ccrs.PlateCarree())
ax4.text(clon + 0.03, clat, city, fontsize=12, transform=ccrs.PlateCarree())
shape.boundary.plot(ax=ax4, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax4.set_title('Event Rainfall - NCUM-R Day 3', fontsize=12)

plt.tight_layout()
plt.show()