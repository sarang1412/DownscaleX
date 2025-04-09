import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from shapely.geometry import mapping
from matplotlib.colors import LinearSegmentedColormap

# Load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc')

# Trim dataset
obs = obs.sel(time=slice('2021-06-01T03:00:00.000000000', '2024-09-30T03:00:00.000000000'))

# Lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))

# Use the first timestep of rf
obs_first_timestep = obs.isel(time=0)

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)

obs_first_timestep = obs_first_timestep.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

plt.figure(figsize=(7, 6))

# Create a custom colormap
hex_colors = ['#feecbe','#dcfecb','#95ff98','#64fffc','#04c4ff','#0066ff','#9364ff','#dc64ff','#ff01fe']
hexa = LinearSegmentedColormap.from_list('custom_gradient', hex_colors)

# Plot for Observed Rainfall (First Timestep)
ax1 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
obs_first_timestep['rf'].plot(ax=ax1, cmap=hexa, transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'rainfall(mm/day)', 'shrink': 0.9}, vmin=0, vmax=52)
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
shape.boundary.plot(ax=ax1, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
ax1.set_title('Rainfall Observation (First Timestep)', fontsize=12)

plt.tight_layout()
plt.show()