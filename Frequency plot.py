import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from shapely.geometry import mapping
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from matplotlib.colors import ListedColormap

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




# Compute frequency of exceedance
threshold = 64.5
obs_freq = (obs['rf'] > threshold).sum(dim='time')  # Count occurrences over time
ncum_r1_freq = (ncum_r1['APCP_surface'] > threshold).sum(dim='time')
ncum_r2_freq = (ncum_r2['APCP_surface'] > threshold).sum(dim='time')
ncum_r3_freq = (ncum_r3['APCP_surface'] > threshold).sum(dim='time')



# Create a custom colormap
hex_colors = ['#04c4ff00', '#a100c9', '#1e3cff','#00c8c8','#00dc00','#f08228','#fa3c3c']
hexa = ListedColormap(hex_colors)

# Plot function
def plot_rainfall_freq(data, title, ax):
    img = ax.pcolormesh(data['lon'], data['lat'], data, cmap=hexa, vmin=0, vmax=70, transform=ccrs.PlateCarree())
    shape.boundary.plot(ax=ax, edgecolor='black', linewidth=0.7, transform=ccrs.PlateCarree())
    ax.set_title(title, fontsize=12)
    return img

fig, axes = plt.subplots(2, 2, figsize=(16, 12), subplot_kw={'projection': ccrs.PlateCarree()})

img1 = plot_rainfall_freq(obs_freq, "Observed Very Heavy Rainfall > 64.5 mm", axes[0, 0])
img2 = plot_rainfall_freq(ncum_r1_freq, "NCUM R Day 1 Very Heavy Rainfall > 64.5 mm", axes[0, 1])
img3 = plot_rainfall_freq(ncum_r2_freq, "NCUM R Day 2 Very Heavy Rainfall > 64.5 mm", axes[1, 0])
img4 = plot_rainfall_freq(ncum_r3_freq, "NCUM R Day 3 Very Heavy Rainfall > 64.5 mm", axes[1, 1])

cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(img1, cax=cbar_ax)
cbar.set_label("Frequency of Rainfall > 64.5 mm")
plt.tight_layout()
plt.show()