import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import geopandas as gpd
from shapely.geometry import mapping
from matplotlib.colors import LinearSegmentedColormap

# Load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc')
ncum_r_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2021-24.nc')
ncum_r_2 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day2rf_jjas2021-24.nc')
ncum_r_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day3rf_jjas2021-24.nc')

# trimming
obs = obs.sel(time=slice('2021-06-01T03:00:00.000000000', '2024-09-30T03:00:00.000000000'))

# Rename coordinates
ncum_r_1 = ncum_r_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_2 = ncum_r_2.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_3 = ncum_r_3.rename({'latitude': 'lat', 'longitude': 'lon'})


# Regrid
ncum_r_1_regridded = ncum_r_1.interp_like(obs, method='nearest')
ncum_r_2_regridded = ncum_r_2.interp_like(obs, method='nearest')
ncum_r_3_regridded = ncum_r_3.interp_like(obs, method='nearest')

# Crop
lat_range, lon_range = slice(6, 41), slice(65, 106)
obs = obs.sel(lat=lat_range, lon=lon_range)
ncum_r_1_regridded = ncum_r_1_regridded.sel(lat=lat_range, lon=lon_range)
ncum_r_2_regridded = ncum_r_2_regridded.sel(lat=lat_range, lon=lon_range)
ncum_r_3_regridded = ncum_r_3_regridded.sel(lat=lat_range, lon=lon_range)

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs = obs.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r_1_regridded = ncum_r_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r_2_regridded = ncum_r_2_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_r_3_regridded = ncum_r_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

# Extract data
obs_data = obs['rf'].mean(dim=['lat', 'lon']).values.ravel()
ncum_r1_data = ncum_r_1_regridded['APCP_surface'].mean(dim=['lat', 'lon']).values.ravel()
ncum_r2_data = ncum_r_2_regridded['APCP_surface'].mean(dim=['lat', 'lon']).values.ravel()
ncum_r3_data = ncum_r_3_regridded['APCP_surface'].mean(dim=['lat', 'lon']).values.ravel()

# Q-Q plot
def qq_plot(model_data, obs_data, label, ax):
    min_val = min(np.min(model_data), np.min(obs_data))
    max_val = max(np.max(model_data), np.max(obs_data))
    quantiles = np.linspace(0, 100, 100)
    
    model_quantiles = np.percentile(model_data, quantiles)
    obs_quantiles = np.percentile(obs_data, quantiles)

    ax.scatter(obs_quantiles, model_quantiles, color='blue', alpha=0.7, label=label)
    ax.plot([min_val, max_val], [min_val, max_val], '--', color='red', label="1:1 Line")
    
    ax.set_xlabel('Observed Quantiles')
    ax.set_ylabel('Model Quantiles')
    ax.set_title(f'Q-Q Plot: {label}')
    ax.legend()

# Plot Q-Q plots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
qq_plot(ncum_r1_data, obs_data, "NCUM R-Day1 vs Observation", axes[0])
qq_plot(ncum_r2_data, obs_data, "NCUM R-Day2 vs Observation", axes[1])
qq_plot(ncum_r3_data, obs_data, "NCUM R-Day3 vs Observation", axes[2])

plt.tight_layout()
plt.show()
