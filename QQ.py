import xarray as xr
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import geopandas as gpd
from shapely.geometry import mapping
import rioxarray

#data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc')
ncum_r_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day1rf_jjas2021-24.nc')
ncum_r_2 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day2rf_jjas2021-24.nc')
ncum_r_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumr_day3rf_jjas2021-24.nc')

# Trim 
obs = obs.sel(time=slice('2021-06-01', '2024-09-30'))

#lat/lon to match
ncum_r_1 = ncum_r_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_2 = ncum_r_2.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_r_3 = ncum_r_3.rename({'latitude': 'lat', 'longitude': 'lon'})

# Regrid 
ncum_r_1 = ncum_r_1.interp_like(obs, method='nearest')
ncum_r_2 = ncum_r_2.interp_like(obs, method='nearest')
ncum_r_3 = ncum_r_3.interp_like(obs, method='nearest')

# slicing
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_1 = ncum_r_1.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_2 = ncum_r_2.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_r_3 = ncum_r_3.sel(lat=slice(6, 41), lon=slice(65, 106))

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)

for ds in [obs, ncum_r_1, ncum_r_2, ncum_r_3]:
    ds = ds.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

# Extract rainfall
def ex_data(data_array, var_name):
    data = data_array[var_name].mean(dim=['time']).values.flatten()
    data = data[~np.isnan(data)]
    data = data[data > 0]
    return data

obs_data = ex_data(obs, 'rf')
ncum_r1_data = ex_data(ncum_r_1, 'APCP_surface')
ncum_r2_data = ex_data(ncum_r_2, 'APCP_surface')
ncum_r3_data = ex_data(ncum_r_3, 'APCP_surface')

#Q-Q plot
def qq(model_data, reference_data, label, ax):
    model_quantiles = np.sort(model_data)
    reference_quantiles = np.sort(reference_data)

    min_len = min(len(model_quantiles), len(reference_quantiles))
    model_quantiles = model_quantiles[:min_len]
    reference_quantiles = reference_quantiles[:min_len]

    ax.scatter(reference_quantiles, model_quantiles, label=label, alpha=0.6)
    ax.plot(reference_quantiles, reference_quantiles, 'r--', label='1:1 Line')
    ax.set_xlabel('Observed Rainfall Quantiles')
    ax.set_ylabel('Predicted Rainfall Quantiles')
    ax.set_title(f'Q-Q Plot: {label}')
    ax.legend()

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
qq(ncum_r1_data, obs_data, "NCUM R Day1 vs OBS", axes[0])
qq(ncum_r2_data, obs_data, "NCUM R Day2 vs OBS", axes[1])
qq(ncum_r3_data, obs_data, "NCUM R Day3 vs OBS", axes[2])

plt.tight_layout()
plt.show()
