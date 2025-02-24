import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import geopandas as gpd
import rioxarray
from shapely.geometry import mapping


#load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc' )
ncum_g_1 = xr.open_dataset('e:\\Dissertation\\data\\ncum_day1rf_jjas2021-24.nc')
ncum_g_2 = xr.open_dataset('e:\\Dissertation\\data\\ncum_day2rf_jjas2021-24.nc')
ncum_g_3 = xr.open_dataset('e:\\Dissertation\\data\\ncum_day3rf_jjas2021-24.nc')

# Trim dataset
obs = obs.sel(time=slice('2021-06-01T03:00:00.000000000', '2024-09-30T03:00:00.000000000'))

#change lat and lon to match
ncum_g_1 = ncum_g_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_3 = ncum_g_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_g_2 = ncum_g_2.rename({'latitude': 'lat', 'longitude': 'lon'})
#regrid
ncum_g_1_regridded = ncum_g_1.interp_like(obs, method='nearest')
ncum_g_3_regridded = ncum_g_3.interp_like(obs, method='nearest')
ncum_g_2_regridded = ncum_g_2.interp_like(obs, method='nearest')
#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_1_regridded = ncum_g_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_3_regridded = ncum_g_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_g_2_regridded = ncum_g_2_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

# mean rainfall
obs =obs['rf']
ncum_g_1 = ncum_g_1_regridded['APCP_surface']
ncum_g_3 = ncum_g_3_regridded['APCP_surface']
ncum_g_2 = ncum_g_2_regridded['APCP_surface']

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs = obs.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_1 = ncum_g_1.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_3 = ncum_g_3.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_g_2 = ncum_g_2.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

# rainfall to binary
threshold = 1
obs_binary = (obs > threshold).astype(int)
ncum_g_1_binary = (ncum_g_1 > threshold).astype(int)
ncum_g_2_binary = (ncum_g_2 > threshold).astype(int)
ncum_g_3_binary = (ncum_g_3 > threshold).astype(int)

def hitORmiss(obs, fcst):
    H = ((obs == 1) & (fcst == 1)).sum().item()
    M = ((obs == 1) & (fcst == 0)).sum().item()
    F = ((obs == 0) & (fcst == 1)).sum().item()
    C = ((obs == 0) & (fcst == 0)).sum().item()
    return H, M, F, C

H1, M1, F1, C1 = hitORmiss(obs_binary, ncum_g_1_binary)
H2, M2, F2, C2 = hitORmiss(obs_binary, ncum_g_2_binary)
H3, M3, F3, C3 = hitORmiss(obs_binary, ncum_g_3_binary)

def compute_scores(H, M, F, C):
    POD = H / (H + M) if (H + M) > 0 else 0
    FAR = F / (H + F) if (H + F) > 0 else 0
    Bias = (H + F) / (H + M) if (H + M) > 0 else 0
    H_random = (H + M) * (H + F) / (H + M + F + C)
    ETS = (H - H_random) / (H + M + F - H_random) if (H + M + F - H_random) > 0 else 0
    return POD, FAR, Bias, ETS

POD1, FAR1, Bias1, ETS1 = compute_scores(H1, M1, F1, C1)
POD2, FAR2, Bias2, ETS2 = compute_scores(H2, M2, F2, C2)
POD3, FAR3, Bias3, ETS3 = compute_scores(H3, M3, F3, C3)

#plot
metrics = ['POD', 'FAR', 'Bias', 'ETS']
values_1 = [POD1, FAR1, Bias1, ETS1]
values_2 = [POD2, FAR2, Bias2, ETS2]
values_3 = [POD3, FAR3, Bias3, ETS3]

x = range(len(metrics))

plt.figure(figsize=(8,5))
plt.bar(x, values_1, width=0.2, label="Day 1", color='red')
plt.bar([i+0.2 for i in x], values_2, width=0.2, label="Day 2", color='blue')
plt.bar([i+0.4 for i in x], values_3, width=0.2, label="Day 3", color='green')

plt.xticks([i+0.2 for i in x], metrics)
plt.ylabel("Score")
plt.title("Forecast Verification Metrics (NCUM-G)")
plt.legend()
plt.grid(axis='y', linestyle='--', alpha=0.7)

for i, v in enumerate(values_1):
    plt.text(i, v + 0.02, f"{v:.2f}", ha='center', fontsize=10, color='red')
for i, v in enumerate(values_2):
    plt.text(i+0.2, v + 0.02, f"{v:.2f}", ha='center', fontsize=10, color='blue')
for i, v in enumerate(values_3):
    plt.text(i+0.4, v + 0.02, f"{v:.2f}", ha='center', fontsize=10, color='green')

plt.show()


