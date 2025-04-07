import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import mapping
import numpy as np
import rioxarray

obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc')

ncum_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc')
ncum_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day3rf_jjas2020-24.nc')
ncum_5 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day5rf_jjas2020-24.nc')

EQM_1 = xr.open_dataset('e:\\Dissertation\\data\\EQM_ncumg_day1rf_jjas2020-24-0p25.nc')
EQM_3 = xr.open_dataset('e:\\Dissertation\\data\\EQM_ncumg_day3rf_jjas2020-24-0p25.nc')
EQM_5 = xr.open_dataset('e:\\Dissertation\\data\\EQM_ncumg_day5rf_jjas2020-24-0p25.nc')

PQM_1 = xr.open_dataset('e:\\Dissertation\\data\\PQM_ncumg_day1rf_jjas2020-24-0p25.nc')
PQM_3 = xr.open_dataset('e:\\Dissertation\\data\\PQM_ncumg_day3rf_jjas2020-24-0p25.nc')
PQM_5 = xr.open_dataset('e:\\Dissertation\\data\\PQM_ncumg_day5rf_jjas2020-24-0p25.nc')

GPQM_1 = xr.open_dataset('e:\\Dissertation\\data\\GPQM_ncumg_day1rf_jjas2020-24-0p25.nc')
GPQM_3 = xr.open_dataset('e:\\Dissertation\\data\\GPQM_ncumg_day3rf_jjas2020-24-0p25.nc')
GPQM_5 = xr.open_dataset('e:\\Dissertation\\data\\GPQM_ncumg_day5rf_jjas2020-24-0p25.nc')

ncum_1 = ncum_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_3 = ncum_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_5 = ncum_5.rename({'latitude': 'lat', 'longitude': 'lon'})

regrid_to_obs = lambda ds: ds.interp_like(obs, method='nearest')

ncum_1 = regrid_to_obs(ncum_1)
ncum_3 = regrid_to_obs(ncum_3)
ncum_5 = regrid_to_obs(ncum_5)

EQM_1 = regrid_to_obs(EQM_1)
EQM_3 = regrid_to_obs(EQM_3)
EQM_5 = regrid_to_obs(EQM_5)

PQM_1 = regrid_to_obs(PQM_1)
PQM_3 = regrid_to_obs(PQM_3)
PQM_5 = regrid_to_obs(PQM_5)

GPQM_1 = regrid_to_obs(GPQM_1)
GPQM_3 = regrid_to_obs(GPQM_3)
GPQM_5 = regrid_to_obs(GPQM_5)

def slice_india(ds):
    return ds.sel(lat=slice(6, 41), lon=slice(65, 106))

datasets = [obs, ncum_1, ncum_3, ncum_5, EQM_1, EQM_3, EQM_5, PQM_1, PQM_3, PQM_5, GPQM_1, GPQM_3, GPQM_5]
datasets = [slice_india(ds) for ds in datasets]

shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)

def clip_to_shape(ds):
    return ds.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

obs, ncum_1, ncum_3, ncum_5, EQM_1, EQM_3, EQM_5, PQM_1, PQM_3, PQM_5, GPQM_1, GPQM_3, GPQM_5 = [clip_to_shape(ds) for ds in datasets]

def clean_flatten(ds, varname):
    arr = ds[varname].mean(dim='time').values.flatten()
    arr = arr[~np.isnan(arr)]
    return arr[arr > 0]

obs_data     = clean_flatten(obs, "rf")
ncum1_data   = clean_flatten(ncum_1, "APCP_surface")
ncum3_data   = clean_flatten(ncum_3, "APCP_surface")
ncum5_data   = clean_flatten(ncum_5, "APCP_surface")
EQM1_data    = clean_flatten(EQM_1, "RAINFALL")
EQM3_data    = clean_flatten(EQM_3, "RAINFALL")
EQM5_data    = clean_flatten(EQM_5, "RAINFALL")
PQM1_data    = clean_flatten(PQM_1, "RAINFALL")
PQM3_data    = clean_flatten(PQM_3, "RAINFALL")
PQM5_data    = clean_flatten(PQM_5, "RAINFALL")
GPQM1_data   = clean_flatten(GPQM_1, "RAINFALL")
GPQM3_data   = clean_flatten(GPQM_3, "RAINFALL")
GPQM5_data   = clean_flatten(GPQM_5, "RAINFALL")

def plot_qq(obs_data, model_data, label, ax):
    # Sort both
    obs_sorted = np.sort(obs_data)
    model_sorted = np.sort(model_data)
    
    # Ensure equal lengths
    min_len = min(len(obs_sorted), len(model_sorted))
    obs_sorted = obs_sorted[:min_len]
    model_sorted = model_sorted[:min_len]

    ax.scatter(obs_sorted, model_sorted, alpha=0.5, s=10, label=label)
    ax.plot([obs_sorted.min(), obs_sorted.max()],
            [obs_sorted.min(), obs_sorted.max()],
            color='red', linestyle='--', lw=1)
    ax.set_xlabel("Observed Rainfall (mm)")
    ax.set_ylabel("Model Rainfall (mm)")
    ax.set_title(f"QQ Plot: {label}")
    ax.grid(True)
    ax.legend()  # Ensure legend is added to the correct axis

models1 = {
    "EQM Day1": EQM1_data,
    "PQM Day1": PQM1_data,
    "GPQM Day1": GPQM1_data,
}

models3 = {
    "EQM Day3": EQM3_data,
    "PQM Day3": PQM3_data,
    "GPQM Day3": GPQM3_data,
}

models5 = {
    "EQM Day5": EQM5_data,
    "PQM Day5": PQM5_data,
    "GPQM Day5": GPQM5_data
}

fig, axs = plt.subplots(3, 3, figsize=(15, 18))
axs = axs.flatten()

for i, (label, model_data) in enumerate(models1.items()):
    plot_qq(obs_data, model_data, label, axs[i])
    plot_qq(obs_data, ncum1_data, "NCUM Day 1", axs[i])
for i, (label, model_data) in enumerate(models3.items()):
    plot_qq(obs_data, model_data, label, axs[3+i])
    plot_qq(obs_data, ncum3_data, "NCUM Day 3", axs[3+i])
for i, (label, model_data) in enumerate(models5.items()):
    plot_qq(obs_data, model_data, label, axs[6+i])
    plot_qq(obs_data, ncum5_data, "NCUM Day 5", axs[6+i])
    
plt.tight_layout()
plt.suptitle("QQ Plots of Models vs Observed Rainfall", fontsize=16, y=1.02)
plt.show()
