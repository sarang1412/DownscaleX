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

var_map = {
    "obs": "rf",              
    "ncum": "APCP_surface",    
    "bias": "RAINFALL"         
}

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

def compute_iqr(data):
    q75, q25 = np.percentile(data, [75, 25])
    return q75 - q25

def extract_monthly_iqr(ds, varname, month):
    monthly = ds.sel(time=ds.time.dt.month == month)
    data = clean_flatten(monthly, varname)
    return compute_iqr(data)

# Store IQRs
months = [6, 7, 8, 9]
month_names = ["June", "July", "August", "September"]
methods = ["OBS", "Raw", "EQM", "PQM", "GPQM"]

iqr_d1 = {m: [] for m in methods}
iqr_d3 = {m: [] for m in methods}
iqr_d5 = {m: [] for m in methods}

for m in months:
    iqr_d1["OBS"].append(extract_monthly_iqr(obs, var_map["obs"], m))
    iqr_d1["Raw"].append(extract_monthly_iqr(ncum_1, var_map["ncum"], m))
    iqr_d1["EQM"].append(extract_monthly_iqr(EQM_1, var_map["bias"], m))
    iqr_d1["PQM"].append(extract_monthly_iqr(PQM_1, var_map["bias"], m))
    iqr_d1["GPQM"].append(extract_monthly_iqr(GPQM_1, var_map["bias"], m))

    iqr_d3["OBS"].append(extract_monthly_iqr(obs, var_map["obs"], m))
    iqr_d3["Raw"].append(extract_monthly_iqr(ncum_3, var_map["ncum"], m))
    iqr_d3["EQM"].append(extract_monthly_iqr(EQM_3, var_map["bias"], m))
    iqr_d3["PQM"].append(extract_monthly_iqr(PQM_3, var_map["bias"], m))
    iqr_d3["GPQM"].append(extract_monthly_iqr(GPQM_3, var_map["bias"], m))

    iqr_d5["OBS"].append(extract_monthly_iqr(obs, var_map["obs"], m))
    iqr_d5["Raw"].append(extract_monthly_iqr(ncum_5, var_map["ncum"], m))
    iqr_d5["EQM"].append(extract_monthly_iqr(EQM_5, var_map["bias"], m))
    iqr_d5["PQM"].append(extract_monthly_iqr(PQM_5, var_map["bias"], m))
    iqr_d5["GPQM"].append(extract_monthly_iqr(GPQM_5, var_map["bias"], m))

fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

forecast_days = ["Day 1", "Day 3", "Day 5"]
iqr_all = [iqr_d1, iqr_d3, iqr_d5]
colors = ["black", "red", "green", "orange", "purple"]

for ax, day, iqr_dict in zip(axes, forecast_days, iqr_all):
    for method, color in zip(methods, colors):
        ax.plot(month_names, iqr_dict[method], marker='o', label=method, color=color)
    ax.set_ylabel("IQR (mm)")
    ax.set_title(f"IQR of Rainfall – {day}")
    ax.grid(True)

axes[-1].set_xlabel("Month")
axes[0].legend(ncol=5, loc='upper center', bbox_to_anchor=(0.5, 1.35))
plt.tight_layout()
plt.show()
