import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import mapping
import numpy as np
import rioxarray  # Make sure this is installed

# === Load datasets ===
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

# === Harmonise dimension names ===
ncum_1 = ncum_1.rename({'latitude': 'lat', 'longitude': 'lon', 'APCP_surface': 'RAINFALL'})
ncum_3 = ncum_3.rename({'latitude': 'lat', 'longitude': 'lon', 'APCP_surface': 'RAINFALL'})
ncum_5 = ncum_5.rename({'latitude': 'lat', 'longitude': 'lon', 'APCP_surface': 'RAINFALL'})
# === Regrid model and corrected data to match obs ===
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

# === Slice all datasets to India domain ===
def slice_india(ds):
    return ds.sel(lat=slice(6, 41), lon=slice(65, 106))

datasets = [obs, ncum_1, ncum_3, ncum_5, EQM_1, EQM_3, EQM_5, PQM_1, PQM_3, PQM_5, GPQM_1, GPQM_3, GPQM_5]
datasets = [slice_india(ds) for ds in datasets]

# === Clip using shapefile ===
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)

def clip_to_shape(ds):
    return ds.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

obs, ncum_1, ncum_3, ncum_5, EQM_1, EQM_3, EQM_5, PQM_1, PQM_3, PQM_5, GPQM_1, GPQM_3, GPQM_5 = [clip_to_shape(ds) for ds in datasets]


def calculate_pod(obs_data, model_data, thresholds):
    pod_values = []
    for threshold in thresholds:
        obs_event = obs_data > threshold
        model_event = model_data > threshold
        hit = ((model_event) & (obs_event)).sum(dim=["time", "lat", "lon"])
        miss = ((~model_event) & (obs_event)).sum(dim=["time", "lat", "lon"])
        pod = hit / (hit + miss)
        pod_values.append(pod.item())
    return pod_values

def calculate_far(obs_data, model_data, thresholds):
    far_values = []
    for threshold in thresholds:
        obs_event = obs_data > threshold
        model_event = model_data > threshold
        false_alarm = ((model_event) & (~obs_event)).sum(dim=["time", "lat", "lon"])
        hit = ((model_event) & (obs_event)).sum(dim=["time", "lat", "lon"])
        far = false_alarm / (false_alarm + hit)
        far_values.append(far.item())
    return far_values

def calculate_ets(obs_data, model_data, thresholds):
    ets_values = []
    for threshold in thresholds:
        obs_event = obs_data > threshold
        model_event = model_data > threshold
        hit = ((model_event) & (obs_event)).sum(dim=["time", "lat", "lon"])
        false_alarm = ((model_event) & (~obs_event)).sum(dim=["time", "lat", "lon"])
        miss = ((~model_event) & (obs_event)).sum(dim=["time", "lat", "lon"])
        random_hits = ((hit + miss) * (hit + false_alarm)) / obs_event.size
        ets = (hit - random_hits) / (hit + miss + false_alarm - random_hits)
        ets_values.append(ets.item())
    return ets_values


thresholds = np.arange(0, 301, 10)  # Thresholds from 0 to 300 with step of 10

datasets_dict = {
    "NCUM": [ncum_1, ncum_3, ncum_5],
    "EQM": [EQM_1, EQM_3, EQM_5],
    "PQM": [PQM_1, PQM_3, PQM_5],
    "GPQM": [GPQM_1, GPQM_3, GPQM_5]
}

metrics = {
    "POD": calculate_pod,
    "FAR": calculate_far,
    "ETS": calculate_ets
}

for day, day_index in zip([1, 3, 5], [0, 1, 2]):  # Day 1, Day 3, Day 5
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))  # 1x3 grid for 3 metrics
    metric_names = list(metrics.keys())
    
    for idx, (metric_name, metric_func) in enumerate(metrics.items()):
        ax = axes[idx]
        for label, datasets in datasets_dict.items():
            dataset = datasets[day_index]
            metric_values = metric_func(obs["rf"], dataset["RAINFALL"], thresholds)
            ax.plot(thresholds, metric_values, label=f"{label}")
        ax.set_title(f"{metric_name} for Day {day}")
        ax.set_xlabel("Threshold (mm)")
        ax.set_ylabel(f"{metric_name} Value")
        ax.legend()
        ax.grid()
    
    plt.suptitle(f"Metrics vs Threshold for Day {day}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit the suptitle
    plt.show()
