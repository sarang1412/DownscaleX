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

def binary_classification(obs_data, model_data, threshold=64.5):
    obs_event = obs_data > threshold
    model_event = model_data > threshold

    hit = ((model_event) & (obs_event)).sum(dim=["time", "lat", "lon"])
    miss = ((~model_event) & (obs_event)).sum(dim=["time", "lat", "lon"])
    false_alarm = ((model_event) & (~obs_event)).sum(dim=["time", "lat", "lon"])
    correct_negative = ((~model_event) & (~obs_event)).sum(dim=["time", "lat", "lon"])
    
    return {
        "Hit": hit.item(),
        "Miss": miss.item(),
        "False Alarm": false_alarm.item(),
        "Correct Negative": correct_negative.item()
    }

def compute_sedi(stats_dict, eps=1e-10):
    hit = stats_dict["Hit"]
    miss = stats_dict["Miss"]
    fa = stats_dict["False Alarm"]
    cn = stats_dict["Correct Negative"]

    H = hit / (hit + miss + eps)
    F = fa / (fa + cn + eps)

    # Bound values away from 0 and 1
    H = np.clip(H, eps, 1 - eps)
    F = np.clip(F, eps, 1 - eps)

    sedi_numerator = np.log(F) - np.log(H) - np.log(1 - F) + np.log(1 - H)
    sedi_denominator = np.log(F) + np.log(H) + np.log(1 - F) + np.log(1 - H)
    sedi = sedi_numerator / sedi_denominator

    return sedi

# === Compute SEDI scores for all datasets ===
datasets_dict = {
    "ncum_1": ncum_1, "ncum_3": ncum_3, "ncum_5": ncum_5,
    "EQM_1": EQM_1, "EQM_3": EQM_3, "EQM_5": EQM_5,
    "PQM_1": PQM_1, "PQM_3": PQM_3, "PQM_5": PQM_5,
    "GPQM_1": GPQM_1, "GPQM_3": GPQM_3, "GPQM_5": GPQM_5
}

sedi_scores = {}

for name, dataset in datasets_dict.items():
    model_data = dataset["RAINFALL"]
    obs_data = obs["rf"]
    stats = binary_classification(obs_data, model_data, threshold=64.5)
    sedi_scores[name] = compute_sedi(stats)

# === Plot SEDI scores as a bar plot ===
plt.figure(figsize=(12, 6))
plt.bar(sedi_scores.keys(), sedi_scores.values(), color='mediumseagreen', edgecolor='black')
plt.xlabel('Datasets', fontsize=12)
plt.ylabel('SEDI Score', fontsize=12)
plt.title('SEDI Scores for Different Datasets', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.yticks(np.arange(0, 0.7, 0.1))  # Add ticks at intervals of 0.1
plt.grid(axis='y', linestyle='--', alpha=0.7)  # Add horizontal grid lines
plt.tight_layout()
plt.show()