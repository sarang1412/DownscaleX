import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import mapping
import numpy as np
import rioxarray  

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


ncum_1 = ncum_1.rename({'latitude': 'lat', 'longitude': 'lon', 'APCP_surface': 'RAINFALL'})
ncum_3 = ncum_3.rename({'latitude': 'lat', 'longitude': 'lon', 'APCP_surface': 'RAINFALL'})
ncum_5 = ncum_5.rename({'latitude': 'lat', 'longitude': 'lon', 'APCP_surface': 'RAINFALL'})

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
    return ds.sel(lat=slice(21, 31), lon=slice(88, 103))

datasets = [obs, ncum_1, ncum_3, ncum_5, EQM_1, EQM_3, EQM_5, PQM_1, PQM_3, PQM_5, GPQM_1, GPQM_3, GPQM_5]
datasets = [slice_india(ds) for ds in datasets]

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

def compute_sedi_spatial(obs, model, threshold=64.5):
    hits = ((obs >= threshold) & (model >= threshold)).sum(dim='time')
    false_alarms = ((obs < threshold) & (model >= threshold)).sum(dim='time')
    misses = ((obs >= threshold) & (model < threshold)).sum(dim='time')
    correct_negatives = ((obs < threshold) & (model < threshold)).sum(dim='time')

    H = hits / (hits + misses)
    F = false_alarms / (false_alarms + correct_negatives)

    H = H.where((H > 0) & (H < 1))
    F = F.where((F > 0) & (F < 1))

    sedi = (np.log(F) - np.log(H) - np.log(1 - H) + np.log(1 - F)) / \
           (np.log(F) + np.log(H) + np.log(1 - H) + np.log(1 - F))
    return sedi

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
    sedi_scores[name] = compute_sedi_spatial(obs_data, model_data, threshold=64.5)

sedi_map = compute_sedi_spatial(obs["rf"], EQM_1["RAINFALL"])

plt.figure(figsize=(10, 6))

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl

mpl.rcParams.update({'font.size': 16}) 

# Plot
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(15, 16))
axes = axes.flatten()
cmap = plt.cm.Spectral
cmap.set_bad(color=(0, 0, 0, 0))  

for i, (name, dataset) in enumerate(datasets_dict.items()):
    ax = axes[i]
    sedi = compute_sedi_spatial(obs["rf"], dataset["RAINFALL"], threshold=64.5)
    sedi = sedi.where(sedi != 0) 

    sedi.plot(ax=ax, cmap=cmap, cbar_kwargs={'label': 'SEDI'})
    shape.boundary.plot(ax=ax, edgecolor='black', linewidth=0.8)

    ax.set_title(f"SEDI - {name.upper()}")
    ax.set_xlabel('')
    ax.set_ylabel('')

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()

plt.show()

