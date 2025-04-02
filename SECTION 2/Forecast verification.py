import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import mapping
import numpy as np
import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import mapping

#load data
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc')
ncum_1 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc')
ncum_3 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day3rf_jjas2020-24.nc')
ncum_5 = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day5rf_jjas2020-24.nc')
EQM_ncum_1 = xr.open_dataset('e:\\Dissertation\\data\\EQM_ncumg_day1rf_jjas2020-24-0p25.nc')
EQM_ncum_3 = xr.open_dataset('e:\\Dissertation\\data\\EQM_ncumg_day3rf_jjas2020-24-0p25.nc')
EQM_ncum_5 = xr.open_dataset('e:\\Dissertation\\data\\EQM_ncumg_day5rf_jjas2020-24-0p25.nc')
PQM_ncum_1 = xr.open_dataset('e:\\Dissertation\\data\\PQM_ncumg_day1rf_jjas2020-24-0p25.nc')
PQM_ncum_3 = xr.open_dataset('e:\\Dissertation\\data\\PQM_ncumg_day3rf_jjas2020-24-0p25.nc')
PQM_ncum_5 = xr.open_dataset('e:\\Dissertation\\data\\PQM_ncumg_day5rf_jjas2020-24-0p25.nc')
GPQM_ncum_1 = xr.open_dataset('e:\\Dissertation\\data\\GPQM_ncumg_day1rf_jjas2020-24-0p25.nc')
GPQM_ncum_3 = xr.open_dataset('e:\\Dissertation\\data\\GPQM_ncumg_day3rf_jjas2020-24-0p25.nc')
GPQM_ncum_5 = xr.open_dataset('e:\\Dissertation\\data\\GPQM_ncumg_day5rf_jjas2020-24-0p25.nc')

# Trim dataset
#obs = obs.sel(time=slice('2021-06-01T03:00:00.000000000', '2024-09-30T03:00:00.000000000'))

#change lat and lon to match
ncum_1 = ncum_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_3 = ncum_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_5 = ncum_5.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
obs = obs.interp_like(ncum_1, method='nearest')
EQM_ncum_1_regridded = EQM_ncum_1.interp_like(obs, method='nearest')
EQM_ncum_3_regridded = EQM_ncum_3.interp_like(obs, method='nearest')
EQM_ncum_5_regridded = EQM_ncum_5.interp_like(obs, method='nearest')
PQM_ncum_1_regridded = PQM_ncum_1.interp_like(obs, method='nearest')
PQM_ncum_3_regridded = PQM_ncum_3.interp_like(obs, method='nearest')
PQM_ncum_5_regridded = PQM_ncum_5.interp_like(obs, method='nearest')
GPQM_ncum_1_regridded = GPQM_ncum_1.interp_like(obs, method='nearest')
GPQM_ncum_3_regridded = GPQM_ncum_3.interp_like(obs, method='nearest')
GPQM_ncum_5_regridded = GPQM_ncum_5.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(65, 106))#62
ncum_1 = ncum_1.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_5 = ncum_5.sel(lat=slice(6, 41), lon=slice(65, 106))
ncum_3 = ncum_3.sel(lat=slice(6, 41), lon=slice(65, 106))
EQM_ncum_1_regridded = EQM_ncum_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
EQM_ncum_5_regridded = EQM_ncum_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
EQM_ncum_3_regridded = EQM_ncum_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
PQM_ncum_1_regridded = PQM_ncum_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
PQM_ncum_5_regridded = PQM_ncum_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
PQM_ncum_3_regridded = PQM_ncum_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
GPQM_ncum_1_regridded = GPQM_ncum_1_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
GPQM_ncum_5_regridded = GPQM_ncum_5_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))
GPQM_ncum_3_regridded = GPQM_ncum_3_regridded.sel(lat=slice(6, 41), lon=slice(65, 106))

# Shapefile
shapefile_path = "e:\\Dissertation\\data\\SHP&DEM\\shp\\India_State_Boundary_02may2020.shp"
shape = gpd.read_file(shapefile_path)
#shape = shape.to_crs(epsg=4326)

obs_data = obs.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_1_data = ncum_1.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_5_data = ncum_5.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
ncum_3_data = ncum_3.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
EQM_ncum_1_data = EQM_ncum_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
EQM_ncum_5_data = EQM_ncum_5_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
EQM_ncum_3_data = EQM_ncum_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
PQM_ncum_1_data = PQM_ncum_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
PQM_ncum_5_data = PQM_ncum_5_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
PQM_ncum_3_data = PQM_ncum_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
GPQM_ncum_1_data = GPQM_ncum_1_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
GPQM_ncum_5_data = GPQM_ncum_5_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)
GPQM_ncum_3_data = GPQM_ncum_3_regridded.rio.write_crs("EPSG:4326").rio.clip(shape.geometry.apply(mapping), shape.crs)

obs_datas = obs_data["rf"]
ncum_1_data = ncum_1_data["APCP_surface"]
ncum_5_data = ncum_5_data["APCP_surface"]
ncum_3_data = ncum_3_data["APCP_surface"]
EQM_ncum_1_data = EQM_ncum_1_data["RAINFALL"]
EQM_ncum_5_data = EQM_ncum_5_data["RAINFALL"]
EQM_ncum_3_data = EQM_ncum_3_data["RAINFALL"]
PQM_ncum_1_data = PQM_ncum_1_data["RAINFALL"]
PQM_ncum_5_data = PQM_ncum_5_data["RAINFALL"]
PQM_ncum_3_data = PQM_ncum_3_data["RAINFALL"]
GPQM_ncum_1_data = GPQM_ncum_1_data["RAINFALL"]
GPQM_ncum_5_data = GPQM_ncum_5_data["RAINFALL"]
GPQM_ncum_3_data = GPQM_ncum_3_data["RAINFALL"]

# Filter data for threshold 0 - 30 mm
obss = obs_datas.where((obs >= 0) & (obs <= 30), drop=True)
ncum_1_data = ncum_1_data.where((ncum_1_data >= 0) & (ncum_1_data <= 30), drop=True)
ncum_3_data = ncum_3_data.where((ncum_3_data >= 0) & (ncum_3_data <= 30), drop=True)
ncum_5_data = ncum_5_data.where((ncum_5_data >= 0) & (ncum_5_data <= 30), drop=True)
EQM_ncum_1_data = EQM_ncum_1_data.where((EQM_ncum_1_data >= 0) & (EQM_ncum_1_data <= 30), drop=True)
EQM_ncum_3_data = EQM_ncum_3_data.where((EQM_ncum_3_data >= 0) & (EQM_ncum_3_data <= 30), drop=True)
EQM_ncum_5_data = EQM_ncum_5_data.where((EQM_ncum_5_data >= 0) & (EQM_ncum_5_data <= 30), drop=True)
PQM_ncum_1_data = PQM_ncum_1_data.where((PQM_ncum_1_data >= 0) & (PQM_ncum_1_data <= 30), drop=True)
PQM_ncum_3_data = PQM_ncum_3_data.where((PQM_ncum_3_data >= 0) & (PQM_ncum_3_data <= 30), drop=True)
PQM_ncum_5_data = PQM_ncum_5_data.where((PQM_ncum_5_data >= 0) & (PQM_ncum_5_data <= 30), drop=True)
GPQM_ncum_1_data = GPQM_ncum_1_data.where((GPQM_ncum_1_data >= 0) & (GPQM_ncum_1_data <= 30), drop=True)
GPQM_ncum_3_data = GPQM_ncum_3_data.where((GPQM_ncum_3_data >= 0) & (GPQM_ncum_3_data <= 30), drop=True)
GPQM_ncum_5_data = GPQM_ncum_5_data.where((GPQM_ncum_5_data >= 0) & (GPQM_ncum_5_data <= 30), drop=True)

# rainfall to binary
threshold = 1
obs_binary = (obss > threshold).astype(int)
obs_binary = obs_binary['rf']
ncum_1_binary = (ncum_1_data > threshold).astype(int)
ncum_3_binary = (ncum_3_data > threshold).astype(int)
ncum_5_binary = (ncum_5_data > threshold).astype(int)
EQM_ncum_1_binary = (EQM_ncum_1_data > threshold).astype(int)
EQM_ncum_3_binary = (EQM_ncum_3_data > threshold).astype(int)
EQM_ncum_5_binary = (EQM_ncum_5_data > threshold).astype(int)
PQM_ncum_1_binary = (PQM_ncum_1_data > threshold).astype(int)
PQM_ncum_3_binary = (PQM_ncum_3_data > threshold).astype(int)
PQM_ncum_5_binary = (PQM_ncum_5_data > threshold).astype(int)
GPQM_ncum_1_binary = (GPQM_ncum_1_data > threshold).astype(int)
GPQM_ncum_3_binary = (GPQM_ncum_3_data > threshold).astype(int)
GPQM_ncum_5_binary = (GPQM_ncum_5_data > threshold).astype(int)

def hitORmiss(obs, fcst):
    H = ((obs == 1) & (fcst == 1)).sum().item()
    M = ((obs == 1) & (fcst == 0)).sum().item()
    F = ((obs == 0) & (fcst == 1)).sum().item()
    C = ((obs == 0) & (fcst == 0)).sum().item()
    return H, M, F, C

H1, M1, F1, C1 = hitORmiss(obs_binary, ncum_1_binary)
H2, M2, F2, C2 = hitORmiss(obs_binary, ncum_3_binary)
H3, M3, F3, C3 = hitORmiss(obs_binary, ncum_5_binary)
H4, M4, F4, C4 = hitORmiss(obs_binary, EQM_ncum_1_binary)
H5, M5, F5, C5 = hitORmiss(obs_binary, EQM_ncum_3_binary)
H6, M6, F6, C6 = hitORmiss(obs_binary, EQM_ncum_5_binary)
H7, M7, F7, C7 = hitORmiss(obs_binary, PQM_ncum_1_binary)
H8, M8, F8, C8 = hitORmiss(obs_binary, PQM_ncum_3_binary)
H9, M9, F9, C9 = hitORmiss(obs_binary, PQM_ncum_5_binary)
H10, M10, F10, C10 = hitORmiss(obs_binary, GPQM_ncum_1_binary)
H11, M11, F11, C11 = hitORmiss(obs_binary, GPQM_ncum_3_binary)
H12, M12, F12, C12 = hitORmiss(obs_binary, GPQM_ncum_5_binary)


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
POD4, FAR4, Bias4, ETS4 = compute_scores(H4, M4, F4, C4)
POD5, FAR5, Bias5, ETS5 = compute_scores(H5, M5, F5, C5)
POD6, FAR6, Bias6, ETS6 = compute_scores(H6, M6, F6, C6)
POD7, FAR7, Bias7, ETS7 = compute_scores(H7, M7, F7, C7)
POD8, FAR8, Bias8, ETS8 = compute_scores(H8, M8, F8, C8)
POD9, FAR9, Bias9, ETS9 = compute_scores(H9, M9, F9, C9)
POD10, FAR10, Bias10, ETS10 = compute_scores(H10, M10, F10, C10)
POD11, FAR11, Bias11, ETS11 = compute_scores(H11, M11, F11, C11)
POD12, FAR12, Bias12, ETS12 = compute_scores(H12, M12, F12, C12)

metrics = ['POD', 'FAR', 'Bias', 'ETS']
values_1 = [POD1, FAR1, Bias1, ETS1]
values_2 = [POD2, FAR2, Bias2, ETS2]
values_3 = [POD3, FAR3, Bias3, ETS3]
values_4 = [POD4, FAR4, Bias4, ETS4]
values_5 = [POD5, FAR5, Bias5, ETS5]
values_6 = [POD6, FAR6, Bias6, ETS6]
values_7 = [POD7, FAR7, Bias7, ETS7]
values_8 = [POD8, FAR8, Bias8, ETS8]
values_9 = [POD9, FAR9, Bias9, ETS9]
values_10 = [POD10, FAR10, Bias10, ETS10]
values_11 = [POD11, FAR11, Bias11, ETS11]
values_12 = [POD12, FAR12, Bias12, ETS12]

x = range(len(metrics))
lead_times = [1, 2, 3]

POD_values = [POD1, POD2, POD3]
POD_values_eqm = [POD4, POD5, POD6]
POD_values_pqm = [POD7, POD8, POD9]
POD_values_gpqm = [POD10, POD11, POD12]
FAR_values = [FAR1, FAR2, FAR3]
FAR_values_eqm = [FAR4, FAR5, FAR6]
FAR_values_pqm = [FAR7, FAR8, FAR9]
FAR_values_gpqm = [FAR10, FAR11, FAR12]
Bias_values = [Bias1, Bias2, Bias3]
Bias_values_eqm = [Bias4, Bias5, Bias6]
Bias_values_pqm = [Bias7, Bias8, Bias9]
Bias_values_gpqm = [Bias10, Bias11, Bias12]
ETS_values = [ETS1, ETS2, ETS3]
ETS_values_eqm = [ETS4, ETS5, ETS6]
ETS_values_pqm = [ETS7, ETS8, ETS9]
ETS_values_gpqm = [ETS10, ETS11, ETS12]

# Combine all plots into a single figure with 2x2 subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))  # Create a 2x2 grid of subplots
axes = axes.flatten()  # Flatten the 2D array of axes for easier indexing

# Create the plots for each dataset
datasets = [
    ("NCUM_G", POD_values, FAR_values, Bias_values, ETS_values),
    ("EQM", POD_values_eqm, FAR_values_eqm, Bias_values_eqm, ETS_values_eqm),
    ("PQM", POD_values_pqm, FAR_values_pqm, Bias_values_pqm, ETS_values_pqm),
    ("GPQM", POD_values_gpqm, FAR_values_gpqm, Bias_values_gpqm, ETS_values_gpqm),
]

for i, (dataset_name, POD, FAR, Bias, ETS) in enumerate(datasets):
    ax = axes[i]  # Select the corresponding subplot
    
    # Plot each metric with different styles
    ax.plot(lead_times, POD, 'bo-', label="POD")  # Blue circles
    ax.plot(lead_times, FAR, 'rs-', label="FAR")  # Red squares
    ax.plot(lead_times, Bias, 'g^-', label="BIAS")  # Green triangles
    ax.plot(lead_times, ETS, 'mp-', label="ETS")  # Purple pentagons
    
    # Formatting
    ax.set_xlabel("Forecast Lead Time (Days)")
    ax.set_ylabel("Score")
    ax.set_title(f"{dataset_name} Metrics")
    ax.legend()
    ax.grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the combined plot
plt.show()

