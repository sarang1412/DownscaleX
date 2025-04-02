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
EQM_ncum_1_regridded = EQM_ncum_1.interp_like(ncum_1, method='nearest')
EQM_ncum_3_regridded = EQM_ncum_3.interp_like(ncum_3, method='nearest')
EQM_ncum_5_regridded = EQM_ncum_5.interp_like(ncum_5, method='nearest')
PQM_ncum_1_regridded = PQM_ncum_1.interp_like(ncum_1, method='nearest')
PQM_ncum_3_regridded = PQM_ncum_3.interp_like(ncum_3, method='nearest')
PQM_ncum_5_regridded = PQM_ncum_5.interp_like(ncum_5, method='nearest')
GPQM_ncum_1_regridded = GPQM_ncum_1.interp_like(ncum_1, method='nearest')
GPQM_ncum_3_regridded = GPQM_ncum_3.interp_like(ncum_3, method='nearest')
GPQM_ncum_5_regridded = GPQM_ncum_5.interp_like(ncum_5, method='nearest')

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

'''
# Extract data
obs_data = obs_data['rf'].mean(dim=['lat', 'lon']).values.ravel()
ncum_1_data = ncum_1_data['APCP_surface'].mean(dim=['lat', 'lon']).values.ravel()
ncum_5_data = ncum_5_data['APCP_surface'].mean(dim=['lat', 'lon']).values.ravel()
ncum_3_data = ncum_3_data['APCP_surface'].mean(dim=['lat', 'lon']).values.ravel()
EQM_ncum_1_data = EQM_ncum_1_data['RAINFALL'].mean(dim=['lat', 'lon']).values.ravel()
EQM_ncum_5_data = EQM_ncum_5_data['RAINFALL'].mean(dim=['lat', 'lon']).values.ravel()
EQM_ncum_3_data = EQM_ncum_3_data['RAINFALL'].mean(dim=['lat', 'lon']).values.ravel()
PQM_ncum_1_data = PQM_ncum_1_data['RAINFALL'].mean(dim=['lat', 'lon']).values.ravel()
PQM_ncum_5_data = PQM_ncum_5_data['RAINFALL'].mean(dim=['lat', 'lon']).values.ravel()
PQM_ncum_3_data = PQM_ncum_3_data['RAINFALL'].mean(dim=['lat', 'lon']).values.ravel()
GPQM_ncum_1_data = GPQM_ncum_1_data['RAINFALL'].mean(dim=['lat', 'lon']).values.ravel()
GPQM_ncum_5_data = GPQM_ncum_5_data['RAINFALL'].mean(dim=['lat', 'lon']).values.ravel()
GPQM_ncum_3_data = GPQM_ncum_3_data['RAINFALL'].mean(dim=['lat', 'lon']).values.ravel()
'''
obs_data = obs_data["rf"].mean(dim=['time']).values.flatten()
obs_data = obs_data[~np.isnan(obs_data)]
obs_data = obs_data[obs_data > 0]

ncum_1_data = ncum_1_data["APCP_surface"].mean(dim=['time']).values.flatten()
ncum_1_data = ncum_1_data[~np.isnan(ncum_1_data)]
ncum_1_data = ncum_1_data[ncum_1_data > 0]
ncum_5_data = ncum_5_data["APCP_surface"].mean(dim=['time']).values.flatten()
ncum_5_data = ncum_5_data[~np.isnan(ncum_5_data)]
ncum_5_data = ncum_5_data[ncum_5_data > 0]
ncum_3_data = ncum_3_data["APCP_surface"].mean(dim=['time']).values.flatten()
ncum_3_data = ncum_3_data[~np.isnan(ncum_3_data)]
ncum_3_data = ncum_3_data[ncum_3_data > 0]

EQM_ncum_1_data = EQM_ncum_1_data["RAINFALL"].mean(dim=['time']).values.flatten()
EQM_ncum_1_data = EQM_ncum_1_data[~np.isnan(EQM_ncum_1_data)]
EQM_ncum_1_data = EQM_ncum_1_data[EQM_ncum_1_data > 0]
EQM_ncum_5_data = EQM_ncum_5_data["RAINFALL"].mean(dim=['time']).values.flatten()
EQM_ncum_5_data = EQM_ncum_5_data[~np.isnan(EQM_ncum_5_data)]
EQM_ncum_5_data = EQM_ncum_5_data[EQM_ncum_5_data > 0]
EQM_ncum_3_data = EQM_ncum_3_data["RAINFALL"].mean(dim=['time']).values.flatten()
EQM_ncum_3_data = EQM_ncum_3_data[~np.isnan(EQM_ncum_3_data)]
EQM_ncum_3_data = EQM_ncum_3_data[EQM_ncum_3_data > 0]

PQM_ncum_1_data = PQM_ncum_1_data["RAINFALL"].mean(dim=['time']).values.flatten()
PQM_ncum_1_data = PQM_ncum_1_data[~np.isnan(PQM_ncum_1_data)]
PQM_ncum_1_data = PQM_ncum_1_data[PQM_ncum_1_data > 0]
PQM_ncum_5_data = PQM_ncum_5_data["RAINFALL"].mean(dim=['time']).values.flatten()
PQM_ncum_5_data = PQM_ncum_5_data[~np.isnan(PQM_ncum_5_data)]
PQM_ncum_5_data = PQM_ncum_5_data[PQM_ncum_5_data > 0]
PQM_ncum_3_data = PQM_ncum_3_data["RAINFALL"].mean(dim=['time']).values.flatten()
PQM_ncum_3_data = PQM_ncum_3_data[~np.isnan(PQM_ncum_3_data)]
PQM_ncum_3_data = PQM_ncum_3_data[PQM_ncum_3_data > 0]

GPQM_ncum_1_data = GPQM_ncum_1_data["RAINFALL"].mean(dim=['time']).values.flatten()
GPQM_ncum_1_data = GPQM_ncum_1_data[~np.isnan(GPQM_ncum_1_data)]
GPQM_ncum_1_data = GPQM_ncum_1_data[GPQM_ncum_1_data > 0]
GPQM_ncum_5_data = GPQM_ncum_5_data["RAINFALL"].mean(dim=['time']).values.flatten()
PQM_ncum_5_data = GPQM_ncum_5_data[~np.isnan(GPQM_ncum_5_data)]
GPQM_ncum_5_data = GPQM_ncum_5_data[GPQM_ncum_5_data > 0]
GPQM_ncum_3_data = GPQM_ncum_3_data["RAINFALL"].mean(dim=['time']).values.flatten()
GPQM_ncum_3_data = GPQM_ncum_3_data[~np.isnan(GPQM_ncum_3_data)]
GPQM_ncum_3_data = GPQM_ncum_3_data[GPQM_ncum_3_data > 0]

reference_quantiles = np.sort(obs_data)
model_quantiles_1 = np.sort(ncum_1_data)
model_quantiles_5 = np.sort(ncum_5_data)
model_quantiles_3 = np.sort(ncum_3_data)
model_quantiles_EQM_1 = np.sort(EQM_ncum_1_data)
model_quantiles_EQM_3 = np.sort(EQM_ncum_3_data)
model_quantiles_EQM_5 = np.sort(EQM_ncum_5_data)
model_quantiles_PQM_1 = np.sort(PQM_ncum_1_data)
model_quantiles_PQM_3 = np.sort(PQM_ncum_3_data)
model_quantiles_PQM_5 = np.sort(PQM_ncum_5_data)
model_quantiles_GPQM_1 = np.sort(GPQM_ncum_1_data)
model_quantiles_GPQM_3 = np.sort(GPQM_ncum_3_data)
model_quantiles_GPQM_5 = np.sort(GPQM_ncum_5_data)

min_len = min(len(model_quantiles_1), len(reference_quantiles))
model_quantiles_1 = model_quantiles_1[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_5), len(reference_quantiles))
model_quantiles_5 = model_quantiles_5[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_3), len(reference_quantiles))
model_quantiles_3 = model_quantiles_3[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_EQM_1), len(reference_quantiles))
model_quantiles_EQM_1 = model_quantiles_EQM_1[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_EQM_5), len(reference_quantiles))
model_quantiles_EQM_2 = model_quantiles_EQM_5[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_EQM_3), len(reference_quantiles))
model_quantiles_EQM_3 = model_quantiles_EQM_3[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_PQM_1), len(reference_quantiles))
model_quantiles_PQM_1 = model_quantiles_PQM_1[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_PQM_5), len(reference_quantiles))
model_quantiles_PQM_5 = model_quantiles_PQM_5[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_PQM_3), len(reference_quantiles))
model_quantiles_PQM_3 = model_quantiles_PQM_3[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_GPQM_1), len(reference_quantiles))
model_quantiles_GPQM_1 = model_quantiles_GPQM_1[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_GPQM_5), len(reference_quantiles))
model_quantiles_GPQM_5 = model_quantiles_GPQM_5[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_GPQM_3), len(reference_quantiles))
model_quantiles_GPQM_3 = model_quantiles_GPQM_3[:min_len]
reference_quantiles = reference_quantiles[:min_len]

fig, axs = plt.subplots(4, 3, figsize=(15, 10))  # Adjusted to 4 rows and 3 columns for 12 plots

# Plot for model_quantiles_g1
min_len = min(len(model_quantiles_1), len(reference_quantiles))
axs[0, 0].scatter(reference_quantiles[:min_len], model_quantiles_1[:min_len], label="G1", alpha=0.6)
axs[0, 0].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[0, 0].legend()
axs[0, 0].set_title('Q-Q plot: NCUM-G DAY 1 vs Observed')

# Plot for model_quantiles_g2
min_len = min(len(model_quantiles_3), len(reference_quantiles))
axs[0, 1].scatter(reference_quantiles[:min_len], model_quantiles_3[:min_len], label="G2", alpha=0.6)
axs[0, 1].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[0, 1].legend()
axs[0, 1].set_title('Q-Q plot: NCUM-G DAY 3 vs Observed')

# Plot for model_quantiles_g3
min_len = min(len(model_quantiles_5), len(reference_quantiles))
axs[0, 2].scatter(reference_quantiles[:min_len], model_quantiles_5[:min_len], label="G3", alpha=0.6)
axs[0, 2].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[0, 2].legend()
axs[0, 2].set_title('Q-Q plot: NCUM-G DAY 5 vs Observed')

# Repeat for other subplots, adjusting the row and column indices
min_len = min(len(model_quantiles_EQM_1), len(reference_quantiles))
axs[1, 0].scatter(reference_quantiles[:min_len], model_quantiles_EQM_1[:min_len], label="G1", alpha=0.6)
axs[1, 0].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[1, 0].legend()
axs[1, 0].set_title('Q-Q plot: NCUM-G DAY 1 (EQM) vs Observed')

min_len = min(len(model_quantiles_EQM_3), len(reference_quantiles))
axs[1, 1].scatter(reference_quantiles[:min_len], model_quantiles_EQM_3[:min_len], label="G2", alpha=0.6)
axs[1, 1].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[1, 1].legend()
axs[1, 1].set_title('Q-Q plot: NCUM-G DAY 3 (EQM) vs Observed')

min_len = min(len(model_quantiles_EQM_5), len(reference_quantiles))
axs[1, 2].scatter(reference_quantiles[:min_len], model_quantiles_EQM_5[:min_len], label="G3", alpha=0.6)
axs[1, 2].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[1, 2].legend()
axs[1, 2].set_title('Q-Q plot: NCUM-G DAY 5 (EQM) vs Observed')

min_len = min(len(model_quantiles_PQM_1), len(reference_quantiles))
axs[2, 0].scatter(reference_quantiles[:min_len], model_quantiles_PQM_1[:min_len], label="G1", alpha=0.6)
axs[2, 0].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[2, 0].legend()
axs[2, 0].set_title('Q-Q plot: NCUM-G DAY 1 (PQM) vs Observed')

min_len = min(len(model_quantiles_PQM_3), len(reference_quantiles))
axs[2, 1].scatter(reference_quantiles[:min_len], model_quantiles_PQM_3[:min_len], label="G2", alpha=0.6)
axs[2, 1].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[2, 1].legend()
axs[2, 1].set_title('Q-Q plot: NCUM-G DAY 3 (PQM) vs Observed')

min_len = min(len(model_quantiles_PQM_5), len(reference_quantiles))
axs[2, 2].scatter(reference_quantiles[:min_len], model_quantiles_PQM_5[:min_len], label="G3", alpha=0.6)
axs[2, 2].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[2, 2].legend()
axs[2, 2].set_title('Q-Q plot: NCUM-G DAY 5 (PQM) vs Observed')

min_len = min(len(model_quantiles_GPQM_1), len(reference_quantiles))
axs[3, 0].scatter(reference_quantiles[:min_len], model_quantiles_GPQM_1[:min_len], label="G1", alpha=0.6)
axs[3, 0].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[3, 0].legend()
axs[3, 0].set_title('Q-Q plot: NCUM-G DAY 1 (GPQM) vs Observed')

min_len = min(len(model_quantiles_GPQM_3), len(reference_quantiles))
axs[3, 1].scatter(reference_quantiles[:min_len], model_quantiles_GPQM_3[:min_len], label="G2", alpha=0.6)
axs[3, 1].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[3, 1].legend()
axs[3, 1].set_title('Q-Q plot: NCUM-G DAY 3 (GPQM) vs Observed')

min_len = min(len(model_quantiles_GPQM_5), len(reference_quantiles))
axs[3, 2].scatter(reference_quantiles[:min_len], model_quantiles_GPQM_5[:min_len], label="G3", alpha=0.6)
axs[3, 2].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[3, 2].legend()
axs[3, 2].set_title('Q-Q plot: NCUM-G DAY 5 (GPQM) vs Observed')

plt.tight_layout()
plt.show()
