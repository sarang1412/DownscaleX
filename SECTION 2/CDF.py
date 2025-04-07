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
ncum_1 = ncum_1.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_3 = ncum_3.rename({'latitude': 'lat', 'longitude': 'lon'})
ncum_5 = ncum_5.rename({'latitude': 'lat', 'longitude': 'lon'})
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

# === Extract and clean data for analysis ===
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

def plot_cdf(data, label, colour):
    sorted_data = np.sort(data)
    cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
    plt.plot(sorted_data, cdf, label=label, color=colour)

plt.figure(figsize=(8, 6))
plot_cdf(obs_data, "Observed (IMD-MSG)", "black")
plot_cdf(ncum1_data, "Raw NCUM-D1", "red")
plot_cdf(EQM1_data, "EQM D1", "blue")
plot_cdf(PQM1_data, "PQM D1", "green")
plot_cdf(GPQM1_data, "GPQM D1", "orange")

plt.xlabel("Rainfall (mm)")
plt.ylabel("Cumulative Probability")
plt.title("CDF Plot - Day 1 Forecast vs Observations")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
