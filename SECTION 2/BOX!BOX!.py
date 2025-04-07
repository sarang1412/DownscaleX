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

# === Separate Plots for Day 1, Day 3, and Day 5 ===
day1_datasets_flat = []
day3_datasets_flat = []
day5_datasets_flat = []

day1_labels = []
day3_labels = []
day5_labels = []

for month_num, month_name in zip(month_nums, month_names):
    print(f"Processing {month_name}...")

    # Select time for each month
    obs_m       = obs.sel(time=obs.time.dt.month == month_num)
    ncum1_m     = ncum_1.sel(time=ncum_1.time.dt.month == month_num)
    ncum3_m     = ncum_3.sel(time=ncum_3.time.dt.month == month_num)
    ncum5_m     = ncum_5.sel(time=ncum_5.time.dt.month == month_num)
    EQM1_m      = EQM_1.sel(time=EQM_1.time.dt.month == month_num)
    EQM3_m      = EQM_3.sel(time=EQM_3.time.dt.month == month_num)
    EQM5_m      = EQM_5.sel(time=EQM_5.time.dt.month == month_num)
    PQM1_m      = PQM_1.sel(time=PQM_1.time.dt.month == month_num)
    PQM3_m      = PQM_3.sel(time=PQM_3.time.dt.month == month_num)
    PQM5_m      = PQM_5.sel(time=PQM_5.time.dt.month == month_num)
    GPQM1_m     = GPQM_1.sel(time=GPQM_1.time.dt.month == month_num)
    GPQM3_m     = GPQM_3.sel(time=GPQM_3.time.dt.month == month_num)
    GPQM5_m     = GPQM_5.sel(time=GPQM_5.time.dt.month == month_num)

    # Flatten datasets
    obs_data     = clean_flatten(obs_m, "rf")
    ncum1_data   = clean_flatten(ncum1_m, "APCP_surface")
    ncum3_data   = clean_flatten(ncum3_m, "APCP_surface")
    ncum5_data   = clean_flatten(ncum5_m, "APCP_surface")
    EQM1_data    = clean_flatten(EQM1_m, "RAINFALL")
    EQM3_data    = clean_flatten(EQM3_m, "RAINFALL")
    EQM5_data    = clean_flatten(EQM5_m, "RAINFALL")
    PQM1_data    = clean_flatten(PQM1_m, "RAINFALL")
    PQM3_data    = clean_flatten(PQM3_m, "RAINFALL")
    PQM5_data    = clean_flatten(PQM5_m, "RAINFALL")
    GPQM1_data   = clean_flatten(GPQM1_m, "RAINFALL")
    GPQM3_data   = clean_flatten(GPQM3_m, "RAINFALL")
    GPQM5_data   = clean_flatten(GPQM5_m, "RAINFALL")

    # Append data for Day 1, Day 3, and Day 5
    day1_datasets_flat.extend([obs_data, ncum1_data, EQM1_data, PQM1_data, GPQM1_data])
    day3_datasets_flat.extend([obs_data, ncum3_data, EQM3_data, PQM3_data, GPQM3_data])
    day5_datasets_flat.extend([obs_data, ncum5_data, EQM5_data, PQM5_data, GPQM5_data])

    day1_labels.extend([
        f"OBS ({month_name})", f"Raw D1 ({month_name})", f"EQM D1 ({month_name})",
        f"PQM D1 ({month_name})", f"GPQM D1 ({month_name})"
    ])
    day3_labels.extend([
        f"OBS ({month_name})", f"Raw D3 ({month_name})", f"EQM D3 ({month_name})",
        f"PQM D3 ({month_name})", f"GPQM D3 ({month_name})"
    ])
    day5_labels.extend([
        f"OBS ({month_name})", f"Raw D5 ({month_name})", f"EQM D5 ({month_name})",
        f"PQM D5 ({month_name})", f"GPQM D5 ({month_name})"
    ])

# Define updated colors for each dataset type
colors = {
    "OBS": "#1f77b4",  # Blue
    "Raw": "#ff7f0e",  # Orange
    "EQM": "#2ca02c",  # Green
    "PQM": "#d62728",  # Red
    "GPQM": "#9467bd"  # Purple
}

# Helper function to create grouped box plots with adjusted text and plot size
def plot_grouped_boxplot(datasets_flat, labels, title):
    plt.figure(figsize=(16, 7))  
    positions = []
    color_map = []
    group_width = 5  # Number of boxes per month
    gap = 1  # Gap between groups

    for i, label in enumerate(labels):
        month_index = i // group_width
        box_index = i % group_width
        positions.append(month_index * (group_width + gap) + box_index)
        if "OBS" in label:
            color_map.append(colors["OBS"])
        elif "Raw" in label:
            color_map.append(colors["Raw"])
        elif "EQM" in label:
            color_map.append(colors["EQM"])
        elif "PQM" in label:
            color_map.append(colors["PQM"])
        elif "GPQM" in label:
            color_map.append(colors["GPQM"])

    bp = plt.boxplot(datasets_flat, positions=positions, patch_artist=True, showfliers=False)
    for patch, color in zip(bp['boxes'], color_map):
        patch.set_facecolor(color)

    # Add legend
    legend_handles = [
        plt.Line2D([0], [0], color=color, lw=4, label=label)
        for label, color in colors.items()
    ]
    plt.legend(handles=legend_handles, title="Dataset Type", loc="upper right", fontsize=14, title_fontsize=16)

    plt.xticks(
        ticks=[(i * (group_width + gap) + (group_width - 1) / 2) for i in range(len(month_names))],
        labels=month_names,
        fontsize=14  # Increase font size for x-axis labels
    )
    plt.ylabel("Daily Rainfall (mm)", fontsize=16)  # Increase font size for y-axis label
    plt.title(title, fontsize=18)  # Increase font size for title
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

# Plot Day 1
plot_grouped_boxplot(
    day1_datasets_flat,
    day1_labels,
    "JJAS Daily Rainfall Distribution - Day 1 (2020–2024, India Region)"
)

# Plot Day 3
plot_grouped_boxplot(
    day3_datasets_flat,
    day3_labels,
    "JJAS Daily Rainfall Distribution - Day 3 (2020–2024, India Region)"
)

# Plot Day 5
plot_grouped_boxplot(
    day5_datasets_flat,
    day5_labels,
    "JJAS Daily Rainfall Distribution - Day 5 (2020–2024, India Region)"
)
