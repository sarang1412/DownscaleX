import xarray as xr
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import geopandas as gpd
from shapely.geometry import mapping
import rioxarray

#data
obs = xr.open_dataset('e:\\Dissertation\\data\\obs_regridded.nc')
ncum_r_1 = xr.open_dataset('e:\\Dissertation\\data\\ncum_r_1_regridded.nc')
ncum_r_2 = xr.open_dataset('e:\\Dissertation\\data\\ncum_r_2_regridded.nc')
ncum_r_3 = xr.open_dataset('e:\\Dissertation\\data\\ncum_r_3_regridded.nc')

obs_data = obs["rf"].mean(dim=['time']).values.flatten()
obs_data = obs_data[~np.isnan(obs_data)]
obs_data = obs_data[obs_data > 0]

ncum_r1_data = ncum_r_1["APCP_surface"].mean(dim=['time']).values.flatten()
ncum_r1_data = ncum_r1_data[~np.isnan(ncum_r1_data)]
ncum_r1_data = ncum_r1_data[ncum_r1_data > 0]
ncum_r2_data = ncum_r_2["APCP_surface"].mean(dim=['time']).values.flatten()
ncum_r2_data = ncum_r2_data[~np.isnan(ncum_r2_data)]
ncum_r2_data = ncum_r2_data[ncum_r2_data > 0]
ncum_r3_data = ncum_r_3["APCP_surface"].mean(dim=['time']).values.flatten()
ncum_r3_data = ncum_r3_data[~np.isnan(ncum_r3_data)]
ncum_r3_data = ncum_r3_data[ncum_r3_data > 0]


print("Maximum ncum_r1", np.max(ncum_r1_data))
print("Maximum ncum_r2", np.max(ncum_r2_data))
print("Maximum ncum_r3", np.max(ncum_r3_data))
print("Maximum obs", np.max(obs_data))

#Q-Q plot
reference_quantiles = np.sort(obs_data)
model_quantiles_r1 = np.sort(ncum_r1_data)
model_quantiles_r2 = np.sort(ncum_r2_data)
model_quantiles_r3 = np.sort(ncum_r3_data)
    
min_len = min(len(model_quantiles_r1), len(reference_quantiles))
model_quantiles_r1 = model_quantiles_r1[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_r2), len(reference_quantiles))
model_quantiles_r2 = model_quantiles_r2[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_r3), len(reference_quantiles))
model_quantiles_r3 = model_quantiles_r3[:min_len]
reference_quantiles = reference_quantiles[:min_len]

# ...existing code...

#Q-Q plot
reference_quantiles = np.sort(obs_data)
model_quantiles_r1 = np.sort(ncum_r1_data)
model_quantiles_r2 = np.sort(ncum_r2_data)
model_quantiles_r3 = np.sort(ncum_r3_data)

# Create subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Plot for model_quantiles_r1
min_len = min(len(model_quantiles_r1), len(reference_quantiles))
axs[0].scatter(reference_quantiles[:min_len], model_quantiles_r1[:min_len], label="R1", alpha=0.6)
axs[0].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[0].legend()
axs[0].set_title('Q-Q plot: Model R1 vs Observed')

# Plot for model_quantiles_r2
min_len = min(len(model_quantiles_r2), len(reference_quantiles))
axs[1].scatter(reference_quantiles[:min_len], model_quantiles_r2[:min_len], label="R2", alpha=0.6)
axs[1].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[1].legend()
axs[1].set_title('Q-Q plot: Model R2 vs Observed')

# Plot for model_quantiles_r3
min_len = min(len(model_quantiles_r3), len(reference_quantiles))
axs[2].scatter(reference_quantiles[:min_len], model_quantiles_r3[:min_len], label="R3", alpha=0.6)
axs[2].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[2].legend()
axs[2].set_title('Q-Q plot: Model R3 vs Observed')

plt.tight_layout()
plt.show()