import xarray as xr
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import geopandas as gpd
from shapely.geometry import mapping
import rioxarray

#data
obs = xr.open_dataset('e:\\Dissertation\\data\\obs_regridded.nc')

ncum_g_1 = xr.open_dataset('e:\\Dissertation\\data\\ncum_gg_1_regridded.nc')
ncum_g_2 = xr.open_dataset('e:\\Dissertation\\data\\ncum_g_2_regridded.nc')
ncum_g_3 = xr.open_dataset('e:\\Dissertation\\data\\ncum_g_3_regridded.nc')

obs_data = obs["rf"].mean(dim=['time']).values.flatten()
obs_data = obs_data[~np.isnan(obs_data)]
obs_data = obs_data[obs_data > 0]

ncum_g1_data = ncum_g_1["APCP_surface"].mean(dim=['time']).values.flatten()
ncum_g1_data = ncum_g1_data[~np.isnan(ncum_g1_data)]
ncum_g1_data = ncum_g1_data[ncum_g1_data > 0]
ncum_g2_data = ncum_g_2["APCP_surface"].mean(dim=['time']).values.flatten()
ncum_g2_data = ncum_g2_data[~np.isnan(ncum_g2_data)]
ncum_g2_data = ncum_g2_data[ncum_g2_data > 0]
ncum_g3_data = ncum_g_3["APCP_surface"].mean(dim=['time']).values.flatten()
ncum_g3_data = ncum_g3_data[~np.isnan(ncum_g3_data)]
ncum_g3_data = ncum_g3_data[ncum_g3_data > 0]

print("Maximum obs", np.max(obs_data))
print("Maximum ncum_g1", np.max(ncum_g1_data))
print("Maximum ncum_g2", np.max(ncum_g2_data))
print("Maximum ncum_g3", np.max(ncum_g3_data))

reference_quantiles = np.sort(obs_data)
model_quantiles_g1 = np.sort(ncum_g1_data)
model_quantiles_g2 = np.sort(ncum_g2_data)
model_quantiles_g3 = np.sort(ncum_g3_data)

min_len = min(len(model_quantiles_g1), len(reference_quantiles))
model_quantiles_g1 = model_quantiles_g1[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_g2), len(reference_quantiles))
model_quantiles_g2 = model_quantiles_g2[:min_len]
reference_quantiles = reference_quantiles[:min_len]

min_len = min(len(model_quantiles_g3), len(reference_quantiles))
model_quantiles_g3 = model_quantiles_g3[:min_len]
reference_quantiles = reference_quantiles[:min_len]

reference_quantiles = np.sort(obs_data)
model_quantiles_g1 = np.sort(ncum_g1_data)
model_quantiles_g2 = np.sort(ncum_g2_data)
model_quantiles_g3 = np.sort(ncum_g3_data)


fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Plot for model_quantiles_g1
min_len = min(len(model_quantiles_g1), len(reference_quantiles))
axs[0].scatter(reference_quantiles[:min_len], model_quantiles_g1[:min_len], label="G1", alpha=0.6)
axs[0].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[0].legend()
axs[0].set_title('Q-Q plot: Model G1 vs Observed')

# Plot for model_quantiles_g2
min_len = min(len(model_quantiles_g2), len(reference_quantiles))
axs[1].scatter(reference_quantiles[:min_len], model_quantiles_g2[:min_len], label="G2", alpha=0.6)
axs[1].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[1].legend()
axs[1].set_title('Q-Q plot: Model G2 vs Observed')

# Plot for model_quantiles_g3
min_len = min(len(model_quantiles_g3), len(reference_quantiles))
axs[2].scatter(reference_quantiles[:min_len], model_quantiles_g3[:min_len], label="G3", alpha=0.6)
axs[2].plot(reference_quantiles[:min_len], reference_quantiles[:min_len], 'r--', label='1:1 Line')
axs[2].legend()
axs[2].set_title('Q-Q plot: Model G3 vs Observed')

plt.tight_layout()
plt.show()