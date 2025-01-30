import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import dask.array as da

#load data
ncum_g = xr.open_dataset('e:\\Dissertation\\data\\ncumg_day1rf_jjas2020-24.nc', chunks={'time': 10})
ncum_r = xr.open_dataset('e:\\Dissertation\\ncumr_day1rf_jjas2023.nc')
obs = xr.open_dataset('e:\\Dissertation\\data\\IMD_MSG-2020-24-jjas.nc', )

#change lat and lon to match
ncum_g = ncum_g.rename({'latitude': 'lat', 'longitude': 'lon'})
#ncum_r = ncum_r.rename({'latitude': 'lat', 'longitude': 'lon'})

#regrid
ncum_g_regridded = ncum_g.interp_like(obs, method='nearest')
ncum_r_regridded = ncum_r.interp_like(obs, method='nearest')

#lat and lon slice
obs = obs.sel(lat=slice(6, 41), lon=slice(62, 106))
ncum_g_regridded = ncum_g_regridded.sel(lat=slice(6, 41), lon=slice(62, 106))
ncum_r_regridded = ncum_r_regridded.sel(lat=slice(6, 41), lon=slice(62, 106))

#proper date format
obs['time'] = pd.to_datetime(obs['time'].values)
ncum_g_regridded['time'] = pd.to_datetime(ncum_g_regridded['time'].values)
ncum_r_regridded['time'] = pd.to_datetime(ncum_r_regridded['time'].values)

# daily averages for each day of JJAS months over the years 2020-2024
daily_avg_obs = obs.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_g = ncum_g_regridded.groupby('time.dayofyear').mean(dim='time')
daily_avg_ncum_r = ncum_r_regridded.groupby('time.dayofyear').mean(dim='time')

#variance
obs_variance = daily_avg_obs['rf'].var(dim=('lat', 'lon'))
#daily_avg_ncum_g = daily_avg_ncum_g.dropna(dim='dayofyear', how='all')
ncum_g_variance = daily_avg_ncum_g['APCP_surface'].var(dim=('lat', 'lon'))
ncum_r_variance = daily_avg_ncum_r['APCP_24'].var(dim=('lat', 'lon'))

# Plot
data_variances = [
    (obs_variance, 'Observed data 10Km'),
    (ncum_g_variance, 'NCUM-G 12Km to 10Km (day 1 forecast)'),
    (ncum_r_variance, 'NCUM-R 4km to 10Km (day 1 forecast)')
]

for data, label in data_variances:
    plt.plot(data['dayofyear'], data, label=label)
    plt.xlabel('Time')
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%B'))
    plt.ylabel('Variance')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()