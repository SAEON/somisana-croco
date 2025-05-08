#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 11:14:30 2024
@author: nkululeko
"""
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

dataset_path = "/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/GLORYS/Validation/"
ds = xr.open_mfdataset(dataset_path + "GLORYS_Y2009M01_Y2013M12.nc", combine='by_coords')

BettysBay_lat = -34.817
BettysBay_lon = 18.922

Doringbaai_lat = -31.817
Doringbaai_lon = 18.231

Time = ds.time.values
BettysBay = ds['temp'].sel(latT=BettysBay_lat, lonT=BettysBay_lon,method="nearest").values[:,0]
Doringbaai = ds['temp'].sel(latT=Doringbaai_lat, lonT=Doringbaai_lon,method="nearest").values[:,0]

plt.figure(figsize=(10, 6))
plt.plot(Time, BettysBay, marker='', linestyle='--', label=f'BettysBay=({BettysBay_lat}:{BettysBay_lon})')
plt.plot(Time, Doringbaai, marker='', linestyle='-', label=f'Doringbaai=({Doringbaai_lat}:{Doringbaai_lon})')
plt.title(f"GLORYS Model Surface Temperature Time Series \n at Betts Bay ({BettysBay_lat} S: {BettysBay_lon} E) vs Doringbaai ({Doringbaai_lat} S: {Doringbaai_lon} E)", fontsize=16)
plt.xlabel("Time")
plt.ylabel("Temperature (°C)")
plt.grid(True)
plt.legend(loc="upper right")
plt.show()

#%% Do the same to croco

# List of filenames
filenames = [
    "Doringbaai_DAFF.nc",
    "Sea Point_SAWS.nc"
]

ERA5_VAL_path = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/AJSMIT_UTR/'
WASA3_VAL_path = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/Validation/AJSMIT_UTR/'

savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/'
# plt.rc('text', usetex=True)

for filename in filenames:
    ds_ERA5 = xr.open_dataset(ERA5_VAL_path + "Validated_" + filename)
    data_obs_model_timeaxis_ERA5 = ds_ERA5.insitu_data_temp.squeeze()    
    time_variable = ds_ERA5.time
    data_model = ds_ERA5.model_data_temp.squeeze()
    # data_obs_model_timeaxis_ERA5 = ds_ERA5.insitu_data_temp.squeeze()

    ds_WASA3 = xr.open_dataset(WASA3_VAL_path + "Validated_" + filename)
    data_obs_model_timeaxis_WASA3 = ds_WASA3.insitu_data_temp.squeeze()    
    time_variable_WASA3 = ds_WASA3.time
    data_model_WASA3 = ds_WASA3.model_data_temp.squeeze()
    # data_obs_model_timeaxis_WASA3 = ds_WASA3.insitu_data_temp.squeeze()

    lat_insitu = ds_ERA5.latitude
    lat_insitu = np.round(lat_insitu, 3)
    lon_insitu = ds_ERA5.longitude
    lon_insitu = np.round(lon_insitu, 3)

    lat_model = ds_ERA5.latitude
    lat_model = np.round(lat_model, 3)
    lon_model = ds_ERA5.longitude
    lon_model = np.round(lon_model, 3)
    
    plt.figure(figsize=(10, 6))
    plt.plot(Time, data_model, marker='', linestyle='--', label=f'SeaPoint')
    plt.plot(Time, data_model, marker='', linestyle='--', label=f'Doringbaai')
    # plt.plot(Time, Doringbaai, marker='', linestyle='-', label=f'Doringbaai=({Doringbaai_lat}:{Doringbaai_lon})')
    plt.title(f"GLORYS Model Surface Temperature Time Series \n at Betts Bay ", fontsize=16)
    plt.xlabel("Time")
    plt.ylabel("Temperature (°C)")
    plt.grid(True)
    plt.legend(loc="upper right")
    plt.show()

#%% 

import xarray as xr
import os
from tqdm import tqdm

# Define the dataset path and output path
dataset_path = "/mnt/e/GLORYS_SA_EEZ_ZOOM/"
model_path = "/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/GLORYS/Validation/"

# Create a list to store all file patterns for the years 2009 to 2013
file_patterns = [dataset_path + f"GLORYS_Y{year}M{month}.nc" for year in range(2009, 2014) for month in range(1, 13)]

# Filter out any non-existent files to avoid errors
existing_files = [f for f in file_patterns if os.path.exists(f)]

# Load all files for the specified years with a progress bar
print("Loading datasets...")
ds = xr.open_mfdataset(tqdm(existing_files), combine='by_coords')

# Extract the temperature data over the specified period
temperature_data = ds.temp

# Save the extracted temperature data to a new NetCDF file
output_file = model_path + "GLORYS_Y2009M01_Y2013M12.nc"
print(f"Saving combined dataset to {output_file}...")
temperature_data.to_netcdf(output_file)

print(f"Combined temperature dataset saved to {output_file}")


