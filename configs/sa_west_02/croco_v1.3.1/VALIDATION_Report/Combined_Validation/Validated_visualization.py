#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 15:32:55 2024

@author: nkululeko
"""
from scipy.interpolate import griddata
import glob
from pandas.plotting import table
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from math import sqrt
import cartopy.crs as ccrs
import calendar
import cartopy.feature as cfeature
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import crocotools_py.validation as val


# filenames_with_shifts = [
#     ["Betty's Bay_DAFF.nc","0","0"],
#     ["Bordjies Deep_DAFF.nc","0","0"],
#     ["Bordjies_DAFF.nc","0","0"], #-4
#     ["Doringbaai_DAFF.nc","0","0"],
#     ["Fish Hoek_SAWS.nc","2","0"],
#     ["Gordons Bay_SAWS.nc","0","0"],
#     ["Kalk Bay_SAWS.nc","2","0"],
#     ["Kommetjie_SAWS.nc","0","0"],
#     ["Lamberts Bay_SAWS.nc","0","0"],
#     ["Muizenberg_SAWS.nc","0","0"],
#     ["Oudekraal_DAFF.nc","0","0"],
#     ["Saldanha Bay_SAWS.nc","0","-2"],
#     ["Sea Point_SAWS.nc","0","0"],
#     ["St Helena Bay_SAWS.nc","0","0"],
#     ["Yzerfontein_SAWS.nc","0","0"]
# ]

filenames_with_shifts = [
    # ["WQM_20m_TS.nc","0","0"],
    ["WQM_70m_TS.nc", "0", "0"]
]

ERA5_VAL_path = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/GPITCHER_SHB/'
WASA3_VAL_path = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/Validation/GPITCHER_SHB/'

# ERA5_VAL_path = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/AJSMIT_UTR/'
# WASA3_VAL_path = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/Validation/AJSMIT_UTR/'
GLORYS_path = "/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/GLORYS/Validation/"
ds_GLORYS = xr.open_mfdataset(
    GLORYS_path + "GLORYS_Y2009M01_Y2013M12.nc",
    combine='by_coords'
)

savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/'
# plt.rc('text', usetex=True)

# Subset the spatial domain
lat_bounds = (-36, -31)
lon_bounds = (14, 20)
ds_GLORYS = ds_GLORYS.sel(latT=slice(
    lat_bounds[0], lat_bounds[1]), lonT=slice(lon_bounds[0], lon_bounds[1]))

# Reduce precision to save memory
ds_GLORYS['temp'] = ds_GLORYS['temp'].astype(np.float64)


for file_info in filenames_with_shifts:
    filename, i_moved, j_moved = file_info
    ds_ERA5 = xr.open_dataset(ERA5_VAL_path + "Validated_" + filename)
    time_variable = ds_ERA5.time
    data_model = ds_ERA5.model_data_temp.squeeze()
    data_obs_model_timeaxis_ERA5 = ds_ERA5.insitu_data_temp.squeeze()

    # Check if the dataset's time range exceeds the desired range
    timeSlice = slice("2009-01-01", "2013-12-31")
    if time_variable.min() < np.datetime64("2009-01-01") or time_variable.max() > np.datetime64("2013-12-31"):
        # Subset the data to the desired time range (2009-01-01 to 2013-12-31)
        data_obs_model_timeaxis_ERA5 = data_obs_model_timeaxis_ERA5.sel(
            time=timeSlice)
        time_variable = ds_ERA5.time.sel(time=timeSlice)
        data_model = data_model.sel(time=timeSlice)
    else:
        # No subset needed, use the original data
        data_obs_model_timeaxis_ERA5 = data_obs_model_timeaxis_ERA5
        time_variable = ds_ERA5.time

    ds_WASA3 = xr.open_dataset(WASA3_VAL_path + "Validated_" + filename)
    data_obs_model_timeaxis_WASA3 = ds_WASA3.insitu_data_temp.squeeze()
    time_variable_WASA3 = ds_WASA3.time
    data_model_WASA3 = ds_WASA3.model_data_temp.squeeze()
    data_obs_model_timeaxis_WASA3 = ds_WASA3.insitu_data_temp.squeeze()

    lat_insitu = ds_ERA5.latitude
    lat_insitu = np.round(lat_insitu, 3)
    lon_insitu = ds_ERA5.longitude
    lon_insitu = np.round(lon_insitu, 3)

    lat_model = ds_ERA5.latitude
    lat_model = np.round(lat_model, 3)
    lon_model = ds_ERA5.longitude
    lon_model = np.round(lon_model, 3)

    idx_i = ds_GLORYS["lonT"][1]-ds_GLORYS["lonT"][0]
    i = idx_i.values*int(i_moved)

    idx_j = ds_GLORYS["latT"][1]-ds_GLORYS["latT"][0]
    j = idx_j.values*int(j_moved)

    Temps_GLORYS = ds_GLORYS['temp'].sel(
        latT=lat_model.values + j, lonT=lon_model.values + i, method="nearest")[:, 0]
    Temps_GLORYS = Temps_GLORYS.squeeze()

    Time_GLORYS = ds_GLORYS.time.values
    data_GLORYS = Temps_GLORYS.squeeze()
    # data_GLORYS = ds_GLORYS.temp[:,0,int(lat_model),int(lon_model)].squeeze()

    model_ano = data_model.groupby(
        "time.month") - data_model.groupby("time.month").mean("time")
    model_clim = data_model.groupby("time.month").mean("time")
    model_ano_WASA3 = data_model_WASA3.groupby(
        "time.month") - data_model_WASA3.groupby("time.month").mean("time")
    model_clim_WASA3 = data_model_WASA3.groupby("time.month").mean("time")

    GLORYS_ano = data_GLORYS.groupby(
        "time.month") - data_GLORYS.groupby("time.month").mean("time")
    GLORYS_clim = data_GLORYS.groupby("time.month").mean("time")

    obs_ano = data_obs_model_timeaxis_ERA5.groupby(
        "time.month") - data_obs_model_timeaxis_ERA5.groupby("time.month").mean("time")
    obs_clim = data_obs_model_timeaxis_ERA5.groupby("time.month").mean("time")
    obs_ano_WASA3 = data_obs_model_timeaxis_WASA3.groupby(
        "time.month") - data_obs_model_timeaxis_WASA3.groupby("time.month").mean("time")
    obs_clim_WASA3 = data_obs_model_timeaxis_WASA3.groupby(
        "time.month").mean("time")

    temp_corr = data_obs_model_timeaxis_ERA5[~np.isnan(
        data_obs_model_timeaxis_ERA5)]
    model_corr = data_model[~np.isnan(data_obs_model_timeaxis_ERA5)]
    GLORYS_corr = data_GLORYS.values[~np.isnan(
        data_obs_model_timeaxis_ERA5.values)]
    # GLORYS_CROCO_corr = data_GLORYS[~np.isnan(data_model)]

    temp_corr_WASA3 = data_obs_model_timeaxis_WASA3[~np.isnan(
        data_obs_model_timeaxis_WASA3)]
    model_corr_WASA3 = data_model_WASA3[~np.isnan(
        data_obs_model_timeaxis_WASA3)]
    GLORYS_corr_WASA3 = data_GLORYS.values[~np.isnan(
        data_obs_model_timeaxis_WASA3.values)]
    # GLORYS_CROCO_corr = data_GLORYS[~np.isnan(data_model)]

    model_ERA5_bias = np.mean(np.abs(model_corr - temp_corr))  # ERA5 bias
    model_in_situ_corr = np.corrcoef(
        np.array(temp_corr), np.array(model_corr.values))[1][0]
    if model_ERA5_bias < 0:
        model_in_situ_rmse = sqrt(mean_squared_error(
            temp_corr, model_corr + model_ERA5_bias))
    else:
        model_in_situ_rmse = sqrt(mean_squared_error(
            temp_corr, model_corr - model_ERA5_bias))

    GLORYS_in_situ_corr = np.corrcoef(
        np.array(temp_corr.values), np.array(GLORYS_corr))[1][0]

    model_GLORYS_bias = np.mean(
        np.abs(GLORYS_corr - temp_corr.values))  # GLORYS bias

    if np.all(np.isnan(GLORYS_corr)):
        GLORYS_in_situ_rmse = np.nan
    elif model_GLORYS_bias < 0:
        GLORYS_in_situ_rmse = sqrt(
            mean_squared_error(temp_corr.values, GLORYS_corr + model_GLORYS_bias))
    else:
        GLORYS_in_situ_rmse = sqrt(
            mean_squared_error(temp_corr.values, GLORYS_corr - model_GLORYS_bias))

    model_WASA3_bias = np.mean(
        np.abs(model_corr_WASA3 - temp_corr_WASA3))  # WASA3 bias
    model_in_situ_corr_WASA3 = np.corrcoef(
        np.array(temp_corr_WASA3), np.array(model_corr_WASA3.values))[1][0]
    if model_WASA3_bias < 0:
        model_in_situ_rmse_WASA3 = sqrt(
            mean_squared_error(temp_corr_WASA3, model_corr_WASA3 + model_WASA3_bias.values))
    else:
        model_in_situ_rmse_WASA3 = sqrt(
            mean_squared_error(temp_corr_WASA3, model_corr_WASA3 - model_WASA3_bias.values))

    # model_in_situ_rmse_WASA3 = sqrt(
    #     mean_squared_error(temp_corr_WASA3, model_corr_WASA3 - model_WASA3_bias.values))
    GLORYS_in_situ_corr_WASA3 = np.corrcoef(
        np.array(temp_corr_WASA3.values), np.array(GLORYS_corr_WASA3))[1][0]

    obs_ano_reshaped = obs_ano.values.reshape(model_ano.shape)

    # Prepare data so that nans are removed from insitu
    temp_ano_corr = obs_ano.values[~np.isnan(obs_ano.values)]
    # Prep data so that CROCO indexes where insitu has gaps are removed
    model_ano_corr = model_ano.values[~np.isnan(obs_ano_reshaped)]
    # Prep data so that GLORYS indexes where insitu has gaps are removed
    GLORYS_ano_corr = GLORYS_ano.values[~np.isnan(obs_ano_reshaped)]
    model_in_situ_ano_corr = np.corrcoef(temp_ano_corr, model_ano_corr)[
        1][0]  # Compare CROCO with insitu correlations
    model_in_situ_ano_rmse = sqrt(
        mean_squared_error(temp_ano_corr, model_ano_corr))
    GLORYS_in_situ_ano_corr = np.corrcoef(temp_ano_corr, GLORYS_ano_corr)[
        1][0]  # Compare GLORYS with insitu correlations

    if np.all(np.isnan(GLORYS_ano_corr)):
        GLORYS_in_situ_ano_rmse = np.nan
    else:
        GLORYS_in_situ_ano_rmse = sqrt(mean_squared_error(
            temp_ano_corr, GLORYS_ano_corr))  # verified to be correct

    # Prepare data so that nans are removed from insitu
    temp_ano_corr_WASA3 = obs_ano.values[~np.isnan(obs_ano.values)]
    # Prep data so that CROCO indexes where insitu has gaps are removed
    model_ano_corr_WASA3 = model_ano_WASA3.values[~np.isnan(obs_ano_reshaped)]
    # Prep data so that GLORYS indexes where insitu has gaps are removed
    GLORYS_ano_corr_WASA3 = GLORYS_ano.values[~np.isnan(obs_ano_reshaped)]
    model_in_situ_ano_corr_WASA3 = np.corrcoef(temp_ano_corr_WASA3, model_ano_corr_WASA3)[
        1][0]  # Compare CROCO with insitu correlations
    model_in_situ_ano_rmse_WASA3 = sqrt(mean_squared_error(
        temp_ano_corr_WASA3, model_ano_corr_WASA3))
    GLORYS_in_situ_ano_corr_WASA3 = np.corrcoef(temp_ano_corr_WASA3, GLORYS_ano_corr_WASA3)[
        1][0]  # Compare GLORYS with insitu correlations

    fig, ax = plt.subplots(figsize=(12, 10), nrows=3, ncols=1)

    month_names = [calendar.month_abbr[i] for i in range(1, 13)]

    if len(obs_clim['month']) != len(month_names):
        months_in_ObsData = (obs_clim['month'] + 1).values.tolist()
        obs_clim['month'] = [month_names[i - 1] for i in months_in_ObsData]
        model_clim['month'] = [month_names[i - 1] for i in months_in_ObsData]
        model_clim_WASA3['month'] = [month_names[i - 1]
                                     for i in months_in_ObsData]
    else:
        obs_clim['month'] = month_names
        model_clim['month'] = month_names
        model_clim_WASA3['month'] = month_names

    plt.title(filename, fontsize=22)

    ax[0].plot(time_variable, data_obs_model_timeaxis_ERA5,
               label=f'{filename[:-3]} in situ', color='blue')
    ax[0].plot(data_model.time, data_model, '--', label='SA_West ERA5_Run',
               color='red', linewidth=2.5, linestyle='-')
    ax[0].plot(data_model_WASA3.time, data_model_WASA3, '--',
               label='SA_West WASA3_Run', color='green', linewidth=2.5, linestyle='-')
    ax[0].plot(Time_GLORYS, Temps_GLORYS, '--', label='GLORYS',
               color='black', linewidth=2.5, linestyle='-')
    # Using multiple text calls to assign different colors
    # Using multiple text calls to assign different colors with a transparent background
    bbox_props = dict(boxstyle="round,pad=0.3",
                      facecolor="white", alpha=0.5, edgecolor="none")

    ax[0].text(0.02, 0.1, '(a)', fontsize=16, color='black',
               transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.06, 0.9, 'r=', fontsize=16, color='black',
               transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.09, 0.9, f'{np.round(model_in_situ_corr, 2)}', fontsize=16,
               color='red', transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.14, 0.9, ';', fontsize=16, color='black',
               transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.16, 0.9, f'{np.round(model_in_situ_corr_WASA3, 2)}',
               fontsize=16, color='green', transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.21, 0.9, ';', fontsize=16, color='black',
               transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.23, 0.9, f'{np.round(GLORYS_in_situ_corr, 2)}', fontsize=16,
               color='black', transform=ax[0].transAxes, bbox=bbox_props)

    ax[0].text(0.01, 0.8, 'RMSE=', fontsize=16, color='black',
               transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.09, 0.8, f'{np.round(model_in_situ_rmse, 2)}', fontsize=16,
               color='red', transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.14, 0.8, ';', fontsize=16, color='black',
               transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.16, 0.8, f'{np.round(model_in_situ_rmse_WASA3, 2)}',
               fontsize=16, color='green', transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.21, 0.8, ';', fontsize=16, color='black',
               transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].text(0.23, 0.8, f'{np.round(GLORYS_in_situ_rmse, 2)}', fontsize=16,
               color='black', transform=ax[0].transAxes, bbox=bbox_props)
    ax[0].set_title('In situ vs Model SST timeseries', fontsize=18)
    ax[0].legend(loc="upper right")
    ax[0].set_xlabel('Time (days)')
    ax[0].set_ylabel('SST (°C)')

    ax[1].plot(time_variable, obs_ano_reshaped,
               label=f'{filename[:-3]} in situ', color='blue', linewidth=2.5, linestyle='--')
    ax[1].plot(data_model.time, model_ano, '--', label='SA_West ERA5_Run',
               color='red', linewidth=2.5, linestyle='--')
    ax[1].plot(data_model_WASA3.time, model_ano_WASA3, '--',
               label='SA_West WASA3_Run', color='green', linewidth=2.5, linestyle='--')
    ax[1].plot(data_model.time, GLORYS_ano, '--', label='GLORYS',
               color='black', linewidth=2.5, linestyle='--')
    # Using multiple text calls to assign different colors with a transparent background
    bbox_props = dict(boxstyle="round,pad=0.3",
                      facecolor="white", alpha=0.5, edgecolor="none")
    ax[1].text(0.02, 0.1, '(b)', fontsize=16, color='black',
               transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.06, 0.90, 'r=', fontsize=16, color='black',
               transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.09, 0.9, f'{np.round(model_in_situ_ano_corr, 2)}',
               fontsize=16, color='red', transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.14, 0.9, ';', fontsize=16, color='black',
               transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.16, 0.9, f'{np.round(model_in_situ_ano_corr_WASA3, 2)}',
               fontsize=16, color='green', transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.21, 0.9, ';', fontsize=16, color='black',
               transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.23, 0.9, f'{np.round(GLORYS_in_situ_ano_corr, 2)}',
               fontsize=16, color='black', transform=ax[1].transAxes, bbox=bbox_props)

    ax[1].text(0.01, 0.8, 'RMSE=', fontsize=16, color='black',
               transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.09, 0.8, f'{np.round(model_in_situ_ano_rmse, 2)}',
               fontsize=16, color='red', transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.14, 0.8, ';', fontsize=16, color='black',
               transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.16, 0.8, f'{np.round(model_in_situ_ano_rmse_WASA3, 2)}',
               fontsize=16, color='green', transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.21, 0.8, ';', fontsize=16, color='black',
               transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].text(0.23, 0.8, f'{np.round(GLORYS_in_situ_ano_rmse, 2)}',
               fontsize=16, color='black', transform=ax[1].transAxes, bbox=bbox_props)
    ax[1].set_title('In situ vs Model SST anomaly', fontsize=18)
    ax[1].legend(loc="upper right")
    ax[1].set_xlabel('Time (days)')
    ax[1].set_ylabel('SST Anomaly')

    ax[2].plot(obs_clim.month, obs_clim,
               label=f'{filename[:-3]} in situ', color='blue', linewidth=2.5, linestyle='-')
    ax[2].plot(model_clim.month, model_clim, '--',
               label='SA_West ERA5_Run', color='red', linewidth=2.5, linestyle='-')
    ax[2].plot(model_clim_WASA3.month, model_clim_WASA3, '--',
               label='SA_West WASA3_Run', color='green', linewidth=2.5, linestyle='-')
    ax[2].plot(model_clim.month, GLORYS_clim, '--', label='GLORYS',
               color='black', linewidth=2.5, linestyle='-')
    ax[2].text(0.02, 0.1, '(c)', fontsize=16, color='black',
               transform=ax[2].transAxes, bbox=bbox_props)
    ax[2].set_title('In situ vs Model SST climatology ', fontsize=18)
    ax[2].legend(loc="upper right")
    ax[2].set_xlabel('Time (months)')
    ax[2].set_ylabel('SST (°C)')

    fig.tight_layout()
    plt.savefig(savepath + "Validated_" + filename[:-3] + '_.png')

    print("---------------------------------------------------------------------------")
    print("---------------------------------------------------------------------------")
    print(
        f'SST validation of the SA_West Model, both the ERA5 and WASA3 wind forced configurations at {filename[:-3]} station located at {lat_model.values}: {lon_model.values}, depth = {-ds_ERA5.depth.values} m and a comparison with GLORYS global model. Subplot (a) shows sst time series, (b) shows anomaly series and (c) shows the climatology of the in situ sst vs the model outputs')
    print("---------------------------------------------------------------------------")
    print("---------------------------------------------------------------------------")

    savename_tbl = 'stats_table.png'

    (insitu_std_G, insitu_correlation_G, insitu_rmse_G, insitu_mean_diff, insitu_min_value_G, insitu_max_value_G,
     model_std_G, model_correlation_G, model_rmse_G, model_mean_diff_G, model_min_value_G, model_max_value_G,
     model_total_bias_G) = val.statistics(Temps_GLORYS, data_obs_model_timeaxis_WASA3)

    (insitu_std_ERA5, insitu_correlation_ERA5, insitu_rmse_ERA5, insitu_mean_diff_ERA5, insitu_min_value_ERA5,
     insitu_max_value_ERA5, model_std_ERA5, model_correlation_ERA5, model_rmse_ERA5, model_mean_diff_ERA5,
     model_min_value_ERA5, model_max_value_ERA5,
     model_total_bias_ERA5) = val.statistics(data_model, data_obs_model_timeaxis_WASA3)

    (insitu_std_WASA3, insitu_correlation_WASA3, insitu_rmse_WASA3, insitu_mean_diff_WASA3,
     insitu_min_value_WASA3, insitu_max_value_WASA3, model_std_WASA3, model_correlation_WASA3,
     model_rmse_WASA3, model_mean_diff_WASA3, model_min_value_WASA3, model_max_value_WASA3,
     model_total_bias_WASA3) = val.statistics(data_model_WASA3, data_obs_model_timeaxis_WASA3)

    # Call the statistics function
    seasonal_stats_GLORYS = val.seasonal_statistics(
        Temps_GLORYS, data_obs_model_timeaxis_WASA3)

    # Extract Summer (DJF) statistics for in situ data
    summer_insitu_std = seasonal_stats_GLORYS['insitu']['DJF']['std']
    summer_insitu_correlation = seasonal_stats_GLORYS['insitu']['DJF']['correlation']
    summer_insitu_rmse = seasonal_stats_GLORYS['insitu']['DJF']['rmse']
    summer_insitu_mean_diff = seasonal_stats_GLORYS['insitu']['DJF']['mean']
    summer_insitu_min_value = seasonal_stats_GLORYS['insitu']['DJF']['min_value']
    summer_insitu_max_value = seasonal_stats_GLORYS['insitu']['DJF']['max_value']

    winter_insitu_std = seasonal_stats_GLORYS['insitu']['JJA']['std']
    winter_insitu_correlation = seasonal_stats_GLORYS['insitu']['JJA']['correlation']
    winter_insitu_rmse = seasonal_stats_GLORYS['insitu']['JJA']['rmse']
    winter_insitu_mean_diff = seasonal_stats_GLORYS['insitu']['JJA']['mean']
    winter_insitu_min_value = seasonal_stats_GLORYS['insitu']['JJA']['min_value']
    winter_insitu_max_value = seasonal_stats_GLORYS['insitu']['JJA']['max_value']

    summer_model_std = seasonal_stats_GLORYS['model']['DJF']['std']
    summer_model_correlation = seasonal_stats_GLORYS['model']['DJF']['correlation']
    summer_model_rmse = seasonal_stats_GLORYS['model']['DJF']['rmse']
    summer_model_mean_diff = seasonal_stats_GLORYS['model']['DJF']['mean']
    summer_model_min_value = seasonal_stats_GLORYS['model']['DJF']['min_value']
    summer_model_max_value = seasonal_stats_GLORYS['model']['DJF']['max_value']
    summer_model_total_bias = seasonal_stats_GLORYS['model']['DJF']['total_bias']

    winter_model_std = seasonal_stats_GLORYS['model']['JJA']['std']
    winter_model_correlation = seasonal_stats_GLORYS['model']['JJA']['correlation']
    winter_model_rmse = seasonal_stats_GLORYS['model']['JJA']['rmse']
    winter_model_mean_diff = seasonal_stats_GLORYS['model']['JJA']['mean']
    winter_model_min_value = seasonal_stats_GLORYS['model']['JJA']['min_value']
    winter_model_max_value = seasonal_stats_GLORYS['model']['JJA']['max_value']
    winter_model_total_bias = seasonal_stats_GLORYS['model']['JJA']['total_bias']

    # Call the statistics function for ERA5
    seasonal_stats_ERA5 = val.seasonal_statistics(
        data_model, data_obs_model_timeaxis_WASA3)

    summer_model_std_ERA5 = seasonal_stats_ERA5['model']['DJF']['std']
    summer_model_correlation_ERA5 = seasonal_stats_ERA5['model']['DJF']['correlation']
    summer_model_rmse_ERA5 = seasonal_stats_ERA5['model']['DJF']['rmse']
    summer_model_mean_diff_ERA5 = seasonal_stats_ERA5['model']['DJF']['mean']
    summer_model_min_value_ERA5 = seasonal_stats_ERA5['model']['DJF']['min_value']
    summer_model_max_value_ERA5 = seasonal_stats_ERA5['model']['DJF']['max_value']
    summer_model_total_bias_ERA5 = seasonal_stats_ERA5['model']['DJF']['total_bias']

    winter_model_std_ERA5 = seasonal_stats_ERA5['model']['JJA']['std']
    winter_model_correlation_ERA5 = seasonal_stats_ERA5['model']['JJA']['correlation']
    winter_model_rmse_ERA5 = seasonal_stats_ERA5['model']['JJA']['rmse']
    winter_model_mean_diff_ERA5 = seasonal_stats_ERA5['model']['JJA']['mean']
    winter_model_min_value_ERA5 = seasonal_stats_ERA5['model']['JJA']['min_value']
    winter_model_max_value_ERA5 = seasonal_stats_ERA5['model']['JJA']['max_value']
    winter_model_total_bias_ERA5 = seasonal_stats_ERA5['model']['JJA']['total_bias']

    # Call the statistics function for WASA3
    seasonal_stats_WASA3 = val.seasonal_statistics(
        data_model_WASA3, data_obs_model_timeaxis_WASA3)

    summer_model_std_WASA3 = seasonal_stats_WASA3['model']['DJF']['std']
    summer_model_correlation_WASA3 = seasonal_stats_WASA3['model']['DJF']['correlation']
    summer_model_rmse_WASA3 = seasonal_stats_WASA3['model']['DJF']['rmse']
    summer_model_mean_diff_WASA3 = seasonal_stats_WASA3['model']['DJF']['mean']
    summer_model_min_value_WASA3 = seasonal_stats_WASA3['model']['DJF']['min_value']
    summer_model_max_value_WASA3 = seasonal_stats_WASA3['model']['DJF']['max_value']
    summer_model_total_bias_WASA3 = seasonal_stats_WASA3['model']['DJF']['total_bias']

    winter_model_std_WASA3 = seasonal_stats_WASA3['model']['JJA']['std']
    winter_model_correlation_WASA3 = seasonal_stats_WASA3['model']['JJA']['correlation']
    winter_model_rmse_WASA3 = seasonal_stats_WASA3['model']['JJA']['rmse']
    winter_model_mean_diff_WASA3 = seasonal_stats_WASA3['model']['JJA']['mean']
    winter_model_min_value_WASA3 = seasonal_stats_WASA3['model']['JJA']['min_value']
    winter_model_max_value_WASA3 = seasonal_stats_WASA3['model']['JJA']['max_value']
    winter_model_total_bias_WASA3 = seasonal_stats_WASA3['model']['JJA']['total_bias']

    # Create a pandas DataFrame for plotting
    table_data = {
        r'$\mathbf{Statistic}$': [
            'In situ depth',  # 'Model depth (h)',
            'Correlation',
            'RMSE ', 'Total bias', 'Std.',
            'Mean', 'Min',
            'Max',
            r'$\mathbf{Summer\ Correlation}$',
            r'$\mathbf{Summer\ RMSE}$',
            r'$\mathbf{Summer\ Bias}$',
            r'$\mathbf{Summer\ Std.}$',
            r'$\mathbf{Summer\ Mean}$',
            r'$\mathbf{Winter\ Correlation}$',
            r'$\mathbf{Winter\ RMSE}$',
            r'$\mathbf{Winter\ Bias}$',
            r'$\mathbf{Winter\ Std.}$',
            r'$\mathbf{Winter\ Mean}$'
        ],

        r'$\mathbf{In\ situ}$': [
            f'{ds_ERA5.depth.values} m',  # '',
            '',
            '',
            '',
            f'{round(Decimal(str(ds_ERA5.insitu_std)), 3)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_mean_difference)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_min_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_max_value)), 2)} °C',

            '',
            '',
            '',
            f'{round(Decimal(str(summer_insitu_std)), 2)} °C',
            f'{round(Decimal(str(summer_insitu_mean_diff)), 2)} °C',

            '',
            '',
            '',
            f'{round(Decimal(str(winter_insitu_std)), 2)} °C',
            f'{round(Decimal(str(winter_insitu_mean_diff)), 2)} °C'
        ],

        r'$\mathbf{ERA5\ Run}$': [
            '',  # f'{-round(ds_ERA5.h)} m',
            f'{round(ds_ERA5.model_correlation_coefficient, 3)} ',
            f'{round(Decimal(str(ds_ERA5.model_rmse)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_total_bias)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_std)), 3)} °C',
            f'{round(Decimal(str(ds_ERA5.model_mean_difference)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_min_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_max_value)), 2)} °C',
            f'{round(Decimal(str(summer_model_correlation_ERA5)), 2)} ',
            f'{round(Decimal(str(summer_model_rmse_ERA5)), 2)} °C',
            f'{round(Decimal(str(summer_model_total_bias_ERA5)), 2)} °C',
            f'{round(Decimal(str(summer_model_std_ERA5)), 2)} °C',
            f'{round(Decimal(str(summer_model_mean_diff_ERA5)), 2)} °C',
            f'{round(Decimal(str(winter_model_correlation_ERA5)), 2)} ',
            f'{round(Decimal(str(winter_model_rmse_ERA5)), 2)} °C',
            f'{round(Decimal(str(winter_model_total_bias_ERA5)), 2)} °C',
            f'{round(Decimal(str(winter_model_std_ERA5)), 2)} °C',
            f'{round(Decimal(str(winter_model_mean_diff_ERA5)), 2)} °C'
        ],

        r'$\mathbf{WASA3\ Run}$': [
            '',  # f'{-round(ds_WASA3.h)} m',
            f'{round(ds_WASA3.model_correlation_coefficient, 3)} ',
            f'{round(Decimal(str(ds_WASA3.model_rmse)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_total_bias)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_std)), 3)} °C',
            f'{round(Decimal(str(ds_WASA3.model_mean_difference)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_min_value)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_max_value)), 2)} °C',
            f'{round(Decimal(str(summer_model_correlation_WASA3)), 2)} ',
            f'{round(Decimal(str(summer_model_rmse_WASA3)), 2)} °C',
            f'{round(Decimal(str(summer_model_total_bias_WASA3)), 2)} °C',
            f'{round(Decimal(str(summer_model_std_WASA3)), 2)} °C',
            f'{round(Decimal(str(summer_model_mean_diff_WASA3)), 2)} °C',
            f'{round(Decimal(str(winter_model_correlation_WASA3)), 2)} ',
            f'{round(Decimal(str(winter_model_rmse_WASA3)), 2)} °C',
            f'{round(Decimal(str(winter_model_total_bias_WASA3)), 2)} °C',
            f'{round(Decimal(str(winter_model_std_WASA3)), 2)} °C',
            f'{round(Decimal(str(winter_model_mean_diff_WASA3)), 2)} °C'
        ],

        r'$\mathbf{GLORYS}$': [
            '',  # f'{-round(ds_ERA5.h)} m',
            f'{round(model_correlation_G, 3)} ',
            f'{round(Decimal(str(model_rmse_G)), 2)} °C',
            f'{round(Decimal(str(model_total_bias_G)), 2)} °C',
            f'{round(Decimal(str( model_std_G)), 3)} °C',
            f'{round(Decimal(str(model_mean_diff_G)), 2)} °C',
            f'{round(Decimal(str(model_min_value_G)), 2)} °C',
            f'{round(Decimal(str(model_max_value_G)), 2)} °C',
            f'{round(Decimal(str(summer_model_correlation)), 2)} ',
            f'{round(Decimal(str(summer_model_rmse)), 2)} °C',
            f'{round(Decimal(str(summer_model_total_bias)), 2)} °C',
            f'{round(Decimal(str(summer_model_std)), 2)} °C',
            f'{round(Decimal(str(summer_model_mean_diff)), 2)} °C',
            f'{round(Decimal(str(winter_model_correlation)), 2)} ',
            f'{round(Decimal(str(winter_model_rmse)), 2)} °C',
            f'{round(Decimal(str(winter_model_total_bias)), 2)} °C',
            f'{round(Decimal(str(winter_model_std)), 2)} °C',
            f'{round(Decimal(str(winter_model_mean_diff)), 2)} °C'
        ]
    }

    table_df = pd.DataFrame(table_data)

    # DataFrame for Excel with both ERA5 and WASA3
    table_data_excel = {
        'Statistic': [
            'In situ depth', 'Model depth (h)', 'Correlation',
            'Std. model', 'Std. observations', 'RMSE ', 'Total bias',
            'Model mean', 'Observations mean', 'Model min', 'Observations min',
            'Model max', 'Observations max', 'JFM Insitu mean',
            'AMJ Insitu mean', 'JAS Insitu mean', 'OND Insitu mean',
            'JFM model mean', 'AMJ model mean', 'JAS model mean',
            'OND model mean'
        ],

        'Insitu': [
            f'{ds_ERA5.depth.values} m', f'{-round(ds_ERA5.h)} m',
            f'{round(ds_ERA5.model_correlation_coefficient, 3)} ',
            f'{round(Decimal(str(ds_ERA5.model_std)), 3)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_std)), 3)} °C',
            f'{round(Decimal(str(ds_ERA5.model_rmse)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_total_bias)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_mean_difference)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_mean_difference)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_min_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_min_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_max_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_max_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_seasonal_mean_JFM)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_seasonal_mean_AMJ)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_seasonal_mean_JAS)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_seasonal_mean_OND)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_seasonal_mean_JFM)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_seasonal_mean_AMJ)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_seasonal_mean_JAS)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_seasonal_mean_OND)), 2)} °C'
        ],

        'WASA3-Run': [
            f'{ds_WASA3.depth.values} m', f'{-round(ds_WASA3.h)} m',
            f'{round(ds_WASA3.model_correlation_coefficient, 3)} ',
            f'{round(Decimal(str(ds_WASA3.model_std)), 3)} °C',
            f'{round(Decimal(str(ds_WASA3.insitu_std)), 3)} °C',
            f'{round(Decimal(str(ds_WASA3.model_rmse)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_total_bias)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_mean_difference)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.insitu_mean_difference)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_min_value)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.insitu_min_value)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_max_value)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.insitu_max_value)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.insitu_seasonal_mean_JFM)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.insitu_seasonal_mean_AMJ)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.insitu_seasonal_mean_JAS)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.insitu_seasonal_mean_OND)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_seasonal_mean_JFM)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_seasonal_mean_AMJ)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_seasonal_mean_JAS)), 2)} °C',
            f'{round(Decimal(str(ds_WASA3.model_seasonal_mean_OND)), 2)} °C'
        ],

        'ERA5-Run': [
            f'{ds_ERA5.depth.values} m', f'{-round(ds_ERA5.h)} m',
            f'{round(ds_ERA5.model_correlation_coefficient, 3)} ',
            f'{round(Decimal(str(ds_ERA5.model_std)), 3)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_std)), 3)} °C',
            f'{round(Decimal(str(ds_ERA5.model_rmse)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_total_bias)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_mean_difference)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_mean_difference)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_min_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_min_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_max_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_max_value)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_seasonal_mean_JFM)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_seasonal_mean_AMJ)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_seasonal_mean_JAS)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.insitu_seasonal_mean_OND)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_seasonal_mean_JFM)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_seasonal_mean_AMJ)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_seasonal_mean_JAS)), 2)} °C',
            f'{round(Decimal(str(ds_ERA5.model_seasonal_mean_OND)), 2)} °C'
        ]
    }

    table_df_excel = pd.DataFrame(table_data_excel)

    # Update the DataFrame index to start at 1
    table_df.index = table_df.index + 1

    # Increase figure size for better readability
    fig, ax = plt.subplots(figsize=(8, 7))

    # Plot the DataFrame as a table with left-aligned text and wider columns
    ax.axis('off')
    tbl = table(ax, table_df, loc='center', cellLoc='left', colWidths=[
                0.3, 0.2, 0.2, 0.2, 0.2], rowLoc='center', fontsize=12)

    # Adjust the row heights
    for i, key in enumerate(tbl.get_celld().keys()):
        cell = tbl[key]
        cell.set_fontsize(12)
        cell.set_height(0.07)  # Adjust the height as needed

    # Save the figure as a PNG file
    # savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/'
    # filename = 'filename'  # Replace with your actual filename
    savename_tbl = 'table.png'  # Replace with your desired table image name
    plt.savefig(savepath + filename[:-3] + " Model (h) ="+f'{round(ds_ERA5.h)} m'+" stats " +
                savename_tbl, bbox_inches='tight', dpi=300)  # Increased DPI for better quality
    plt.show()

    # Save the DataFrame as an Excel file with both ERA5 and WASA3 data
    excel_savepath = savepath + filename + '_stats_table.xlsx'
    table_df_excel.to_excel(excel_savepath, index=False)


# %%Ploting the sea state map with stations based on the CROCO MODEL domain boundaries


# Initialize lists to store latitude and longitude values
all_lat_insitu = []
all_lon_insitu = []
all_lat_model = []
all_lon_model = []

# Combine filenames into a single list with corresponding paths
files_and_paths = {
    "AJSMIT_UTR": [
        "Betty's Bay_DAFF.nc",
        "Bordjies Deep_DAFF.nc",
        "Bordjies_DAFF.nc",
        "Doringbaai_DAFF.nc",
        "Fish Hoek_SAWS.nc",
        "Gordons Bay_SAWS.nc",
        "Kalk Bay_SAWS.nc",
        "Kommetjie_SAWS.nc",
        "Lamberts Bay_SAWS.nc",
        "Muizenberg_SAWS.nc",
        "Oudekraal_DAFF.nc",
        "Saldanha Bay_SAWS.nc",
        "Sea Point_SAWS.nc",
        "St Helena Bay_SAWS.nc",
        "Yzerfontein_SAWS.nc"
    ],
    "GPITCHER_SHB": [
        "WQM_20m_TS.nc",
        "WQM_70m_TS.nc"
    ],
    "ATAP": [
        "FalseBay_FB001.nc",
        "CapePoint_CP001.nc"
    ]
}

# Paths for each group
paths = {
    "AJSMIT_UTR": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/AJSMIT_UTR/',
    "GPITCHER_SHB": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/GPITCHER_SHB/',
    "ATAP": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/ATAP/'
}

# Single loop to process all files and accumulate latitude and longitude data
for group, filenames in files_and_paths.items():
    for filename in filenames:
        ds_ERA5 = xr.open_dataset(paths[group] + "Validated_" + filename)

        lat_insitu = ds_ERA5.latitude.values
        lon_insitu = ds_ERA5.longitude.values
        lat_model = ds_ERA5.latitude.values
        lon_model = ds_ERA5.longitude.values

        all_lat_insitu.extend(lat_insitu)
        all_lon_insitu.extend(lon_insitu)
        all_lat_model.extend(lat_model)
        all_lon_model.extend(lon_model)

# Function to get CROCO boundary


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values
        depth = ds_ERA5.h.values  # depth is stored in 'h' variable
        # Compute mean temperature at surface level
        temp = ds_ERA5.temp.mean(dim='time').isel(s_rho=29)

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, temp


# Assuming fname is defined and points to the CROCO file
fname = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/output/croco_avg_Y2009M07.nc"

# Get the boundaries from the CROCO grid
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, temp = get_croco_boundary(
    fname)

# Main map setup
fig = plt.figure(figsize=(12, 10), facecolor='white')
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot mean temperature
temp_levels = np.linspace(8, 30, 20)
cf = ax.contourf(lon_rho, lat_rho, temp, levels=temp_levels,
                 cmap='turbo', extend='both')  # turbo is an alternative cmap
cb = plt.colorbar(cf, ax=ax, orientation='vertical',
                  label='Mean Temperature (°C)')
cb.set_label('Surface Temperature (°C)', fontsize=14)
cb.ax.tick_params(labelsize=14)

# Plot depth contours
contour_levels = [100, 200, 500]
cs = ax.contour(lon_rho, lat_rho, depth, levels=contour_levels,
                colors='k', linestyles='--')
ax.clabel(cs, inline=True, fontsize=8, fmt='%d m')

# Overlay land features to conceal overlaps
ax.add_feature(cfeature.LAND, zorder=1)
ax.coastlines(zorder=2)

# Plot in situ and model data points
ax.plot(all_lon_insitu, all_lat_insitu, 'or', label='In situ', zorder=3,
        markeredgecolor='black', markeredgewidth=1, markersize=11)

# Plot the boundary as separate lines
ax.plot(lon_west, lat_west, 'b-', linewidth=2,
        label='CROCO Boundary', zorder=3)
ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=3)
ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=3)

# Add custom legend entry
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles, labels, loc='upper right',
           bbox_to_anchor=(1, 1), fontsize=9)

plt.title('Position of all sampling stations and croco ERA5 run model \n 1 July 2009', fontsize=16)

# Set the extent using the boundary coordinates
all_lon = np.hstack((lon_west, lon_south, lon_north))
all_lat = np.hstack((lat_west, lat_south, lat_north))
ax.set_extent([np.min(all_lon), np.max(all_lon),
              np.min(all_lat), np.max(all_lat)])

gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0,
                  color='black', draw_labels=True)
gl.left_labels = True
gl.top_labels = False
gl.right_labels = False
gl.xlines = True
gl.ylines = True

# Subplot for Africa map
ax_inset = fig.add_axes([0.245, 0.14, 0.2, 0.2], projection=ccrs.PlateCarree())
ax_inset.set_extent([-20, 60, -40, 40], crs=ccrs.PlateCarree())
ax_inset.add_feature(cfeature.LAND, edgecolor='black')
ax_inset.add_feature(cfeature.COASTLINE)
ax_inset.add_feature(cfeature.BORDERS, linestyle='-')

# Highlight the study domain with a red box
study_lon_min, study_lon_max = np.min(all_lon), np.max(all_lon)
study_lat_min, study_lat_max = np.min(all_lat), np.max(all_lat)
rect = plt.Rectangle((study_lon_min, study_lat_min),
                     study_lon_max - study_lon_min,
                     study_lat_max - study_lat_min,
                     linewidth=2, edgecolor='red', facecolor='none', transform=ccrs.PlateCarree())
ax_inset.add_patch(rect)

# Remove labels and ticks from the inset map
ax_inset.set_xticks([])
ax_inset.set_yticks([])

# Save and show the plot
plt.savefig('Map_ERA5_all_stations_and_temp_1Jul2009.png')
plt.show()

# %% Different station datasets, different polygons


# Combine filenames into a single list with corresponding paths
files_and_paths = {
    "AJSMIT_UTR": [
        "Betty's Bay_DAFF.nc",
        "Bordjies Deep_DAFF.nc",
        "Bordjies_DAFF.nc",
        "Doringbaai_DAFF.nc",
        "Fish Hoek_SAWS.nc",
        "Gordons Bay_SAWS.nc",
        "Kalk Bay_SAWS.nc",
        "Kommetjie_SAWS.nc",
        "Lamberts Bay_SAWS.nc",
        "Muizenberg_SAWS.nc",
        "Oudekraal_DAFF.nc",
        "Saldanha Bay_SAWS.nc",
        "Sea Point_SAWS.nc",
        "St Helena Bay_SAWS.nc",
        "Yzerfontein_SAWS.nc"
    ],
    "GPITCHER_SHB": [
        "WQM_20m_TS.nc",
        "WQM_70m_TS.nc"
    ],
    "ATAP": [
        "FalseBay_FB001.nc",
        "CapePoint_CP001.nc"
    ],
    "Wirewalker_Mooring": [
        "wirewalker_mooring_1.nc"
    ],
    "ADCP_Mooring": [
        "adcp_mooring_1.nc"
    ]
}

# Paths for each group
paths = {
    "AJSMIT_UTR": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/AJSMIT_UTR/',
    "GPITCHER_SHB": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/GPITCHER_SHB/',
    "ATAP": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/ATAP/',
    "Wirewalker_Mooring": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/Wirewalker/',
    "ADCP_Mooring": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/adcp_mooring/'
}

savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/'

# Initialize lists to store latitude and longitude values
all_lat_insitu = []
all_lon_insitu = []

# Separate lists for each group
gpitcher_lat_insitu = []
gpitcher_lon_insitu = []
atap_lat_insitu = []
atap_lon_insitu = []
wirewalker_lat_insitu = []
wirewalker_lon_insitu = []
adcp_lat_insitu = []
adcp_lon_insitu = []

# Single loop to process all files and accumulate latitude and longitude data
for group, filenames in files_and_paths.items():
    for filename in filenames:
        ds_ERA5 = xr.open_dataset(paths[group] + "Validated_" + filename)

        lat_insitu = ds_ERA5.latitude.values
        lon_insitu = ds_ERA5.longitude.values

        if group == "GPITCHER_SHB":
            # Store GPITCHER_SHB data separately
            gpitcher_lat_insitu.extend(lat_insitu)
            gpitcher_lon_insitu.extend(lon_insitu)
        elif group == "ATAP":
            # Store ATAP data separately
            atap_lat_insitu.extend(lat_insitu)
            atap_lon_insitu.extend(lon_insitu)
        elif group == "Wirewalker_Mooring":
            # Store Wirewalker_Mooring data separately
            wirewalker_lat_insitu.extend(lat_insitu)
            wirewalker_lon_insitu.extend(lon_insitu)
        elif group == "ADCP_Mooring":
            # Store ADCP_Mooring data separately
            adcp_lat_insitu.extend(lat_insitu)
            adcp_lon_insitu.extend(lon_insitu)
        else:
            # Store other data in the main lists
            all_lat_insitu.extend(lat_insitu)
            all_lon_insitu.extend(lon_insitu)

# Function to get CROCO boundary


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values
        depth = ds_ERA5.h.values  # depth is stored in 'h' variable
        # Compute mean temperature at surface level
        temp = ds_ERA5.temp.mean(dim='time').isel(s_rho=29)

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, temp


fname = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/output/croco_avg_Y2009M07.nc"

# Get the boundaries from the CROCO grid
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, temp = get_croco_boundary(
    fname)

# Plot the data on the map
fig = plt.figure(figsize=(12, 10), facecolor='white')
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot mean temperature
temp_levels = np.linspace(8, 31, 20)
cf = ax.contourf(lon_rho, lat_rho, temp, levels=temp_levels,
                 cmap='turbo', extend='both')
cb = plt.colorbar(cf, ax=ax, orientation='vertical',
                  label='Mean Temperature (°C)')
cb.set_label('Surface Temperature (°C)', fontsize=14)
cb.ax.tick_params(labelsize=14)

# Plot depth contours
cs = ax.contour(lon_rho, lat_rho, depth, levels=contour_levels,
                colors='k', linestyles='--')
ax.clabel(cs, inline=True, fontsize=8, fmt='%d m')

# Overlay land features to conceal overlaps
ax.add_feature(cfeature.LAND, zorder=1)
ax.coastlines(zorder=2)

# Plot in situ data points (red circles)
ax.plot(all_lon_insitu, all_lat_insitu, 'or', label='In situ (AJSMIT_UTR)', zorder=3,
        markeredgecolor='black', markeredgewidth=1, markersize=11)

# Plot GPITCHER_SHB stations (orange stars)
ax.plot(gpitcher_lon_insitu, gpitcher_lat_insitu, '*', color='orange', label='GPITCHER_SHB Stations', zorder=3,
        markeredgecolor='black', markeredgewidth=1, markersize=20)

# Plot ATAP stations (yellow boxes)
ax.plot(atap_lon_insitu, atap_lat_insitu, 's', color='yellow', label='ATAP Stations', zorder=3,
        markeredgecolor='black', markeredgewidth=1, markersize=9)

# Plot Wirewalker_Mooring stations (black diamonds)
ax.plot(wirewalker_lon_insitu, wirewalker_lat_insitu, '^', color='white', label='Wirewalker Mooring', zorder=3,
        markeredgecolor='black', markeredgewidth=1, markersize=11)

# Plot ADCP_Mooring stations (pink triangles)
ax.plot(adcp_lon_insitu, adcp_lat_insitu, '^', color='white', label='ADCP Mooring', zorder=3,
        markeredgecolor='black', markeredgewidth=1, markersize=11)

# Plot the boundary as separate lines
ax.plot(lon_west, lat_west, 'b-', linewidth=2,
        label='CROCO Boundary', zorder=3)
ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=3)
ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=3)

# Add custom legend entry
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles, labels, loc='upper right',
           bbox_to_anchor=(1, 1), fontsize=9)

plt.title('Position of all sampling stations and croco ERA5 run model \n 1 July 2009', fontsize=16)

# Set the extent using the boundary coordinates
all_lon = np.hstack((lon_west, lon_south, lon_north))
all_lat = np.hstack((lat_west, lat_south, lat_north))
ax.set_extent([np.min(all_lon), np.max(all_lon),
              np.min(all_lat), np.max(all_lat)])

gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0,
                  color='black', draw_labels=True)
gl.left_labels = True
gl.top_labels = False
gl.right_labels = False
gl.xlines = True
gl.ylines = True

# Subplot for Africa map
ax_inset = fig.add_axes([0.245, 0.14, 0.2, 0.2], projection=ccrs.PlateCarree())
ax_inset.set_extent([-20, 60, -40, 40], crs=ccrs.PlateCarree())
ax_inset.add_feature(cfeature.LAND, edgecolor='black')
ax_inset.add_feature(cfeature.COASTLINE)
ax_inset.add_feature(cfeature.BORDERS, linestyle='-')

# Highlight the study domain with a red box
rect = plt.Rectangle((study_lon_min, study_lat_min),
                     study_lon_max - study_lon_min,
                     study_lat_max - study_lat_min,
                     linewidth=2, edgecolor='red', facecolor='none', transform=ccrs.PlateCarree())
ax_inset.add_patch(rect)

# Remove labels and ticks from the inset map
ax_inset.set_xticks([])
ax_inset.set_yticks([])

# Save and show the plot
plt.savefig(savepath + 'Map_ERA5_all_stations_and_temp_1Jul2009.png')
plt.show()


# %% Mean Sea State over the 5 years.


all_lat_insitu = []
all_lon_insitu = []
all_lat_model = []
all_lon_model = []

for filename in filenames:
    ds_ERA5 = xr.open_dataset(ERA5_VAL_path + "Validated_" + filename)

    lat_insitu = ds_ERA5.latitude.values
    lon_insitu = ds_ERA5.longitude.values
    lat_model = ds_ERA5.latitude.values
    lon_model = ds_ERA5.longitude.values

    all_lat_insitu.extend(lat_insitu)
    all_lon_insitu.extend(lon_insitu)
    all_lat_model.extend(lat_model)
    all_lon_model.extend(lon_model)


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values
        depth = ds_ERA5.h.values  # Assuming depth is stored in 'h' variable
        # Compute mean temperature at surface level
        temp = ds_ERA5.temp.mean(dim='time').isel(s_rho=-1)

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, temp


# Collect all relevant .nc files using a wildcard
filenames = glob.glob(
    '/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y20*.nc')
# fname               = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y2013M01.nc"  #'path_to_your_croco_file.nc'

# Compute the mean sea state across all files
ds_croco = xr.open_mfdataset(filenames, combine='by_coords')
mean_temp = ds_croco.temp.mean(dim='time').isel(
    s_rho=0)  # Surface temperature as an example

# Get the boundaries from one of the CROCO grid files
fname = filenames[0]  # Using the first file for boundary extraction
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, _ = get_croco_boundary(
    fname)

# Main map setup
fig = plt.figure(figsize=(12, 10), facecolor='white')
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot mean temperature
temp_levels = np.linspace(mean_temp.min(), mean_temp.max()+5, 20)
cf = ax.contourf(lon_rho, lat_rho, mean_temp, levels=temp_levels,
                 cmap='turbo', extend='both')  # 'rainbow' is an alternative cmap i'll try
cb = plt.colorbar(cf, ax=ax, orientation='vertical',
                  label='Mean Sea State Temperature (°C)')

# Plot depth contours
contour_levels = [100, 200, 500]
cs = ax.contour(lon_rho, lat_rho, depth, levels=contour_levels,
                colors='k', linestyles='--')
ax.clabel(cs, inline=True, fontsize=8, fmt='%d m')

# Overlay land features to conceal overlaps
ax.add_feature(cfeature.LAND, zorder=1)
ax.coastlines(zorder=2)

# Plot in situ and model data points
ax.plot(all_lon_insitu, all_lat_insitu, 'or', label='In situ', zorder=3,
        markeredgecolor='black', markeredgewidth=1, markersize=11)

# Plot the boundary as separate lines
ax.plot(lon_west, lat_west, 'b-', linewidth=2,
        label='CROCO Boundary', zorder=3)
ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=3)
ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=3)

# Add custom legend entry
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles, labels, loc='upper right',
           bbox_to_anchor=(1, 1), fontsize=9)

plt.title('Mean Sea State over the Entire Record')

# Set the extent using the boundary coordinates
all_lon = np.hstack((lon_west, lon_south, lon_north))
all_lat = np.hstack((lat_west, lat_south, lat_north))
ax.set_extent([np.min(all_lon), np.max(all_lon),
              np.min(all_lat), np.max(all_lat)])

gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0,
                  color='black', draw_labels=True)
gl.left_labels = True
gl.top_labels = False
gl.right_labels = False
gl.xlines = True
gl.ylines = True

# Subplot for Africa map
ax_inset = fig.add_axes([0.245, 0.14, 0.2, 0.2], projection=ccrs.PlateCarree())
ax_inset.set_extent([-20, 60, -40, 40], crs=ccrs.PlateCarree())
ax_inset.add_feature(cfeature.LAND, edgecolor='black')
ax_inset.add_feature(cfeature.COASTLINE)
ax_inset.add_feature(cfeature.BORDERS, linestyle='-')

# Highlight the study domain with a red box
study_lon_min, study_lon_max = np.min(all_lon), np.max(all_lon)
study_lat_min, study_lat_max = np.min(all_lat), np.max(all_lat)
rect = plt.Rectangle((study_lon_min, study_lat_min),
                     study_lon_max - study_lon_min,
                     study_lat_max - study_lat_min,
                     linewidth=2, edgecolor='red', facecolor='none', transform=ccrs.PlateCarree())
ax_inset.add_patch(rect)

# Remove labels and ticks from the inset map
ax_inset.set_xticks([])
ax_inset.set_yticks([])

plt.savefig(savepath + 'Mean_Sea_State_2009-2013.png')
plt.show()

# %% SEASONS


all_lat_insitu = []
all_lon_insitu = []
all_lat_model = []
all_lon_model = []

for filename in filenames:
    ds_ERA5 = xr.open_dataset(ERA5_VAL_path + "Validated_" + filename)

    lat_insitu = ds_ERA5.latitude.values
    lon_insitu = ds_ERA5.longitude.values
    lat_model = ds_ERA5.latitude.values
    lon_model = ds_ERA5.longitude.values

    all_lat_insitu.extend(lat_insitu)
    all_lon_insitu.extend(lon_insitu)
    all_lat_model.extend(lat_model)
    all_lon_model.extend(lon_model)


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values
        depth = ds_ERA5.h.values  # Assuming depth is stored in 'h' variable
        # Compute mean temperature at surface level
        temp = ds_ERA5.temp.mean(dim='time').isel(s_rho=0)

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, temp


# Collect all relevant .nc files using a wildcard
filenames = glob.glob(
    '/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y20*.nc')
# fname               = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y2013M01.nc"  #'path_to_your_croco_file.nc'

# Compute the mean sea state across all files
ds_croco = xr.open_mfdataset(filenames, combine='by_coords')

# Define the seasons
seasons = {
    'DJF': ds_croco.sel(time=ds_croco['time.season'] == 'DJF').temp.mean(dim='time').isel(s_rho=0),
    'MAM': ds_croco.sel(time=ds_croco['time.season'] == 'MAM').temp.mean(dim='time').isel(s_rho=0),
    'JJA': ds_croco.sel(time=ds_croco['time.season'] == 'JJA').temp.mean(dim='time').isel(s_rho=0),
    'SON': ds_croco.sel(time=ds_croco['time.season'] == 'SON').temp.mean(dim='time').isel(s_rho=0)
}

# Get the boundaries from one of the CROCO grid files
fname = filenames[0]  # Using the first file for boundary extraction
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, _ = get_croco_boundary(
    fname)

# Main map setup with 4 subplots
fig, axs = plt.subplots(2, 2, figsize=(16, 12), subplot_kw={
                        'projection': ccrs.PlateCarree()})
fig.suptitle('Seasonal Mean Sea State Temperature (°C)', fontsize=16)

season_titles = ['DJF (Winter)', 'MAM (Spring)',
                 'JJA (Summer)', 'SON (Autumn)']

for ax, (season, mean_temp) in zip(axs.flatten(), seasons.items()):
    # Plot mean temperature
    temp_levels = np.linspace(mean_temp.min(), mean_temp.max()+5, 20)
    cf = ax.contourf(lon_rho, lat_rho, mean_temp,
                     levels=temp_levels, cmap='rainbow', extend='both')

    # Plot depth contours
    contour_levels = [100, 200, 500]
    cs = ax.contour(lon_rho, lat_rho, depth,
                    levels=contour_levels, colors='k', linestyles='--')
    ax.clabel(cs, inline=True, fontsize=8, fmt='%d m')

    # Overlay land features to conceal overlaps
    ax.add_feature(cfeature.LAND, zorder=1)
    ax.coastlines(zorder=2)

    # Plot the boundary as separate lines
    ax.plot(lon_west, lat_west, 'b-', linewidth=2,
            label='CROCO Boundary', zorder=3)
    ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=3)
    ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=3)

    # Set the extent using the boundary coordinates
    all_lon = np.hstack((lon_west, lon_south, lon_north))
    all_lat = np.hstack((lat_west, lat_south, lat_north))
    ax.set_extent([np.min(all_lon), np.max(all_lon),
                  np.min(all_lat), np.max(all_lat)])

    ax.set_title(season_titles.pop(0))

    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0,
                      color='black', draw_labels=True)
    gl.left_labels = True
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = True
    gl.ylines = True

# Add a colorbar for the entire figure
cbar = fig.colorbar(cf, ax=axs, orientation='vertical',
                    fraction=0.02, pad=0.04)
cbar.set_label('Mean Sea State Temperature (°C)')

plt.savefig(savepath + 'Seasonal_Mean_Sea_State_2009-2013.png')
plt.show()


# %% GLORYS SST MAP with insitu observations


# Load GLORYS data
glorys_fname = "/mnt/e/GLORYS_SA_EEZ_ZOOM/GLORYS_Y2009M1.nc"
GLORYS_ds = xr.open_dataset(glorys_fname)

# Extract surface temperature (assuming first depth level is surface)
surface_temp_glorys = GLORYS_ds['temp'].isel(depth=20).mean(dim='time')

# Define CROCO file path
croco_fname = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y2013M01.nc"

# Function to get CROCO grid boundary and variables


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values
        depth = ds_ERA5.h.values  # depth is stored in 'h' variable
        # Compute mean temperature at surface level
        temp = ds_ERA5.temp.mean(dim='time').isel(s_rho=29)

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, temp


# Get the boundaries from the CROCO grid
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho, depth, temp_croco = get_croco_boundary(
    croco_fname)

# Regrid the GLORYS surface temperature to the CROCO grid
glorys_lon, glorys_lat = np.meshgrid(GLORYS_ds.lonT, GLORYS_ds.latT)
temp_glorys_on_croco_grid = griddata(
    (glorys_lon.ravel(), glorys_lat.ravel()),
    surface_temp_glorys.values.ravel(),
    (lon_rho, lat_rho),
    method='linear'
)

# Initialize lists for in situ and model data points
all_lat_insitu = []
all_lon_insitu = []
all_lat_model = []
all_lon_model = []

# List of filenames and ERA5 validation path
filenames = [
    "Betty's Bay_DAFF.nc",
    "Bordjies Deep_DAFF.nc",
    "Bordjies_DAFF.nc",
    "Doringbaai_DAFF.nc",
    "Fish Hoek_SAWS.nc",
    "Gordons Bay_SAWS.nc",
    "Kalk Bay_SAWS.nc",
    "Kommetjie_SAWS.nc",
    "Lamberts Bay_SAWS.nc",
    "Muizenberg_SAWS.nc",
    "Oudekraal_DAFF.nc",
    "Saldanha Bay_SAWS.nc",
    "Sea Point_SAWS.nc",
    "St Helena Bay_SAWS.nc",
    "Yzerfontein_SAWS.nc"
]
ERA5_VAL_path = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/AJSMIT_UTR/'

for filename in filenames:
    ds_ERA5 = xr.open_dataset(ERA5_VAL_path + "Validated_" + filename)

    lat_insitu = ds_ERA5.latitude.values
    lon_insitu = ds_ERA5.longitude.values
    lat_model = ds_ERA5.latitude.values
    lon_model = ds_ERA5.longitude.values

    all_lat_insitu.extend(lat_insitu)
    all_lon_insitu.extend(lon_insitu)
    all_lat_model.extend(lat_model)
    all_lon_model.extend(lon_model)

# Main map setup
fig = plt.figure(figsize=(12, 10), facecolor='white')
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot GLORYS surface temperature
temp_levels = np.linspace(temp_glorys_on_croco_grid.min(),
                          temp_glorys_on_croco_grid.max(), 20)
cf = ax.contourf(lon_rho, lat_rho, temp_glorys_on_croco_grid,
                 levels=temp_levels, cmap='turbo', extend='both')
cb = plt.colorbar(cf, ax=ax, orientation='vertical',
                  label='Surface Temperature (°C)')

# Plot depth contours (if needed)
contour_levels = [100, 200, 500]
cs = ax.contour(lon_rho, lat_rho, depth, levels=contour_levels,
                colors='k', linestyles='--')
ax.clabel(cs, inline=True, fontsize=8, fmt='%d m')

# Overlay land features to conceal overlaps
ax.add_feature(cfeature.LAND, zorder=1)
ax.coastlines(zorder=2)

# Plot in situ and model data points
ax.plot(all_lon_insitu, all_lat_insitu, 'or', label='In situ', zorder=3,
        markeredgecolor='black', markeredgewidth=1, markersize=11)

# Plot the CROCO boundary as separate lines
ax.plot(lon_west, lat_west, 'b-', linewidth=2,
        label='CROCO Boundary', zorder=3)
ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=3)
ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=3)

# Add custom legend entry
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles, labels, loc='upper right',
           bbox_to_anchor=(1, 1), fontsize=9)

plt.title('GLORYS Surface Temperature within CROCO Boundaries')

# Set the extent using the boundary coordinates
all_lon = np.hstack((lon_west, lon_south, lon_north))
all_lat = np.hstack((lat_west, lat_south, lat_north))

# Correctly calculate the study domain boundaries
study_lon_min = np.min(all_lon)
study_lon_max = np.max(all_lon)
study_lat_min = np.min(all_lat)
study_lat_max = np.max(all_lat)

ax.set_extent([study_lon_min, study_lon_max, study_lat_min, study_lat_max])

# Gridlines and inset map
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0,
                  color='black', draw_labels=True)
gl.left_labels = True
gl.top_labels = False
gl.right_labels = False
gl.xlines = True
gl.ylines = True

# Subplot for Africa map
ax_inset = fig.add_axes([0.245, 0.14, 0.2, 0.2], projection=ccrs.PlateCarree())
ax_inset.set_extent([-20, 60, -40, 40], crs=ccrs.PlateCarree())
ax_inset.add_feature(cfeature.LAND, edgecolor='black')
ax_inset.add_feature(cfeature.COASTLINE)
ax_inset.add_feature(cfeature.BORDERS, linestyle='-')

# Highlight the study domain with a red box
rect = plt.Rectangle((study_lon_min, study_lat_min),
                     study_lon_max - study_lon_min,
                     study_lat_max - study_lat_min,
                     linewidth=2, edgecolor='red', facecolor='none', transform=ccrs.PlateCarree())
ax_inset.add_patch(rect)

# Remove labels and ticks from the inset map
ax_inset.set_xticks([])
ax_inset.set_yticks([])

# Save and display the figure
savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/'
plt.savefig(savepath + 'Map_with_GLORYS_surface_temp_with_inset.png')
plt.show()

# %% glorys full domain map


# Load GLORYS data
glorys_fname = "/mnt/e/GLORYS_SA_EEZ_ZOOM/GLORYS_Y2009M1.nc"
GLORYS_ds = xr.open_dataset(glorys_fname)

# Extract the surface temperature (first depth level, first time step)
surface_temp = GLORYS_ds['temp'].isel(depth=0, time=0)

# Plot the surface temperature
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
contour = surface_temp.plot(
    ax=ax,
    cmap='coolwarm',
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8},
    vmin=0,
    vmax=35
)

# Add coastlines, land mask, and borders
ax.add_feature(cfeature.LAND, facecolor='sandybrown',
               zorder=10)  # Light brown land mask
ax.coastlines(resolution='10m', color='black',
              linewidth=1, zorder=11)  # Coastlines
ax.add_feature(cfeature.BORDERS, edgecolor='black',
               linewidth=0.5, zorder=12)  # Admin borders

# Set plot title and labels
plt.title('GLORYS Surface Temperature - January 1, 2009')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Save the plot to the specified path
savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/glorys_surface_temp_jan_2009.png'
plt.savefig(savepath, dpi=300, bbox_inches='tight')

plt.show()


# %% glorys full domain map with CROCO boundary


# Load GLORYS data
glorys_fname = "/mnt/e/GLORYS_SA_EEZ_ZOOM/GLORYS_Y2009M1.nc"
GLORYS_ds = xr.open_dataset(glorys_fname)

# Extract the surface temperature (first depth level, first time step)
surface_temp = GLORYS_ds['temp'].isel(depth=0, time=0)

# Function to get CROCO grid boundary and variables


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho


# Define CROCO file path and extract boundaries
croco_fname = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y2013M01.nc"
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho = get_croco_boundary(
    croco_fname)

# Plot the surface temperature
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
contour = surface_temp.plot(
    ax=ax,
    cmap='coolwarm',
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8},
    vmin=0,
    vmax=35
)

# Overlay the CROCO boundary as separate lines
ax.plot(lon_west, lat_west, 'b-', linewidth=2,
        label='CROCO Boundary', zorder=13)
ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=13)
ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=13)

# Add coastlines, land mask, and borders
ax.add_feature(cfeature.LAND, facecolor='sandybrown',
               zorder=10)  # Light brown land mask
ax.coastlines(resolution='10m', color='black',
              linewidth=1, zorder=11)  # Coastlines
ax.add_feature(cfeature.BORDERS, edgecolor='black',
               linewidth=0.5, zorder=12)  # Admin borders

# Set plot title and labels
plt.title('GLORYS Surface Temperature with CROCO Boundary - January 1, 2009')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Save the plot to the specified path
savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/glorys_surface_temp_jan_2009_with_croco_boundary.png'
plt.savefig(savepath, dpi=300, bbox_inches='tight')

plt.show()

# %% glorys full domain map with CROCO boundary and a ZOOM


# Load GLORYS data
glorys_fname = "/mnt/e/GLORYS_SA_EEZ_ZOOM/GLORYS_Y2009M7.nc"
GLORYS_ds = xr.open_dataset(glorys_fname)

# Extract the surface temperature (first depth level, first time step)
surface_temp = GLORYS_ds['temp'].isel(depth=0, time=0)

# Function to get CROCO grid boundary and variables


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho


# Define CROCO file path and extract boundaries
croco_fname = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y2013M01.nc"
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho = get_croco_boundary(
    croco_fname)

# Plot the surface temperature (Full GLORYS domain)
fig, axs = plt.subplots(1, 2, figsize=(18, 8), subplot_kw={
                        'projection': ccrs.PlateCarree()})

# Full domain plot
ax = axs[0]
contour_full = surface_temp.plot(
    ax=ax,
    cmap='coolwarm',
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8},
    vmin=0,
    vmax=35,
    levels=np.linspace(0, 40, 21)  # 20 intervals between 0 and 35 (21 levels)
)

cbar = contour_full.colorbar
cbar.set_label('Surface Temperature (°C)', fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Overlay the CROCO boundary on the full domain plot
ax.plot(lon_west, lat_west, 'b-', linewidth=2,
        label='CROCO Boundary', zorder=13)
ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=13)
ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=13)

# Add coastlines, land mask, and borders
ax.add_feature(cfeature.LAND, facecolor='sandybrown', zorder=10)
ax.coastlines(resolution='10m', color='black', linewidth=1, zorder=11)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5, zorder=12)

# Add gridlines and labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Set title, labels, and aspect ratio
ax.set_title(
    'GLORYS Surface Temperature - Full Domain \n 01 July 2009', fontsize=16)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Zoomed plot (within CROCO boundary)
ax_zoom = axs[1]

# Use correct coordinate names ('lonT' and 'latT')
croco_mask = (surface_temp['lonT'] >= lon_rho.min()) & (surface_temp['lonT'] <= lon_rho.max()) & \
             (surface_temp['latT'] >= lat_rho.min()) & (
                 surface_temp['latT'] <= lat_rho.max())

# Restrict the data to the CROCO boundary region
surface_temp_zoom = surface_temp.where(croco_mask, drop=True)

contour_zoom = surface_temp_zoom.plot(
    ax=ax_zoom,
    cmap='coolwarm',
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8, },
    vmin=0,
    vmax=35,
    levels=np.linspace(0, 40, 21)  # 20 intervals between 0 and 35 (21 levels)
)

cbar = contour_zoom.colorbar
cbar.set_label('Surface Temperature (°C)', fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Overlay the CROCO boundary on the zoomed plot
ax_zoom.plot(lon_west, lat_west, 'b-', linewidth=2,
             label='CROCO Boundary', zorder=13)
ax_zoom.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=13)
ax_zoom.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=13)

# Add coastlines, land mask, and borders
ax_zoom.add_feature(cfeature.LAND, facecolor='sandybrown', zorder=10)
ax_zoom.coastlines(resolution='10m', color='black', linewidth=1, zorder=11)
ax_zoom.add_feature(cfeature.BORDERS, edgecolor='black',
                    linewidth=0.5, zorder=12)

# Add gridlines and labels
gl_zoom = ax_zoom.gridlines(crs=ccrs.PlateCarree(
), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl_zoom.top_labels = False
gl_zoom.right_labels = False

# Zoom in on the CROCO region
ax_zoom.set_extent([lon_rho.min(), lon_rho.max(),
                   lat_rho.min(), lat_rho.max()], crs=ccrs.PlateCarree())

# Set title, labels, and aspect ratio
ax_zoom.set_title(
    'GLORYS Surface Temperature - Zoomed to CROCO Boundary \n 01 July 2009', fontsize=16)
ax_zoom.set_xlabel('Longitude')
ax_zoom.set_ylabel('Latitude')

# Adjust layout to ensure both plots are the same size
plt.tight_layout()

# Save the combined plot
savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/glorys_surface_temp_with_zoom_jan_2009.png'
plt.savefig(savepath, dpi=300, bbox_inches='tight')

plt.show()


# %%


# Load GLORYS data
glorys_fname = "/mnt/e/GLORYS_SA_EEZ_ZOOM/GLORYS_Y2009M1.nc"
GLORYS_ds = xr.open_dataset(glorys_fname)

# Extract the surface temperature (first depth level, first time step)
surface_temp = GLORYS_ds['temp'].isel(depth=0, time=0)

# Function to get CROCO grid boundary and variables


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho


# Define CROCO file path and extract boundaries
croco_fname = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y2013M01.nc"
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho = get_croco_boundary(
    croco_fname)

# Plot the surface temperature (Full GLORYS domain)
fig, axs = plt.subplots(1, 2, figsize=(18, 8), subplot_kw={
                        'projection': ccrs.PlateCarree()})

# Full domain plot
ax = axs[0]
contour_full = surface_temp.plot(
    ax=ax,
    cmap='turbo',
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8},
    vmin=0,
    vmax=35,
    levels=np.linspace(0, 40, 21)  # 20 intervals between 0 and 40 (21 levels)
)

cbar = contour_full.colorbar
cbar.set_label('Surface Temperature (°C)', fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Overlay the CROCO boundary on the full domain plot
ax.plot(lon_west, lat_west, 'b-', linewidth=2,
        label='CROCO Boundary', zorder=13)
ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=13)
ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=13)

# Add coastlines, land mask, and borders using add_feature
ax.add_feature(cfeature.LAND, zorder=10)
ax.coastlines(resolution='10m', color='black', linewidth=1, zorder=11)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5, zorder=12)

# Add gridlines and labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Set title, labels, and aspect ratio
ax.set_title(
    'GLORYS Surface Temperature - Full Domain \n 01 January 2009', fontsize=16)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Zoomed plot (within CROCO boundary)
ax_zoom = axs[1]

# Use correct coordinate names ('lonT' and 'latT')
croco_mask = (surface_temp['lonT'] >= lon_rho.min()) & (surface_temp['lonT'] <= lon_rho.max()) & \
             (surface_temp['latT'] >= lat_rho.min()) & (
                 surface_temp['latT'] <= lat_rho.max())

# Restrict the data to the CROCO boundary region
surface_temp_zoom = surface_temp.where(croco_mask, drop=True)

contour_zoom = surface_temp_zoom.plot(
    ax=ax_zoom,
    cmap='turbo',  # coolwarm
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8},
    vmin=0,
    vmax=35,
    levels=np.linspace(0, 40, 21)  # 20 intervals between 0 and 40 (21 levels)
)

cbar = contour_zoom.colorbar
cbar.set_label('Surface Temperature (°C)', fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Overlay the CROCO boundary on the zoomed plot
ax_zoom.plot(lon_west, lat_west, 'b-', linewidth=2,
             label='CROCO Boundary', zorder=13)
ax_zoom.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=13)
ax_zoom.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=13)

# Add coastlines, land mask, and borders using add_feature
ax_zoom.add_feature(cfeature.LAND, zorder=10)
ax_zoom.coastlines(resolution='10m', color='black', linewidth=1, zorder=11)
ax_zoom.add_feature(cfeature.BORDERS, edgecolor='black',
                    linewidth=0.5, zorder=12)

# Add gridlines and labels
gl_zoom = ax_zoom.gridlines(crs=ccrs.PlateCarree(
), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl_zoom.top_labels = False
gl_zoom.right_labels = False

# Zoom in on the CROCO region
ax_zoom.set_extent([lon_rho.min(), lon_rho.max(),
                   lat_rho.min(), lat_rho.max()], crs=ccrs.PlateCarree())

# Set title, labels, and aspect ratio
ax_zoom.set_title(
    'GLORYS Surface Temperature - Zoomed to CROCO Boundary \n 01 January 2009', fontsize=16)
ax_zoom.set_xlabel('Longitude')
ax_zoom.set_ylabel('Latitude')

# Adjust layout to ensure both plots are the same size
plt.tight_layout()

# Save the combined plot
savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/glorys_surface_temp_with_zoom_jan_2009.png'
plt.savefig(savepath, dpi=300, bbox_inches='tight')

plt.show()

# %%


# Load GLORYS data
glorys_fname = "/mnt/e/GLORYS_SA_EEZ_ZOOM/GLORYS_Y2009M1.nc"
GLORYS_ds = xr.open_dataset(glorys_fname)

# Extract the surface temperature (first depth level, first time step)
surface_temp = GLORYS_ds['temp'].isel(depth=0, time=0)

# Extract the sea surface height (SSH) as a proxy for depth contours
ssh = GLORYS_ds['ssh'].isel(time=0)

# Function to get CROCO grid boundary and variables


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho


# Define CROCO file path and extract boundaries
croco_fname = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y2013M01.nc"
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho = get_croco_boundary(
    croco_fname)

# In situ data file paths and filenames
filenames = [
    "Betty's Bay_DAFF.nc",
    "Bordjies Deep_DAFF.nc",
    "Bordjies_DAFF.nc",
    "Doringbaai_DAFF.nc",
    "Fish Hoek_SAWS.nc",
    "Gordons Bay_SAWS.nc",
    "Kalk Bay_SAWS.nc",
    "Kommetjie_SAWS.nc",
    "Lamberts Bay_SAWS.nc",
    "Muizenberg_SAWS.nc",
    "Oudekraal_DAFF.nc",
    "Saldanha Bay_SAWS.nc",
    "Sea Point_SAWS.nc",
    "St Helena Bay_SAWS.nc",
    "Yzerfontein_SAWS.nc"
]

ERA5_VAL_path = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/AJSMIT_UTR/'

# Load in situ data
all_lat_insitu = []
all_lon_insitu = []

for filename in filenames:
    ds_ERA5 = xr.open_dataset(ERA5_VAL_path + "Validated_" + filename)

    lat_insitu = ds_ERA5.latitude.values
    lon_insitu = ds_ERA5.longitude.values

    all_lat_insitu.extend(lat_insitu)
    all_lon_insitu.extend(lon_insitu)

# Plot the surface temperature (Full GLORYS domain)
fig, axs = plt.subplots(1, 2, figsize=(18, 8), subplot_kw={
                        'projection': ccrs.PlateCarree()})

# Full domain plot
ax = axs[0]
contour_full = surface_temp.plot(
    ax=ax,
    cmap='turbo',
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8},
    vmin=0,
    vmax=35,
    levels=np.linspace(0, 40, 21)  # 20 intervals between 0 and 40 (21 levels)
)

cbar = contour_full.colorbar
cbar.set_label('Surface Temperature (°C)', fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Overlay the CROCO boundary on the full domain plot
ax.plot(lon_west, lat_west, 'b-', linewidth=2,
        label='CROCO Boundary', zorder=13)
ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=13)
ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=13)

# Plot in situ data points
# ax.plot(all_lon_insitu, all_lat_insitu, 'or', label='In situ', zorder=14, markeredgecolor='black', markeredgewidth=1, markersize=11)

# Add coastlines, land mask, and borders using add_feature
ax.add_feature(cfeature.LAND, zorder=10)
ax.coastlines(resolution='10m', color='black', linewidth=1, zorder=11)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5, zorder=12)

# Add gridlines and labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Set title, labels, and aspect ratio
ax.set_title(
    'GLORYS Surface Temperature - Full Domain \n 01 January 2009', fontsize=16)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Zoomed plot (within CROCO boundary)
ax_zoom = axs[1]

# Use correct coordinate names ('lonT' and 'latT')
croco_mask = (surface_temp['lonT'] >= lon_rho.min()) & (surface_temp['lonT'] <= lon_rho.max()) & \
             (surface_temp['latT'] >= lat_rho.min()) & (
                 surface_temp['latT'] <= lat_rho.max())

# Restrict the data to the CROCO boundary region
surface_temp_zoom = surface_temp.where(croco_mask, drop=True)

contour_zoom = surface_temp_zoom.plot(
    ax=ax_zoom,
    cmap='turbo',  # coolwarm
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8},
    vmin=0,
    vmax=35,
    levels=np.linspace(8, 25, 21)  # 20 intervals between 0 and 40 (21 levels)
)

cbar = contour_zoom.colorbar
cbar.set_label('Surface Temperature (°C)', fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Overlay the CROCO boundary on the zoomed plot
ax_zoom.plot(lon_west, lat_west, 'b-', linewidth=2,
             label='CROCO Boundary', zorder=13)
ax_zoom.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=13)
ax_zoom.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=13)

# Plot in situ data points on the zoomed plot
ax_zoom.plot(all_lon_insitu, all_lat_insitu, 'or', label='In situ',
             zorder=14, markeredgecolor='black', markeredgewidth=1, markersize=11)

# Add depth contours at 50, 100, and 500 meters
depth_contours = []
cs = ax_zoom.contour(GLORYS_ds['lonT'], GLORYS_ds['latT'], ssh,
                     levels=depth_contours, colors='k', linestyles='--', zorder=15)
ax_zoom.clabel(cs, inline=True, fontsize=10, fmt='%d m')

# Add coastlines, land mask, and borders using add_feature
ax_zoom.add_feature(cfeature.LAND, zorder=10)
ax_zoom.coastlines(resolution='10m', color='black', linewidth=1, zorder=11)
ax_zoom.add_feature(cfeature.BORDERS, edgecolor='black',
                    linewidth=0.5, zorder=12)

# Add gridlines and labels
gl_zoom = ax_zoom.gridlines(crs=ccrs.PlateCarree(
), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl_zoom.top_labels = False
gl_zoom.right_labels = False

# Zoom in on the CROCO region
ax_zoom.set_extent([lon_rho.min(), lon_rho.max(),
                   lat_rho.min(), lat_rho.max()], crs=ccrs.PlateCarree())

# Add custom legend entry
handles, labels = plt.gca().get_legend_handles_labels()
ax_zoom.legend(handles, labels, loc='upper right',
               bbox_to_anchor=(1, 1), fontsize=9)

# Set title, labels, and aspect ratio
ax_zoom.set_title(
    'GLORYS Surface Temperature - Zoomed to CROCO Boundary \n 01 January 2009', fontsize=16)
ax_zoom.set_xlabel('Longitude')
ax_zoom.set_ylabel('Latitude')

# Adjust layout to ensure both plots are the same size
plt.tight_layout()

# Save the combined plot
savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/glorys_surface_temp_with_zoom_jan_2009.png'
plt.savefig(savepath, dpi=300, bbox_inches='tight')

plt.show()

# %% Same as above but included all stations and datasets.


# Load GLORYS data
glorys_fname = "/mnt/e/GLORYS_SA_EEZ_ZOOM/GLORYS_Y2009M1.nc"
GLORYS_ds = xr.open_dataset(glorys_fname)

# Extract the surface temperature (first depth level, first time step)
surface_temp = GLORYS_ds['temp'].isel(depth=0, time=0)

# Extract the sea surface height (SSH) as a proxy for depth contours
ssh = GLORYS_ds['ssh'].isel(time=0)

# Function to get CROCO grid boundary and variables


def get_croco_boundary(fname):
    '''
    Return lon, lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds_ERA5:
        lon_rho = ds_ERA5.lon_rho.values
        lat_rho = ds_ERA5.lat_rho.values

    lon_west = lon_rho[:, 0]
    lat_west = lat_rho[:, 0]

    lon_south = lon_rho[-1, 1:-1]
    lat_south = lat_rho[-1, 1:-1]

    lon_north = lon_rho[0, 1:-1]
    lat_north = lat_rho[0, 1:-1]

    return lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho


# Define CROCO file path and extract boundaries
croco_fname = "/mnt/d/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_WASA3/output/croco_avg_Y2013M07.nc"
lon_west, lat_west, lon_south, lat_south, lon_north, lat_north, lon_rho, lat_rho = get_croco_boundary(
    croco_fname)

# Paths for each group
paths = {
    "AJSMIT_UTR": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/AJSMIT_UTR/',
    "GPITCHER_SHB": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/GPITCHER_SHB/',
    "ATAP": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/ATAP/',
    "Wirewalker_Mooring": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/Wirewalker/',
    "ADCP_Mooring": '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/adcp_mooring/'
}

savepath = '/home/nkululeko/somisana-croco/configs/swcape_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/'

# Initialize lists to store latitude and longitude values
all_lat_insitu = []
all_lon_insitu = []

# Separate lists for each group
gpitcher_lat_insitu = []
gpitcher_lon_insitu = []
atap_lat_insitu = []
atap_lon_insitu = []
wirewalker_lat_insitu = []
wirewalker_lon_insitu = []
adcp_lat_insitu = []
adcp_lon_insitu = []

# Single loop to process all files and accumulate latitude and longitude data
files_and_paths = {
    "AJSMIT_UTR": [
        "Betty's Bay_DAFF.nc",
        "Bordjies Deep_DAFF.nc",
        "Bordjies_DAFF.nc",
        "Doringbaai_DAFF.nc",
        "Fish Hoek_SAWS.nc",
        "Gordons Bay_SAWS.nc",
        "Kalk Bay_SAWS.nc",
        "Kommetjie_SAWS.nc",
        "Lamberts Bay_SAWS.nc",
        "Muizenberg_SAWS.nc",
        "Oudekraal_DAFF.nc",
        "Saldanha Bay_SAWS.nc",
        "Sea Point_SAWS.nc",
        "St Helena Bay_SAWS.nc",
        "Yzerfontein_SAWS.nc"
    ],
    "GPITCHER_SHB": [
        "WQM_20m_TS.nc",
        "WQM_70m_TS.nc"
    ],
    "ATAP": [
        "FalseBay_FB001.nc",
        "CapePoint_CP001.nc"
    ],
    "Wirewalker_Mooring": [
        "wirewalker_mooring_1.nc"
    ],
    "ADCP_Mooring": [
        "adcp_mooring_1.nc"
    ]
}

for group, filenames in files_and_paths.items():
    for filename in filenames:
        ds_ERA5 = xr.open_dataset(paths[group] + "Validated_" + filename)

        lat_insitu = ds_ERA5.latitude.values
        lon_insitu = ds_ERA5.longitude.values

        if group == "GPITCHER_SHB":
            # Store GPITCHER_SHB data separately
            gpitcher_lat_insitu.extend(lat_insitu)
            gpitcher_lon_insitu.extend(lon_insitu)
        elif group == "ATAP":
            # Store ATAP data separately
            atap_lat_insitu.extend(lat_insitu)
            atap_lon_insitu.extend(lon_insitu)
        elif group == "Wirewalker_Mooring":
            # Store Wirewalker_Mooring data separately
            wirewalker_lat_insitu.extend(lat_insitu)
            wirewalker_lon_insitu.extend(lon_insitu)
        elif group == "ADCP_Mooring":
            # Store ADCP_Mooring data separately
            adcp_lat_insitu.extend(lat_insitu)
            adcp_lon_insitu.extend(lon_insitu)
        else:
            # Store other data in the main lists
            all_lat_insitu.extend(lat_insitu)
            all_lon_insitu.extend(lon_insitu)

# Plot the surface temperature (Full GLORYS domain)
fig, axs = plt.subplots(1, 2, figsize=(18, 8), subplot_kw={
                        'projection': ccrs.PlateCarree()})

# Full domain plot
ax = axs[0]
contour_full = surface_temp.plot(
    ax=ax,
    cmap='turbo',
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8},
    vmin=0,
    vmax=35,
    levels=np.linspace(0, 40, 21)  # 20 intervals between 0 and 40 (21 levels)
)

cbar = contour_full.colorbar
cbar.set_label('Surface Temperature (°C)', fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Overlay the CROCO boundary on the full domain plot
ax.plot(lon_west, lat_west, 'b-', linewidth=2,
        label='CROCO Boundary', zorder=13)
ax.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=13)
ax.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=13)

# Plot in situ data points
ax.plot(all_lon_insitu, all_lat_insitu, 'or', label='In situ',
        zorder=14, markeredgecolor='black', markeredgewidth=1, markersize=11)

# Add coastlines, land mask, and borders using add_feature
ax.add_feature(cfeature.LAND, zorder=10)
ax.coastlines(resolution='10m', color='black', linewidth=1, zorder=11)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5, zorder=12)

# Add gridlines and labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Set title, labels, and aspect ratio
ax.set_title(
    'GLORYS Surface Temperature - Full Domain \n 01 July 2009', fontsize=16)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Zoomed plot (within CROCO boundary)
ax_zoom = axs[1]

# Use correct coordinate names ('lonT' and 'latT')
croco_mask = (surface_temp['lonT'] >= lon_rho.min()) & (surface_temp['lonT'] <= lon_rho.max()) & \
             (surface_temp['latT'] >= lat_rho.min()) & (
                 surface_temp['latT'] <= lat_rho.max())

# Restrict the data to the CROCO boundary region
surface_temp_zoom = surface_temp.where(croco_mask, drop=True)

contour_zoom = surface_temp_zoom.plot(
    ax=ax_zoom,
    cmap='turbo',
    cbar_kwargs={'label': 'Surface Temperature (°C)', 'shrink': 0.8},
    vmin=0,
    vmax=35,
    levels=np.linspace(8, 32, 21)  # 20 intervals between 0 and 40 (21 levels)
)

cbar = contour_zoom.colorbar
cbar.set_label('Surface Temperature (°C)', fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Overlay the CROCO boundary on the zoomed plot
ax_zoom.plot(lon_west, lat_west, 'b-', linewidth=2,
             label='CROCO Boundary', zorder=13)
ax_zoom.plot(lon_south, lat_south, 'b-', linewidth=2, zorder=13)
ax_zoom.plot(lon_north[::-1], lat_north[::-1], 'b-', linewidth=2, zorder=13)

# Plot in situ data points on the zoomed plot
ax_zoom.plot(all_lon_insitu, all_lat_insitu, 'or', label='In situ',
             zorder=14, markeredgecolor='black', markeredgewidth=1, markersize=11)

# Plot GPITCHER_SHB stations (orange stars)
ax_zoom.plot(gpitcher_lon_insitu, gpitcher_lat_insitu, '*', color='orange',
             label='GPITCHER_SHB Stations', zorder=18, markeredgecolor='black', markeredgewidth=1, markersize=20)

# Plot ATAP stations (yellow squares)
ax_zoom.plot(atap_lon_insitu, atap_lat_insitu, 's', color='yellow', label='ATAP Stations',
             zorder=16, markeredgecolor='black', markeredgewidth=1, markersize=9)

# Plot Wirewalker_Mooring and ADCP_Mooring stations (white triangles)
ax_zoom.plot(wirewalker_lon_insitu, wirewalker_lat_insitu, '^', color='white',
             label='Wirewalker Mooring', zorder=17, markeredgecolor='black', markeredgewidth=1, markersize=11)
ax_zoom.plot(adcp_lon_insitu, adcp_lat_insitu, '^', color='white', label='ADCP Mooring',
             zorder=18, markeredgecolor='black', markeredgewidth=1, markersize=11)

# Add depth contours at 50, 100, and 500 meters
depth_contours = [50, 100, 500]
cs = ax_zoom.contour(GLORYS_ds['lonT'], GLORYS_ds['latT'], ssh,
                     levels=depth_contours, colors='k', linestyles='--', zorder=19)
ax_zoom.clabel(cs, inline=True, fontsize=10, fmt='%d m')

# Add coastlines, land mask, and borders using add_feature
ax_zoom.add_feature(cfeature.LAND, zorder=10)
ax_zoom.coastlines(resolution='10m', color='black', linewidth=1, zorder=11)
ax_zoom.add_feature(cfeature.BORDERS, edgecolor='black',
                    linewidth=0.5, zorder=12)

# Add gridlines and labels
gl_zoom = ax_zoom.gridlines(crs=ccrs.PlateCarree(
), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl_zoom.top_labels = False
gl_zoom.right_labels = False

# Zoom in on the CROCO region
ax_zoom.set_extent([lon_rho.min(), lon_rho.max(),
                   lat_rho.min(), lat_rho.max()], crs=ccrs.PlateCarree())

# Add custom legend entry
handles, labels = plt.gca().get_legend_handles_labels()
ax_zoom.legend(handles, labels, loc='upper right',
               bbox_to_anchor=(1, 1), fontsize=9)

# Set title, labels, and aspect ratio
ax_zoom.set_title(
    'GLORYS Surface Temperature within CROCO Boundary \n 01 July 2009', fontsize=16)
ax_zoom.set_xlabel('Longitude')
ax_zoom.set_ylabel('Latitude')

# Adjust layout to ensure both plots are the same size
plt.tight_layout()

# Save the combined plot
plt.savefig(savepath + 'glorys_surface_temp_with_zoom_jan_2009.png',
            dpi=300, bbox_inches='tight')

plt.show()
