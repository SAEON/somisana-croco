import numpy as np
import xarray as xr
import crocotools_py.postprocess as post
# from glob import glob
from datetime import datetime
# import netCDF4 as nc
# import pandas as pd
# import os
# import numpy as np
# import cftime
from colorama import Fore, Style

def statistics(model_data,data_obs_model_timeaxis):
    """
            This Function calculates all model validation statistics
            
            Parameters:
            - model_data:                       Dataset of the CROCO model
            - data_obs_model_timeaxis:          Dataset of the insitu station already modigied to the 
                                                model time axis and depth levels are matching
            Returns:
            - insitu_correlation, insitu_rmse, insitu_mean_diff, insitu_min_value, insitu_max_value,
                    model_correlation, model_rmse, model_mean_diff, model_min_value, model_max_value, model_total_bias
        """
    # Convert xarray DataArrays to numpy arrays        
    model_data = model_data.values
    insitu_data = data_obs_model_timeaxis.values   
    
    print(f'model_data:{model_data}')
    print(f'insitu_data:{insitu_data}')
          
    # Find rows where all values are not NaN in insitu_data & Filter both model_data and insitu_data using the non-NaN indices
    non_nan_indices = ~np.isnan(insitu_data)
    print(f'non_nan_indices: {non_nan_indices}')
    model_data = model_data[non_nan_indices]
    insitu_data = insitu_data[non_nan_indices]
    
    # Calculate statistics for insitu_data
    insitu_correlation = np.corrcoef(insitu_data.flatten(), model_data.flatten())[0, 1]
    insitu_squared_diff = (insitu_data - model_data) ** 2
    insitu_rmse = np.sqrt(np.mean(insitu_squared_diff))
    insitu_mean_diff = np.mean(insitu_data)
    insitu_min_value = np.min(insitu_data) if len(insitu_data) > 0 else np.nan # This handles the cases where min is calculated from an empty dataset
    insitu_max_value = np.max(insitu_data) if len(insitu_data) > 0 else np.nan # This handles the cases where max is calculated from an empty dataset
    
    # Calculate statistics for model_data
    model_correlation = np.corrcoef(model_data.flatten(), insitu_data.flatten())[0, 1]
    model_squared_diff = (model_data - insitu_data) ** 2
    model_rmse = np.sqrt(np.mean(model_squared_diff))
    model_mean_diff = np.mean(model_data)
    model_min_value = np.min(model_data) if len(model_data) > 0 else np.nan # This handles the cases where min is calculated from an empty dataset
    model_max_value = np.max(model_data) if len(model_data) > 0 else np.nan # This handles the cases where max is calculated from an empty dataset
    model_total_bias = np.mean(np.abs(model_data - insitu_data))
    
    print(f'insitu_correlation:{insitu_correlation}')
    print(f'insitu_rmse,:{insitu_rmse}')
    print(f'insitu_mean_diff,:{insitu_mean_diff}') 
     # , insitu_min_value, insitu_max_value,
            # model_correlation, model_rmse, model_mean_diff, model_min_value, model_max_value, model_total_bias')
    
    return (insitu_correlation, insitu_rmse, insitu_mean_diff, insitu_min_value, insitu_max_value,
            model_correlation, model_rmse, model_mean_diff, model_min_value, model_max_value, model_total_bias)

def obs_2_model_timeaxis(ds_obs, ds_mod):
    """
    get observations (or other) time onto the model time axis
    One could of course use xarray's built-in interp function, which includes time averaging over predefined windows
    but there is always the risk of not getting the observation time axis identical to the model, even if you get the same number of time-steps out
    This function at least ensures the times line-up perfectly    
    
    inputs:
        ds_obs: dataset/dataarray of the observations
        ds_mod: dataset/dataarray of the model
    
    returns an xarray dataset/dataarray of the observations on the model time axis, including time averaging over each model time-step                                            
    """
    
    model_time = ds_mod['time']
    
    # Compute bin edges of the model times
    model_dt = model_time[1] - model_time[0] # we can assume the model timestep is constant 
    bin_edges = model_time + model_dt/2
    bin_edges = xr.concat(
        [bin_edges[0] - model_dt, bin_edges],
        dim="time"
    )
    
    # Use groupby_bins to average ds_obs values over each bin
    ds_obs_model_timeaxis = ds_obs.groupby_bins("time", bin_edges, labels=model_time.values).mean(skipna=True)
    
    ds_obs_model_timeaxis = ds_obs_model_timeaxis.rename({"time_bins": "time"})
        
    print('________')
    print(f'ds_obs_model_timeaxis:{ds_obs_model_timeaxis}')
    print('________')
    
    return ds_obs_model_timeaxis

def extract_lat_lon(ds):
    """
            This Function extracts longitudes and latitudes from a dataset. The reason why it is important is 
            because you never know what numecluture structure was used for writing the grid coordinates
            it could be full names or lat and lon, or first letter could be capitalized if not all letters.

            Parameters:
            - dataset:          Dataset usually the insitu ds as this is where for validation we extract
                                the location parameters of the validation station.

            Returns:
            - long_obs, lat_obs
        """
    try:
        long_obs = ds.lon.values
        lat_obs = ds.lat.values

        if long_obs.size > 0 and (long_obs == [0]).any():
            raise ValueError(
                "Undesired output detected in 'longitude': array([0])")
        if lat_obs.size > 0 and (lat_obs == [0]).any():
            raise ValueError(
                "Undesired output detected in 'latitude': array([0])")

    except ValueError as ve:
        print(f"Caught an exception: {ve}")
        # Handle the exception, provide default values
        try:
            long_obs = ds.lon.values
            lat_obs = ds.lat.values
        except ValueError as ve:
            print(f"Caught an exception: {ve}")
            # Handle the exception, provide default values
            try:
                long_obs = ds.Lon.values
                lat_obs = ds.Lat.values
            except ValueError as ve:
                print(f"Caught an exception: {ve}")
                # Handle the exception, provide default values
                try:
                    long_obs = ds.Longitude.values
                    lat_obs = ds.Latitude.values
                except ValueError as ve:
                    print(f"Caught an exception: {ve}")
                    # Handle the exception, provide default values
                    long_obs = ds.lon.item()
                    lat_obs = ds.lat.item()
    else:
        print("Longitude and latitude successfully retrieved.")
    return long_obs, lat_obs

def get_model_obs_ts(fname, fname_obs, output_path, var, depth=-1, i_shifted=0, j_shifted=0, ref_date=None, lon_extract=None):
    """
           This Function is for retrieving all datasets from the other functions and package them into a netCDF 
           file on the same time axis: 

            Parameters:
            - fname             :filename of the model
            - fname_obs         :filename of observations
            - output_path       :filename of fname_out
            - var               :variable input by the user. should be the same in both the modeland obs netCDFs.
            - lat               :lat read from the insitu input
            - lon               :lon read from the insitu input
            - ref_date          :static user input
            - depth             :user input based on insitu sensor depth
            - i_shifted         :optional user inputs if the z level of the model does not match the input lon depth
            - j_shifted         :optional user inputs if the z level of the model does not match the input lat depth
            - time_lims         :limits are computed based on the length of the insitu data that matches model span
    """
    print(f'fname: {fname}')
    print(f'fname_obs: {fname_obs}')
    print(f'output_path: {output_path}')
    print(f'var: {var}')
    print(f'depth: {depth}')
    print(f'i_shifted: {i_shifted}')
    print(f'ref_date: {ref_date}')
    print(f'lon_extract: {lon_extract}')
    
    # Load the NetCDF dataset
    ds_obs = xr.open_dataset(fname_obs).sel(depth=abs(depth))
    
    # Some observations have a different time format, that makes tem compatible
    ds_obs['time'] = ds_obs['time'].astype('datetime64[ns]') 

    # Define the time limits
    # time_lims = [datetime(2011, 3, 1), datetime(2011, 3, 31)]
    # time_lims = [datetime(2013, 1, 1), datetime(2011, 12, 31)]

    # Slice the dataset based on the time limits
    # ds_obs = ds_obs.sel(time=slice(*time_lims))
    
        
    # The following if statement handles the u,v grid. The else part handles the rho grid points computation. 
    # if var=='u' or var=='v':
    if (var != 'temp') or (var != 'salt'):
        long_obs, lat_obs = extract_lat_lon(ds_obs)
        model_data = post.get_ts_uv(fname, long_obs, lat_obs, ref_date, 
                        i_shift=0, j_shift=0, 
                        level=depth,
                        )
        
        da_obs_u = ds_obs.u
        da_obs_v = ds_obs.v
        # data_obs_model_timeaxis = obs_2_model_timeaxis(da_obs, model_data) 
        
        insitu_u_orig = obs_2_model_timeaxis(da_obs_u, model_data)
        insitu_v_orig = obs_2_model_timeaxis(da_obs_v, model_data)

        # Convert xarray DataArrays to numpy arrays      
        model_u_orig = model_data.u.values
        model_v_orig = model_data.v.values
        
        # longitude ,latitude = extract_lat_lon(ds_obs)
        longitude = ds_obs.lon.values
        latitude = ds_obs.lat.values

        # depth_levels = insitu_u_orig.depth.values
        
        # Assuming the number of depth levels is `n_depth_levels`
        n_depth_levels = len([depth]) #len(depth_levels)
        
        # Initialize empty arrays for each statistical measure. This will append on the for loop runs.
        insitu_correlation_u = np.empty(n_depth_levels)
        insitu_rmse_u = np.empty(n_depth_levels)
        insitu_mean_diff_u = np.empty(n_depth_levels)
        insitu_min_value_u = np.empty(n_depth_levels)
        insitu_max_value_u = np.empty(n_depth_levels)
        
        model_correlation_u = np.empty(n_depth_levels)
        model_rmse_u = np.empty(n_depth_levels)
        model_mean_diff_u = np.empty(n_depth_levels)
        model_min_value_u = np.empty(n_depth_levels)
        model_max_value_u = np.empty(n_depth_levels)
        model_total_bias_u = np.empty(n_depth_levels)     
        
        insitu_correlation_v = np.empty(n_depth_levels)
        insitu_rmse_v = np.empty(n_depth_levels)
        insitu_mean_diff_v = np.empty(n_depth_levels)
        insitu_min_value_v = np.empty(n_depth_levels)
        insitu_max_value_v = np.empty(n_depth_levels)
        
        model_correlation_v = np.empty(n_depth_levels)
        model_rmse_v = np.empty(n_depth_levels)
        model_mean_diff_v = np.empty(n_depth_levels)
        model_min_value_v = np.empty(n_depth_levels)
        model_max_value_v = np.empty(n_depth_levels)
        model_total_bias_v = np.empty(n_depth_levels)  
        
        print(f'model_data: {model_data}')
        # Loop through depth levels
        print([depth])
        if len(depth) == 1:
            model_data_depth = model_data[var][:]
            obs_data_depth = insitu_u_orig[:]
                        
            (insitu_correlation_u, insitu_rmse_u, insitu_mean_diff_u,
             insitu_min_value_u, insitu_max_value_u,
             model_correlation_u, model_rmse_u, model_mean_diff_u,
             model_min_value_u, model_max_value_u, model_total_bias_u
             ) = statistics(model_data_depth, obs_data_depth)
        
        # Loop through depth levels
        # for idx_depth, depth_level in enumerate(depth_levels):
            model_data_depth = model_data[var][:]
            obs_data_depth = insitu_v_orig[:]
                        
            (insitu_correlation_v, insitu_rmse_v, insitu_mean_diff_v,
             insitu_min_value_v, insitu_max_value_v,
             model_correlation_v, model_rmse_v, model_mean_diff_v,
             model_min_value_v, model_max_value_v, model_total_bias_v
             ) = statistics(model_data_depth, obs_data_depth)
        
        else:
            for idx_depth, depth_level in enumerate([depth]):
                print('')
                print(f'model_data[var][:, idx_depth]: {model_data[var][:, idx_depth]}')
                print('')
                model_data_depth = model_data[var][:, idx_depth]
                obs_data_depth = insitu_u_orig[:, idx_depth]
                            
                (insitu_correlation_u[idx_depth], insitu_rmse_u[idx_depth], insitu_mean_diff_u[idx_depth],
                 insitu_min_value_u[idx_depth], insitu_max_value_u[idx_depth],
                 model_correlation_u[idx_depth], model_rmse_u[idx_depth], model_mean_diff_u[idx_depth],
                 model_min_value_u[idx_depth], model_max_value_u[idx_depth], model_total_bias_u[idx_depth]
                 ) = statistics(model_data_depth, obs_data_depth)
            
            # Loop through depth levels
            # for idx_depth, depth_level in enumerate(depth_levels):
                model_data_depth = model_data[var][:, idx_depth]
                obs_data_depth = insitu_v_orig[:, idx_depth]
                            
                (insitu_correlation_v[idx_depth], insitu_rmse_v[idx_depth], insitu_mean_diff_v[idx_depth],
                 insitu_min_value_v[idx_depth], insitu_max_value_v[idx_depth],
                 model_correlation_v[idx_depth], model_rmse_v[idx_depth], model_mean_diff_v[idx_depth],
                 model_min_value_v[idx_depth], model_max_value_v[idx_depth], model_total_bias_v[idx_depth]
                 ) = statistics(model_data_depth, obs_data_depth)

        # model_data = model_data[var].values
        insitu_data = insitu_u_orig.values
        # depth = insitu_data.shape #np.arange(insitu_data.shape)
        
        # Create a new time array based on the insitu data timeseries
        time = insitu_u_orig.time.values
        
        # eta_rho = insitu_u_orig.eta_rho.values
        # xi_rho  = insitu_u_orig.xi_rho.values
        # lon_rho = insitu_u_orig.lon_rho.values
        # lat_rho = insitu_u_orig.lat_rho.values
        
        eta_rho = model_data.eta_rho.values
        xi_rho  = model_data.xi_rho.values
        lon_rho = model_data.lon_rho.values
        lat_rho = model_data.lat_rho.values
        
        print(f'time: {time}')
        print(f'depth: {depth}')
        print(f'eta_rho: {eta_rho}')
        print(f'xi_rho: {xi_rho}')
        print(f'lon_rho: {lon_rho}')
        print(f'lat_rho: {lat_rho}')

        # # Create xarray DataArrays for the entire dataset  
        # model_da = xr.DataArray(model_data[var].values, 
        #                         dims=("time", "depth"), 
        #                         coords={"time": time, "depth": depth,"eta_rho":eta_rho,"xi_rho":xi_rho , "lon_rho":lon_rho,"lat_rho":lat_rho})
        
        
        # insitu_da = xr.DataArray(insitu_data, 
        #                          dims=("time", "depth"), 
        #                          coords={"time": time, "depth": depth,"eta_rho":eta_rho,"xi_rho":xi_rho , "lon_rho":lon_rho,"lat_rho":lat_rho})
        
        # # Create a Dataset with the DataArrays
        # ds = xr.Dataset({f"insitu_{var}": insitu_da, f"model_{var}": model_da})
        
        # ds = ds.assign_coords({"longitude":longitude})
        # ds = ds.assign_coords({"latitude":latitude})
        
        # #Global attributes        
        # ds.attrs["longitude"] = longitude
        # ds.attrs["latitude"] = latitude   
        # ds.attrs["h"] = model_data.h.values
        
    
        
        # depth = np.array([depth])  # or depth = [-30]

        # Create xarray DataArrays for u component
        # depth = insitu_u_orig.shape
        time = model_data.time.values
        # insitu_da_u = xr.DataArray(insitu_u_orig.values[:,np.newaxis].T, dims=("time", "depth"), coords={"time": time, "depth": depth})
        # model_da_u = xr.DataArray(model_u_orig[:,np.newaxis].T, dims=("time", "depth"), coords={"time": time, "depth": depth})

        # # insitu_da_u = xr.DataArray(insitu_u_orig, dims=("time", "depth"), coords={"time": time, "depth": depth})
        # # model_da_u = xr.DataArray(model_u_orig, dims=("time", "depth"), coords={"time": time, "depth": depth})

        
        # # Create xarray DataArrays for v component
        # insitu_da_v = xr.DataArray(insitu_v_orig.values[:,np.newaxis].T, dims=("time", "depth"), coords={"time": time, "depth": depth})
        # model_da_v = xr.DataArray(model_v_orig[:,np.newaxis].T, dims=("time", "depth"), coords={"time": time, "depth": depth})
        
        # # Create a Dataset with the DataArrays for u and v components
        # ds = xr.Dataset({"insitu_u": insitu_da_u, "insitu_v": insitu_da_v, "model_u": model_da_u, "model_v": model_da_v})
        
        
        if np.size(depth)==1:    
            insitu_da_u = xr.DataArray(insitu_u_orig.values[:, np.newaxis], dims=("time", "depth"), coords={"time": time, "depth": abs(depth)})
            model_da_u  = xr.DataArray(model_u_orig[:, np.newaxis], dims=("time", "depth"), coords={"time": time, "depth": abs(depth)})
            
            insitu_da_v = xr.DataArray(insitu_v_orig.values[:, np.newaxis], dims=("time", "depth"), coords={"time": time, "depth": abs(depth)})
            model_da_v  = xr.DataArray(model_v_orig[:, np.newaxis], dims=("time", "depth"), coords={"time": time, "depth": abs(depth)})
            
            ds = xr.Dataset({
                "insitu_u": insitu_da_u,
                "insitu_v": insitu_da_v,
                "model_u": model_da_u,
                "model_v": model_da_v
            })
        elif np.size(depth)>1:
            insitu_da_u = xr.DataArray(insitu_u_orig.values, dims=("time", "depth"), coords={"time": time, "depth": abs(depth)})
            model_da_u  = xr.DataArray(model_u_orig, dims=("time", "depth"), coords={"time": time, "depth": abs(depth)})
            
            insitu_da_v = xr.DataArray(insitu_v_orig.values, dims=("time", "depth"), coords={"time": time, "depth": abs(depth)})
            model_da_v  = xr.DataArray(model_v_orig, dims=("time", "depth"), coords={"time": time, "depth": abs(depth)})
            
            ds = xr.Dataset({
                "insitu_u": insitu_da_u,
                "insitu_v": insitu_da_v,
                "model_u": model_da_u,
                "model_v": model_da_v
            })        

        print("ds",ds)
        # Calculate seasonal means for insitu_data
        insitu_seasonal_means = {}
        for season, months in {'JFM': [1, 2, 3], 'AMJ': [4, 5, 6], 'JAS': [7, 8, 9], 'OND': [10, 11, 12]}.items():
            seasonal_data_u = insitu_da_u.sel(time=insitu_da_u.time.dt.month.isin(months))
            seasonal_data_v = insitu_da_v.sel(time=insitu_da_v.time.dt.month.isin(months))
            seasonal_data_magnitude = np.sqrt(seasonal_data_u**2 + seasonal_data_v**2)
            insitu_seasonal_means[season] = seasonal_data_magnitude.mean(dim=["time", "depth"]).item()
        
        # Calculate seasonal means for model_data
        model_seasonal_means = {}
        for season, months in {'JFM': [1, 2, 3], 'AMJ': [4, 5, 6], 'JAS': [7, 8, 9], 'OND': [10, 11, 12]}.items():
            seasonal_data_u = model_da_u.sel(time=model_da_u.time.dt.month.isin(months))
            seasonal_data_v = model_da_v.sel(time=model_da_v.time.dt.month.isin(months))
            seasonal_data_magnitude = np.sqrt(seasonal_data_u**2 + seasonal_data_v**2)
            model_seasonal_means[season] = seasonal_data_magnitude.mean(dim=["time", "depth"]).item()
        
        # Add attributes for seasonal means to the Dataset
        ds.attrs.update({
            f"insitu_seasonal_mean_{season}": mean_value for season, mean_value in insitu_seasonal_means.items()
        })
        ds.attrs.update({
            f"model_seasonal_mean_{season}": mean_value for season, mean_value in model_seasonal_means.items()
        })
        
        #Global attributes        
        ds.attrs["longitude"] = longitude
        ds.attrs["latitude"] = latitude   
        ds.attrs["h"] = model_data.h.values
        
        # Add statistics as variables to the Dataset
        ds["insitu_correlation_u"] = insitu_correlation_u
        ds["insitu_rmse_u"] = insitu_rmse_u
        ds["insitu_mean_difference_u"] = insitu_mean_diff_u
        ds["insitu_min_value_u"] = insitu_min_value_u
        ds["insitu_max_value_u"] = insitu_max_value_u
        
        ds["model_correlation_u"] = model_correlation_u
        ds["model_rmse_u"] = model_rmse_u
        ds["model_mean_difference_u"] = model_mean_diff_u
        ds["model_min_value_u"] = model_min_value_u
        ds["model_max_value_u"] = model_max_value_u
        ds["model_total_bias_u"] = model_total_bias_u
        
        # Add statistics as variables to the Dataset
        ds["insitu_correlation_v"] = insitu_correlation_v
        ds["insitu_rmse_v"] = insitu_rmse_v
        ds["insitu_mean_difference_v"] = insitu_mean_diff_v
        ds["insitu_min_value_v"] = insitu_min_value_v
        ds["insitu_max_value_v"] = insitu_max_value_v
        
        ds["model_correlation_v"] = model_correlation_v
        ds["model_rmse_v"] = model_rmse_v
        ds["model_mean_difference_v"] = model_mean_diff_v
        ds["model_min_value_v"] = model_min_value_v
        ds["model_max_value_v"] = model_max_value_v
        ds["model_total_bias_V"] = model_total_bias_v
        
        print("ds",ds)
        
        # Save the Dataset to a NetCDF file
        ds.to_netcdf(output_path)
        
        print("______________________________________________________________")
        print("______________________________________________________________")            
        print("Your output file will look like this",ds)
        print("Longitude",longitude)
        print("Latitude",latitude)

        print("______________________________________________________________")
        print("______________________________________________________________")
        print("NetCDF file saved successfully: This is a u,v grid, full water column dataset i.e. ADCP.")
        print("Output is at: ", Style.BRIGHT + Fore.GREEN + output_path + Style.RESET_ALL)
        
        
        
        print("______________________________________________________________")
        print("______________________________________________________________")            
        print("Your output file will look like this",ds.attrs)    
        print("______________________________________________________________")
        print("______________________________________________________________")
        print("NetCDF file saved successfully: This is a u,v grid, full water column dataset i.e. ADCP.")
        print("Output is at: ", Style.BRIGHT + Fore.GREEN + output_path + Style.RESET_ALL)
    
    #Not a u-v dataset    
    else:
        if "depth" in ds_obs.dims and len(ds_obs.coords["depth"]) > 1:
            long_obs = ds_obs.lon.values
            lat_obs = ds_obs.lat.values
        else:
            long_obs,lat_obs = extract_lat_lon(ds_obs)
            
        model_data = post.get_ts(fname, var, long_obs, lat_obs, ref_date,
                             i_shift=0, j_shift=0,
                             time_lims=slice(None),
                             depths=depth,  # slice(None),
                             default_to_bottom=True
                             )

        da_obs = ds_obs[var][:]            
        data_obs_model_timeaxis  = obs_2_model_timeaxis(da_obs, model_data)    
        
        #
        if "depth" in ds_obs.dims and len(ds_obs.coords["depth"]) > 1:
            longitude = ds_obs.lon.values
            latitude = ds_obs.lat.values
        else:
            longitude,latitude =extract_lat_lon(ds_obs)
        
        #if the only dimension is time, then this is a single depth dataset like ATAP
        if data_obs_model_timeaxis.dims == ("time",):            
            model_data = model_data[var][:]
            
            # Calling the statistics function
            (insitu_correlation, insitu_rmse, insitu_mean_diff,
             insitu_min_value, insitu_max_value,
             model_correlation, model_rmse, model_mean_diff,
             model_min_value, model_max_value, model_total_bias
            ) = statistics(model_data, data_obs_model_timeaxis)
            
            model_data = model_data.values
            insitu_data = data_obs_model_timeaxis.values            
            time = data_obs_model_timeaxis.time.values
            insitu_da = xr.DataArray(insitu_data, dims=("time"), coords={"time": time})
            model_da = xr.DataArray(model_data, dims=("time"), coords={"time": time})
            
            # Create a Dataset with the DataArrays
            ds = xr.Dataset({f"insitu_data_{var}": insitu_da, f"model_data_{var}": model_da})
            
            #Global attributes        
            ds.attrs["longitude"] = longitude
            ds.attrs["latitude"] = latitude 
            
            ds = ds.assign_coords({"longitude":longitude})
            ds = ds.assign_coords({"latitude":latitude})
            ds = ds.assign_coords({"depth":depth})
            # ds = ds.assign_coords({"h":h.values})
            
            # Add attributes for statistics to the Dataset
            ds.attrs["insitu_correlation_coefficient"] = insitu_correlation
            ds.attrs["insitu_rmse"] = insitu_rmse
            ds.attrs["insitu_mean_difference"] = insitu_mean_diff
            ds.attrs["insitu_min_value"] = insitu_min_value
            ds.attrs["insitu_max_value"] = insitu_max_value
            
            ds.attrs["model_correlation_coefficient"] = model_correlation
            ds.attrs["model_rmse"] = model_rmse
            ds.attrs["model_mean_difference"] = model_mean_diff
            ds.attrs["model_min_value"] = model_min_value
            ds.attrs["model_max_value"] = model_max_value
            ds.attrs["model_total_bias"] = model_total_bias
            ds.attrs["h"] = ds.h.values
            
            # Calculate seasonal means for insitu_data
            insitu_seasonal_means = {}
            for season, months in {'JFM': [1, 2, 3], 'AMJ': [4, 5, 6], 'JAS': [7, 8, 9], 'OND': [10, 11, 12]}.items():
                seasonal_data = insitu_da.sel(time=insitu_da.time.dt.month.isin(months))
                insitu_seasonal_means[season] = seasonal_data.mean(dim=["time"]).item()
            
            # Calculate seasonal means for model_data
            model_seasonal_means = {}
            for season, months in {'JFM': [1, 2, 3], 'AMJ': [4, 5, 6], 'JAS': [7, 8, 9], 'OND': [10, 11, 12]}.items():
                seasonal_data = model_da.sel(time=model_da.time.dt.month.isin(months))
                model_seasonal_means[season] = seasonal_data.mean(dim=["time"]).item()
            
            # Add attributes for seasonal means to the Dataset
            ds.attrs.update({
                f"insitu_seasonal_mean_{season}": mean_value for season, mean_value in insitu_seasonal_means.items()
            })
            ds.attrs.update({
                f"model_seasonal_mean_{season}": mean_value for season, mean_value in model_seasonal_means.items()
            })
            
            # Save the Dataset to a NetCDF file
            ds.to_netcdf(output_path)
            
            print("______________________________________________________________")
            print("______________________________________________________________")            
            print("Your output file will look like this",ds)        
            print("______________________________________________________________")
            print("______________________________________________________________")
            print("NetCDF file saved successfully: This is a rho grid, single water depth dataset i.e. ATAP.")
            print("Output is at: ", Style.BRIGHT + Fore.GREEN + output_path + Style.RESET_ALL)
        
        #else if the data has time and depth dimentions, it is a water column data like wirewalker     
        else:

            depth_levels = data_obs_model_timeaxis.depth.values
            
            # Assuming the number of depth levels is `n_depth_levels`
            n_depth_levels = len(depth_levels)
            
            # Initialize empty arrays for each statistical measure. This will append on the for loop runs.
            insitu_correlation = np.empty(n_depth_levels)
            insitu_rmse = np.empty(n_depth_levels)
            insitu_mean_diff = np.empty(n_depth_levels)
            insitu_min_value = np.empty(n_depth_levels)
            insitu_max_value = np.empty(n_depth_levels)
            
            model_correlation = np.empty(n_depth_levels)
            model_rmse = np.empty(n_depth_levels)
            model_mean_diff = np.empty(n_depth_levels)
            model_min_value = np.empty(n_depth_levels)
            model_max_value = np.empty(n_depth_levels)
            model_total_bias = np.empty(n_depth_levels)          
            
            # Loop through depth levels
            for idx_depth, depth_level in enumerate(depth_levels):
                model_data_depth = model_data[var][:, idx_depth]
                obs_data_depth = data_obs_model_timeaxis[:, idx_depth]
                            
                (insitu_correlation[idx_depth], insitu_rmse[idx_depth], insitu_mean_diff[idx_depth],
                 insitu_min_value[idx_depth], insitu_max_value[idx_depth],
                 model_correlation[idx_depth], model_rmse[idx_depth], model_mean_diff[idx_depth],
                 model_min_value[idx_depth], model_max_value[idx_depth], model_total_bias[idx_depth]
                 ) = statistics(model_data_depth, obs_data_depth)

            # model_data = model_data[var].values
            insitu_data = data_obs_model_timeaxis.values
            depth = np.arange(insitu_data.shape[1])
            
            # Create a new time array based on the insitu data timeseries
            time = data_obs_model_timeaxis.time.values
            
            eta_rho = data_obs_model_timeaxis.eta_rho.values
            xi_rho  = data_obs_model_timeaxis.xi_rho.values
            lon_rho = data_obs_model_timeaxis.lon_rho.values
            lat_rho = data_obs_model_timeaxis.lat_rho.values

            # Create xarray DataArrays for the entire dataset  
            model_da = xr.DataArray(model_data[var].values, dims=("time", "depth"), coords={"time": time, "depth": depth,"eta_rho":eta_rho,"xi_rho":xi_rho , "lon_rho":lon_rho,"lat_rho":lat_rho})
            insitu_da = xr.DataArray(insitu_data, dims=("time", "depth"), coords={"time": time, "depth": depth,"eta_rho":eta_rho,"xi_rho":xi_rho , "lon_rho":lon_rho,"lat_rho":lat_rho})
            
            # Create a Dataset with the DataArrays
            ds = xr.Dataset({f"insitu_{var}": insitu_da, f"model_{var}": model_da})
            
            ds = ds.assign_coords({"longitude":longitude})
            ds = ds.assign_coords({"latitude":latitude})
            
            #Global attributes        
            ds.attrs["longitude"] = longitude
            ds.attrs["latitude"] = latitude   
            ds.attrs["h"] = model_data.h.values
            
            # Add statistics as variables to the Dataset
            ds["insitu_correlation"] = insitu_correlation
            ds["insitu_rmse"] = insitu_rmse
            ds["insitu_mean_difference"] = insitu_mean_diff
            ds["insitu_min_value"] = insitu_min_value
            ds["insitu_max_value"] = insitu_max_value
            
            ds["model_correlation"] = model_correlation
            ds["model_rmse"] = model_rmse
            ds["model_mean_difference"] = model_mean_diff
            ds["model_min_value"] = model_min_value
            ds["model_max_value"] = model_max_value
            ds["model_total_bias"] = model_total_bias
            
            # Calculate seasonal means for insitu_data & model_data
            insitu_seasonal_means = {}
            for season, months in {'JFM': [1, 2, 3], 'AMJ': [4, 5, 6], 'JAS': [7, 8, 9], 'OND': [10, 11, 12]}.items():
                seasonal_data = insitu_da.sel(time=insitu_da.time.dt.month.isin(months))
                insitu_seasonal_means[season] = seasonal_data.mean(dim=["time"]).values

            model_seasonal_means = {}
            for season, months in {'JFM': [1, 2, 3], 'AMJ': [4, 5, 6], 'JAS': [7, 8, 9], 'OND': [10, 11, 12]}.items():
                seasonal_data = model_da.sel(time=model_da.time.dt.month.isin(months))
                model_seasonal_means[season] = seasonal_data.mean(dim=["time"]).values
            
            # Add attributes for seasonal means to the Dataset
            ds.attrs.update({
                f"insitu_seasonal_mean_{season}": mean_value for season, mean_value in insitu_seasonal_means.items()
            })
            ds.attrs.update({
                f"model_seasonal_mean_{season}": mean_value for season, mean_value in model_seasonal_means.items()
            })
            
            # Save the Dataset to a NetCDF file
            ds.to_netcdf(output_path)
            
            print("______________________________________________________________")
            print("______________________________________________________________")            
            print("Your output file will look like this",ds)    
            print("______________________________________________________________")
            print("______________________________________________________________")
            print("NetCDF file saved successfully. i.e Wirewalker")
            print("Output is at: ", Style.BRIGHT + Fore.GREEN + output_path + Style.RESET_ALL)
    
    print("Good Job!!!.")
    
