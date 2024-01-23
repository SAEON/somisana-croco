import numpy as np
import xarray as xr 
from glob import glob
from datetime import datetime
import postprocess as post
import netCDF4 as nc
import pandas as pd
import os

def get_ts_obs(fname_obs, var):
    """
           The next Function is for loading IN SITU time series:
               1. obs       : uses the given filename to open the file using xarray funtion
               2. data_obs  : uses obs[var] extract from postprocess stored in somisana toolkit; 
                              it retrieves obs data array
               3. time_obs  : uses the traditional dot method to retrieve time array from obs and corrects it from
                              a datetime64 datatype to a datetime.datetime object to ensure compliance with model time 
                              astype methos can be found at:
                              https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.astype.html
               4. lat_obs   : uses the dot method to extract latitude values. the same applieas to long_obs
               
            Parameters:
            - fname_obs         :filename of the observations
            - var               :variable input by the user. should be the same in both the modeland obs netCDFs.

            Returns:
            - time_obs, data_obs, long_obs, lat_obs
               
    """
    print("1. Im in get_ts_obs")
    # A generic function like this will work if we process all the observations
    obs = xr.open_dataset(fname_obs)
    data_obs = np.squeeze(obs[var].values)
    time_obs = obs.time.values
    time_obs = time_obs.astype('datetime64[s]').astype(datetime)
    long_obs = obs.longitude.values
    lat_obs = obs.latitude.values

    return time_obs, data_obs, long_obs, lat_obs

def obs_2_model_timeaxis(fname_obs,time_model,model_frequency,var):
    """
           The next Function is for matching time axis of both the insitu and model datasets:
               1. obs_2_model_timeaxis function is designed to put the observation array on the model time axis.
                   This is achieved by reading the obs dataset using 
                   xarray (https://docs.xarray.dev/en/stable/generated/xarray.open_dataset.html) and then 
                   resample that array by taking daily means at base=12, which means daily means of the hourly 
                   obs dataset at 12:00 midday. here is the resample function source:
                   https://www.geeksforgeeks.org/python-pandas-dataframe-resample/
               2. the operationof converting the datetime64 to datetime.datetime is performed to get the obs 
                   time array onto the same data type as the model time array.
               3. the data_obs_model_timeaxis is first created as an empty array of null values of length 
                   time_model. this empty list is populated by the for loop that appends to each index at each
                   model index (index_mod) based on each time of the model index (time_model_now)
               4. the if statement selects conditions at which the for loop applies which are: at every 
                   time step of the time model that also exists in the observations time array.
               5. when the above condition is met, a new index_obs is computed as the index when time obs and
                   time model are the same. That index output represents a location of the obervation that 
                   is to be placed on the model time axis on that time axis array. the x.index method is a 
                   list index extarctor method explained in https://www.w3schools.com/python/ref_list_index.asp
              6. the data_obs_model_timeaxis[index_mod] = data_obs[index_obs] simply assigns the observations to 
                  the above mentioned indices. therefore the data_obs_model_timeaxis is now observations at a model
                  time axis as it is named. 
                  
            Parameters:
            - fname_obs:        filename of observations
            - time_model:       numpy array or list, corresponding time values
            - model_frequency:  user input, it is the model frequency if the model is an average model output
            - var:              variable input by the user. should be the same in both the modeland obs netCDFs.

            Returns:
            - data_obs_model_timeaxis
    """
    print("4. Im in obs_2_new_timeaxis")
    obs = xr.open_dataset(fname_obs)
    obs_formatted = obs.resample(time=model_frequency, base=int(model_frequency[:-1])/2).mean()
    print("base=",int(model_frequency[:-1])/2)
    data_obs = np.squeeze(obs_formatted[var].values)
    time_obs = obs_formatted.time.values
    time_obs = time_obs.astype('datetime64[s]').astype(datetime)
    data_obs_model_timeaxis = [None for i in range(len(time_model))]
    # Convert the NumPy array to a list of datetime objects
    formatted_time_obs = time_obs.tolist()

    for index_mod, time_model_now in enumerate(time_model):
        # Step: Check the time component in time_obs and in each record if it is contained in time_obs then.
        if time_model_now in formatted_time_obs:            
        # Step: If contained then replace that time value in with a value in time_model in the same index.
            index_obs = formatted_time_obs.index(time_model_now)
            data_obs_model_timeaxis[index_mod] = data_obs[index_obs]
       
    return data_obs_model_timeaxis
  
# %% Statistical analysis section
    """
    Calculate model statistical elements: model_mean,obs_mean,model_min,obs_min,model_max,obs_max all 
    rounded to 3 decimal places

    Parameters:
    - model_data: numpy array or list, var values
    - obs data: numpy array or list, of observation values

    Returns:
    - model_mean,obs_mean,model_min,obs_min,model_max,obs_max
    """    

def calculate_rmse(obs_data, model_data, decimal = 3):
    # Calculate RMSE, ignoring NaN values
    rmse = round(np.sqrt(np.nanmean((obs_data - model_data)**2)),decimal)
    return rmse

def calculate_correlation(obs_data, model_data, decimal = 3):
    # Calculate correlation, ignoring NaN values. This function is not really working but I am on it.
    model_in_situ_corr = round(np.corrcoef(obs_data, model_data)[0, 1],decimal) 
    
    # insitu_no_nan = obs_data[~np.isnan(obs_data)]
    # model_no_nan = model_data[~np.isnan(obs_data)]    
    # model_in_situ_corr = round(np.corrcoef(np.array(insitu_no_nan),np.array(model_no_nan)),decimal)[1][0]
    
    return model_in_situ_corr

def calculate_std_dev(data, decimal = 3):
    # Calculate standard deviation, ignoring NaN values
    std_dev = round(np.nanstd(data),decimal)
    return std_dev

def calculate_total_bias(obs_data, model_data, decimal = 3):
    # Calculate total bias, ignoring NaN values
    bias = round(np.nanmean(obs_data - model_data),decimal)
    return bias

def calculate_min_mean_max(obs_data, model_data, decimal = 3):
    # Calculate model, obs mean,min and max ignoring NaN values
    model_mean = round(np.nanmean(model_data),decimal)
    obs_mean = round(np.nanmean(obs_data),decimal)
    model_min = round(np.min(model_data),decimal)
    obs_min = round(np.nanmin(obs_data),decimal)
    model_max = round(np.max(model_data),decimal)
    obs_max = round(np.nanmax(obs_data),decimal)
    
    return model_mean,obs_mean,model_min,obs_min,model_max,obs_max

def calculate_seasonal_means(model_data, time_model,var):
    """
    Calculate seasonal means (JFM, AMJ, JAS, OND) from variable and time arrays.

    Parameters:
    - model_data: numpy array or list, var values
    - time_model: numpy array or list, corresponding time values
    - var: variable input by the user. should be the same in both the modeland obs netCDFs.

    Returns:
    - seasonal means for JFM, AMJ, JAS, OND
    """
    # Convert time_model to datetime format
    time_model = pd.to_datetime(time_model)

    # Create a DataFrame with time and var columns
    df = pd.DataFrame({'Time': time_model, var: model_data})

    # Set the time column as the index
    df.set_index('Time', inplace=True)

    # Resample data to get monthly means
    monthly_mean = df.resample('M').mean()

    # Extract specific columns for each season
    JFM_mean = monthly_mean[monthly_mean.index.month.isin([1, 2, 3])][var].mean()
    AMJ_mean = monthly_mean[monthly_mean.index.month.isin([4, 5, 6])][var].mean()
    JAS_mean = monthly_mean[monthly_mean.index.month.isin([7, 8, 9])][var].mean()
    OND_mean = monthly_mean[monthly_mean.index.month.isin([10, 11, 12])][var].mean()

    return JFM_mean, AMJ_mean, JAS_mean, OND_mean


# %%
def get_model_obs_ts(fname,fname_obs,output_path,model_frequency,var,depth=-1,i_shifted=0,j_shifted=0,ref_date=None,lon_extract=None):
    """
           The next Function is for retrieving all datasets from the previous functions and package them into a netCDF file: 
               1. initialize the function by calling all the parametres that will apply to subsequent functions.
               2. Create a NetCDF file after getting returns from all the functions in the middle run
               
            Parameters:
            - fname             :filename of the model
            - fname_obs         :filename of observations
            - output_path       :filename of fname_out
            - var               :variable input by the user. should be the same in both the modeland obs netCDFs.
            - lat               :lat read from the insitu file
            - lon               :lon read from the insitu file 
            - ref_date          :static user input
            - depth             :user input based on insitu sensor depth
            - model_frequency   :user input, it is the model frequency if the model is an average model output
            - i_shifted         :optional user inputs if the z level of the model does not match the insitu depth
            - j_shifted         :optional user inputs if the z level of the model does not match the insitu depth
            - time_lims         :limits are computed based on the length of the insitu data that matches model span
    """
    print("5. Im in get_model_obs_ts")
    # the output of this function is a netcdf file 'output_path'
    # which will have the model and observations on the same time axis
    
    # get the observations time-series
    time_obs, data_obs, long_obs, lat_obs = get_ts_obs(fname_obs,var)   

    # get the model time-series
    time_model, data_model,lat_mod,lon_mod = post.get_ts(fname,var,long_obs,lat_obs,ref_date,depth=depth,i_shifted=i_shifted,j_shifted=j_shifted,time_lims=[time_obs[0],time_obs[-1]]) # Change the 10 back to -1

    # get the observations onto the model time axis
    data_obs_model_timeaxis = obs_2_model_timeaxis(fname_obs,time_model, model_frequency,var)

    # Create a NetCDF file
    with nc.Dataset(output_path, 'w', format='NETCDF4') as nc_file:
        # Create dimensions
        nc_file.createDimension('time', len(time_model))

        if j_shifted != 0:
            lat_mod = np.array([lat_mod], dtype=np.float32)
            nc_file.createDimension('latitude', len(lat_mod))
        else:
            nc_file.createDimension('latitude', len(lat_obs))

        if i_shifted != 0:
            nc_file.createDimension('longitude', len(lon_mod))
        else:
            nc_file.createDimension('longitude', len(long_obs))
        
        print(f"6.1 NetCDF file created at: {output_path}")
        # Create variables
        time_var = nc_file.createVariable('time', 'f8', ('time'))
        lat_var = nc_file.createVariable('latitude_insitu', 'f4', ('latitude'))
        lon_var = nc_file.createVariable('longitude_insitu', 'f4', ('longitude'))        
        lat_mod_var = nc_file.createVariable('latitude_on_model', 'f4', ('latitude'))
        lon_mod_var = nc_file.createVariable('longitude_on_model', 'f4', ('longitude'))        
        model_var = nc_file.createVariable('data_model', 'f4', ('time', 'latitude', 'longitude'))
        obs_model_var = nc_file.createVariable('data_obs_model_timeaxis', 'f4', ('time', 'latitude', 'longitude'))
    
        # Convert datetime objects to Unix timestamps (floats)
        float_time_model = np.array([dt.timestamp() for dt in time_model], dtype=int)
        
        # Convert each float timestamp to datetime
        for dt in time_model:
            print(dt)
        
        # Assign data to variables
        time_var[:] = float_time_model
        lat_var[:] = lat_obs
        lon_var[:] = long_obs        
        lat_mod_var[:] = lat_mod
        lon_mod_var[:] = lon_mod        
        model_var[:, :, :] = data_model
        obs_model_var[:, :, :] = data_obs_model_timeaxis
 
        # Add attributes if needed
        time_var.units = 'seconds since 1970-01-01 00:00:00'    
        time_var.calendar = 'standard'        
        time_var.long_name = 'time'
        # time_var.units = 'days'
        lat_var.units = 'latitude'
        lon_var.units = 'longitude'
        lat_mod_var.units = 'latitude' 
        lon_mod_var.units= 'longitude'
        model_var.units = 'degrees Celsius' # SHOULDN'T BE HARD CODED
        obs_model_var.units = 'degrees Celsius'
                
        # calculate model seasonal means:
        # seasonal_mean = calculate_seasonal_means(data_model, time_model,var)
        # nc_file.setncattr('seasonal_mean', seasonal_mean)
        JFM_mean, AMJ_mean, JAS_mean, OND_mean = calculate_seasonal_means(data_model, time_model,var)
        nc_file.setncattr('JFM', round(JFM_mean,3))
        nc_file.setncattr('AMJ', round(AMJ_mean,3))
        nc_file.setncattr('JAS', round(JAS_mean,3))
        nc_file.setncattr('OND', round(OND_mean,3))
            
        nc_file.setncattr('depth',depth)
        nc_file.setncattr('i_shift',i_shifted)
        nc_file.setncattr('j_shift',j_shifted)
        
        # Calculate and add correlations as attributes
        model_in_situ_corr = calculate_correlation(data_obs_model_timeaxis, data_model)
        nc_file.setncattr('correlation_model_obs', model_in_situ_corr)
        
        # Calculate and add standard deviations as attributes
        std_dev_model = calculate_std_dev(data_model)
        std_dev_obs_model = calculate_std_dev(data_obs_model_timeaxis)
        nc_file.setncattr('std_dev_model', std_dev_model)
        nc_file.setncattr('std_dev_obs_model', std_dev_obs_model)
        
        # Calculate RMSE and add it as an attribute
        rmse_model_obs = calculate_rmse(data_obs_model_timeaxis, data_model)
        nc_file.setncattr('rmse_model_obs', rmse_model_obs)
        
        # Calculate and add total bias as an attribute
        total_bias = calculate_total_bias(data_obs_model_timeaxis, data_model)
        nc_file.setncattr('total_bias', total_bias)
               
        # Calculate and add total bias as an attribute
        model_mean,obs_mean,model_min,obs_min,model_max,obs_max = calculate_min_mean_max(data_obs_model_timeaxis, data_model)
        nc_file.setncattr('model_mean', model_mean)
        nc_file.setncattr('obs_mean', obs_mean)
        nc_file.setncattr('model_min', model_min)
        nc_file.setncattr('obs_min', obs_min)
        nc_file.setncattr('model_max', model_max)
        nc_file.setncattr('obs_max', obs_max)
         
# %%
if __name__ == "__main__":
    # Define the input parameters
    dir_model = '/mnt/d/Run_False_Bay_2008_2018_SANHO/croco_avg_Y2016*.nc.1'
    fname_obs = '/mnt/d/DATA-20231010T133411Z-003/DATA/ATAP/Processed/Processed_Station_Files/CapePoint_CP002.nc'
    fname_out = 'Validation_'+'CapePoint_CP002.nc' #'CapePoint_CP002.nc'  'FalseBay_FB001.nc'

    # Output file name and directory
    output_directory = '/mnt/d/Run_False_Bay_2008_2018_SANHO/Validation/ATAP/model_validation/'
    fname_out = os.path.join(output_directory, fname_out)

    # Other parameters
    model_frequency='24H'
    var = 'temp'
    depth=-40
    ref_date = datetime(1990, 1, 1, 0, 0, 0)
    
    get_model_obs_ts(dir_model,fname_obs,
                      fname_out,model_frequency=model_frequency,
                      var=var,
                      ref_date=ref_date,
                      depth=depth, 
                      i_shifted=0,j_shifted=-2      
                      # ,lat_extract = -34.4        
                      )
    
