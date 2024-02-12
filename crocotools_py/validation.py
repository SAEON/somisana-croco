import numpy as np
import xarray as xr 
from glob import glob
from datetime import datetime
import crocotools_py.postprocess as post
import netCDF4 as nc
import pandas as pd
import os
from colorama import Fore, Style

def get_ts_obs(fname_obs, var):
    """
           The next Function is for loading IN SITU time series:

            Parameters:
            - fname_obs         :filename of the observations
            - var               :variable input by the user. should be the same in both the modeland obs netCDFs.

            Returns:
            - time_obs, data_obs, long_obs, lat_obs               
    """
    # print("1. Im in get_ts_obs")
    # obs: uses the given filename to open the file using xarray funtion
    obs = xr.open_dataset(fname_obs)
    # data_obs: uses obs[var] it retrieves obs data array
    data_obs = np.squeeze(obs[var].values)
    var_units = obs[var].units
    time_obs = obs.time.values
    # time_obs  : uses the traditional dot method to retrieve time array from obs and corrects it from
    # a datetime64 datatype to a datetime.datetime object to ensure compliance with model time 
    # astype methos can be found at: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.astype.html
    time_obs = time_obs.astype('datetime64[s]').astype(datetime)
    # lat_obs: uses the dot method to extract latitude values. the same applieas to long_obs
    #But First it must handle an exception for an example, if the obs data defines lats and lons in full names    
    try:
        long_obs = obs.longitude.values
        lat_obs = obs.latitude.values
    
        if (long_obs == [0]).any():
            raise ValueError("Undesired output detected in 'longitude': array([0])")
        if (lat_obs == [0]).any():
            raise ValueError("Undesired output detected in 'latitude': array([0])")
        
    except ValueError as ve:
        print(f"Caught an exception: {ve}")
        # Handle the exception, provide default values
        long_obs = obs.lon.values
        lat_obs = obs.lat.values
    
    else:
        print("Longitude and latitude successfully retrieved.")

    return time_obs, data_obs, long_obs, lat_obs, var_units

def obs_2_model_timeaxis(fname_obs,time_model,model_frequency,var):
    """
           This Function is for matching time axis of both the insitu and model datasets:
            obs_2_model_timeaxis function is designed to put the observation array on the model time axis.
            This is achieved by reading the obs dataset using 
              
            Parameters:
            - fname_obs:        filename of observations
            - time_model:       numpy array or list, corresponding time values
            - model_frequency:  user input, it is the model frequency if the model is an average model output
            - var:              variable input by the user. should be the same in both the modeland obs netCDFs.

            Returns:
            - data_obs_model_timeaxis
    """
    # print("4. Im in obs_2_new_timeaxis")
    # xarray (https://docs.xarray.dev/en/stable/generated/xarray.open_dataset.html) imported as xr 
    obs = xr.open_dataset(fname_obs)
    # resample that array by taking daily means at base=12, which means daily means of the hourly 
    # obs dataset at 12:00 midday. here is the resample function source:
    # https://www.geeksforgeeks.org/python-pandas-dataframe-resample/
    obs_formatted = obs.resample(time=model_frequency, base=int(model_frequency[:-1])/2).mean()
    #print("base=",int(model_frequency[:-1])/2)
    data_obs = np.squeeze(obs_formatted[var].values)
    time_obs = obs_formatted.time.values
    # the operationof converting the datetime64 to datetime.datetime is performed to get the obs 
    # time array onto the same data type as the model time array.
    time_obs = time_obs.astype('datetime64[s]').astype(datetime)
    # the data_obs_model_timeaxis is first created as an empty array of null values of length 
    # time_model. this empty list is populated by the for loop that appends to each index at each
    # model index (index_mod) based on each time of the model index (time_model_now)
    data_obs_model_timeaxis = [None for i in range(len(time_model))]
    # Convert the NumPy array to a list of datetime objects
    formatted_time_obs = time_obs.tolist()
    
    if var == 'u' or 'v':
        for index_mod, time_model_now in enumerate(time_model):
            if time_model_now in formatted_time_obs:            
                index_obs = formatted_time_obs.index(time_model_now)
                data_obs_model_timeaxis[index_mod] = data_obs[index_obs]
    else:
        for index_mod, time_model_now in enumerate(time_model):
            # the if statement selects conditions at which the for loop applies which are: at every 
            # time step of the time model that also exists in the observations time array.
            # Step through index: Checks the time component in time_obs and in each record if it is contained in time_obs then.
            if time_model_now in formatted_time_obs:            
            # Step: If contained then replace that time value in with a value in time_model in the same index.
            # when the above condition is met, a new index_obs is computed as the index when time obs and
            # time model are identical times. That index output represents a location of the obervation that 
            # is to be placed on the model time axis on that time axis array. the x.index method is a in https://www.w3schools.com/python/ref_list_index.asp
                
                index_obs = formatted_time_obs.index(time_model_now)
                data_obs_model_timeaxis[index_mod] = data_obs[index_obs]
       
    return np.array(data_obs_model_timeaxis)
  
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
            - model_frequency   :user input, it is the model frequency if the model is an average model output
            - i_shifted         :optional user inputs if the z level of the model does not match the input lon depth
            - j_shifted         :optional user inputs if the z level of the model does not match the input lat depth
            - time_lims         :limits are computed based on the length of the insitu data that matches model span
    """
    # print("5. Im in get_model_obs_ts")
    # the output of this function is a netcdf file 'output_path'
    # which will have the model and observations on the same time axis
    
    # get the observations time-series
    time_obs, data_obs, long_obs, lat_obs, var_units = get_ts_obs(fname_obs,var)   

    # get the model time-series
    time_model, data_model,lat_mod,lon_mod,h = post.get_ts(fname,var,long_obs,lat_obs,ref_date,depth=depth,i_shifted=i_shifted,j_shifted=j_shifted,time_lims=[time_obs[0],time_obs[-1]]) # Change the 10 back to -1

    # get the observations onto the model time axis
    data_obs_model_timeaxis = obs_2_model_timeaxis(fname_obs,time_model, model_frequency,var)
    
    insitu_no_nan =  data_obs_model_timeaxis[~np.isnan( data_obs_model_timeaxis)]  
    model_no_nan = data_model[~np.isnan(data_obs_model_timeaxis)]  

    # Create a NetCDF file
    with nc.Dataset(output_path, 'w', format='NETCDF4') as nc_file:
        # Create dimensions
        nc_file.createDimension('time', len(time_model))

        if j_shifted != 0:
            lat_mod = np.array([lat_mod], dtype=np.float32)
            nc_file.createDimension('latitude', len(lat_mod))
        else:
            lat_obs = np.array([lat_obs], dtype=np.float32)
            nc_file.createDimension('latitude', len(lat_obs))

        if i_shifted != 0:
            lon_mod = np.array([lon_mod], dtype=np.float32)
            nc_file.createDimension('longitude', len(lon_mod))
        else:
            long_obs = np.array([long_obs], dtype=np.float32)
            nc_file.createDimension('longitude', len(long_obs))
        
        # print(f"5.1 NetCDF file created at: {output_path}")
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
        # for dt in time_model:
        #     print(dt)
        
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
        lat_var.units = 'degrees_north'
        lon_var.units = 'degrees_east'
        lat_mod_var.units = 'degrees_north'
        lon_mod_var.units= 'degrees_east'
        model_var.units = var_units     # SHOULDN'T BE HARD CODED
        obs_model_var.units = var_units
        
        
                
        # calculate model seasonal means:
        # seasonal_mean = calculate_seasonal_means(data_model, time_model,var)
        # nc_file.setncattr('seasonal_mean', seasonal_mean)
        JFM_mean, AMJ_mean, JAS_mean, OND_mean = calculate_seasonal_means(data_model, time_model,var)
        nc_file.setncattr('JFM', round(JFM_mean,3))
        nc_file.setncattr('AMJ', round(AMJ_mean,3))
        nc_file.setncattr('JAS', round(JAS_mean,3))
        nc_file.setncattr('OND', round(OND_mean,3))
         
        nc_file.setncattr('depth',depth)
        nc_file.setncattr('h',h)
        nc_file.setncattr('i_shift',i_shifted)
        nc_file.setncattr('j_shift',j_shifted)
        
        decimal = 3
        
        # Calculate and add correlations as attributes
        model_in_situ_corr = np.corrcoef(np.array(insitu_no_nan),np.array(model_no_nan))[1][0]
        model_in_situ_corr = round(model_in_situ_corr, decimal)
        nc_file.setncattr('correlation_model_obs', model_in_situ_corr)
        
        # Calculate and add standard deviations as attributes
        std_dev_model = round(np.std(model_no_nan),decimal)
        std_dev_obs_model = round(np.std(insitu_no_nan),decimal)        
        nc_file.setncattr('std_dev_model', std_dev_model)
        nc_file.setncattr('std_dev_obs_model', std_dev_obs_model)
        
        # Calculate RMSE and add it as an attribute
        rmse_model_obs = round(np.sqrt(np.mean((insitu_no_nan - model_no_nan)**2)),decimal)
        nc_file.setncattr('rmse_model_obs', rmse_model_obs)
        
        # Calculate and add total bias as an attribute
        total_bias = round(np.mean(insitu_no_nan - model_no_nan),decimal)
        nc_file.setncattr('total_bias', total_bias)
        
        # Calculate model, obs mean,min and max ignoring NaN values
        model_mean = round(np.mean(model_no_nan),decimal)
        obs_mean = round(np.mean(insitu_no_nan),decimal)
        model_min = round(np.min(model_no_nan),decimal)
        obs_min = round(np.min(insitu_no_nan),decimal)
        model_max = round(np.max(model_no_nan),decimal)
        obs_max = round(np.max(insitu_no_nan),decimal)
                       
        # Write total bias as an attribute
        nc_file.setncattr('model_mean', model_mean)
        nc_file.setncattr('obs_mean', obs_mean)
        nc_file.setncattr('model_min', model_min)
        nc_file.setncattr('obs_min', obs_min)
        nc_file.setncattr('model_max', model_max)
        nc_file.setncattr('obs_max', obs_max)
        
        # Console feedback to the user
        if -depth<h:
            print('')
            print('----------------------------------------------------------------')
            print(f'{Fore.GREEN}{Style.BRIGHT} Successfully Done!!!{Style.RESET_ALL}')
            print('----------------------------------------------------------------')
            print(f'Find your evaluation output netCDF file at: {output_directory}')
        elif i_shifted == 0 and j_shifted == 0 and h<-depth:
            print('')
            print('----------------------------------------------------------------')
            print(f'{Fore.GREEN}{Style.BRIGHT} Successfully Done!!!{Style.RESET_ALL}')
            print('----------------------------------------------------------------')
            print(f'{Fore.YELLOW}{Style.BRIGHT} WARNING: {Style.RESET_ALL}')
            print(f'{Fore.CYAN}{Style.BRIGHT} The following issue is caused by lower resolution bathymetry data used in the model, not very well suited for coastal areas* {Style.RESET_ALL}')
            print(f'{Fore.YELLOW}{Style.BRIGHT} The model is shallower than the in situ data due to bathymetric inconsistencies, your model height is (h = {-round(h)} m), {Style.RESET_ALL}')
            print(f'{Fore.YELLOW}{Style.BRIGHT} while in situ data is at (depth = {depth} m). Your model output is therefore taken at the deepest model point since model is shallower.  {Style.RESET_ALL}')
            print(f'{Fore.YELLOW}{Style.BRIGHT} Alternatively; if prefered, the user can consider changing/shifting i_shift, j_shift or both from zero to deeper parts of the ocean {Style.RESET_ALL}')
            print(f'{Fore.YELLOW}{Style.BRIGHT} REMEMBER:  {Style.RESET_ALL}')
            print(f'{Fore.YELLOW}{Style.BRIGHT} i_shifted :optional user input for shifting input lon/xi by a grid point at a time on function get_model_obs_ts() {Style.RESET_ALL}')
            print(f'{Fore.YELLOW}{Style.BRIGHT} j_shifted :optional user input for shifting input lat/eta by a grid point at a time on function get_model_obs_ts() {Style.RESET_ALL}')
            print('')
            print(f'Find your evaluation output netCDF file at: {output_directory}')            
        else:
            print('')
            print('----------------------------------------------------------------')
            print(f'{Fore.RED}{Style.BRIGHT} Evaluation FAILED!!!{Style.RESET_ALL}')
            print('----------------------------------------------------------------')
            print(f"{Fore.MAGENTA}{Style.BRIGHT} But don't panic, here's how you can fix it: {Style.RESET_ALL}")
            print('Consider changing i_shift, j_shift or both from zero to deeper parts of the ocean or change the station')
            print(f'the given input is of observation depth ={Fore.GREEN}{Style.BRIGHT} {-depth} m {Style.RESET_ALL}> h = -{round(h)} m, this means you are evaluating below sea-floor; remember h is model sea-height')
            print('')
            print('Remember:')
            print('i_shifted :optional user inputs useful for shifting input lon or specifically xi by a grid point at a time')
            print('j_shifted :optional user inputs useful for shifting input lat or specifically eta by a grid point at a time')
            
               
# %%
    
if __name__ == "__main__":

    # Define the input parameters:  
    # fname_out is the file name you want your output netCDF file to be labelled as.
    fname_out = 'CapePoint_CP002.nc' # or 'CapePoint_CP003.nc'  'FalseBay_FB001.nc'
    # dir_model is the directory at which your model files are located.
    dir_model = '/mnt/d/Run_False_Bay_2008_2018_SANHO/croco_avg_Y201*.nc.1'

    # fname_obs is the file name of your observations
    fname_obs = f'/mnt/d/DATA-20231010T133411Z-003/DATA/ATAP/Processed/Processed_Station_Files/{fname_out}'

    # Output file name and directory:
    # Be sure to change the following directory to the one in which you are working
    output_directory = '/mnt/d/Run_False_Bay_2008_2018_SANHO/Validation/ATAP/model_validation/'#'i-j_Shifted/'
    fname_out = os.path.join(output_directory, 'Validation_'+fname_out )

    # Other parameters:
    # The model frequency is the rate at which the model produces average outputs. usualy an attribute in the model
    model_frequency='24H'
    # The variable of interest during validation. Should idealy be extracted accurately from the observations file
    var = 'temp'
    # The depth you are providing should be the in situ station depth
    depth=-50
    # The following ref_date corrasponds to the refdate of croco SWCC model, 
    # change it if you model has a different ref date
    ref_date = datetime(1990, 1, 1, 0, 0, 0)
    
    get_model_obs_ts(dir_model,fname_obs,
                      fname_out,model_frequency=model_frequency,
                      var=var,
                      ref_date=ref_date,
                      depth=depth,
                      # The following are optional user inputs. Should be left as 0 on both, otherwise if the output
                      # netCDF file returns empty arrays, then what has happened is that the model z level is 
                      # modelled to be beneath bathymetry for the above statn depth. in this case the 
                      #i_shifted and j_shifted can be used to machanically shift the obs station to a deeper ocean.
                       i_shifted=0,j_shifted=0   #CapePoint_CP001.nc , depth=-32
                       # i_shifted=0,j_shifted=-3  #CapePoint_CP002.nc, depth=-50
                       # i_shifted=3,j_shifted=-3  #CapePoint_CP003.nc, depth=-58
                      # i_shifted=0,j_shifted=0   #FalseBay_FB001.nc , depth=-40
                       # i_shifted=-1,j_shifted=0  #FalseBay_FB002.nc, depth=-50
                       # i_shifted=-2,j_shifted=0  #FalseBay_FB003.nc, depth=-58
                      # i_shifted=0,j_shifted=0  #WalkerBay_WB003.nc, depth=-23
                      )
    
