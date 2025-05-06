import xarray as xr
import cftime
import pandas as pd
import os
from datetime import timedelta
import numpy as np
import requests
from glob import glob

def is_server_reachable(url):
    try:
        response = requests.get(url, timeout=5, stream=True)  # Get headers only
        return response.status_code == 200
    except requests.RequestException:
        return False

def is_file_valid(file_path):
    """
    Function to check that the file exists and that its not empty.
    """
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        return True
    else:
        print(f"File {file_path} does not exist or is empty.")
        return False

def is_netcdf_valid(file_path):
    """
    Function to check that the netCDF file can open.
    """
    try:
        with xr.open_dataset(file_path) as ds:
            return True
    except Exception as e:
        print(f"Error opening file {file_path}: {e}")
        return False

def check_variables(file_path, expected_vars):
    """
    Function to check that the netCDF file has the expected variable.
    """
    try:
        with xr.open_dataset(file_path) as ds:
            if expected_vars not in ds.variables:
                print(f"Variable {expected_vars} is not in {file_path}.")
                return False
            else:
                return True
    except Exception as e:
        print(f"Error checking variable in file {file_path}: {e}")
        return False

def check_time_range(file_path, expected_start, expected_end):
    """
    Function to check the time range of the file.
    """
    try:
        with xr.open_dataset(file_path) as ds:
            if 'time' in ds.coords:
                file_times = pd.DatetimeIndex(ds['time'].values)
                if (file_times.min() == expected_start) and (file_times.max() == expected_end):
                    return True
                else:
                    print(f"Time range in {file_path} is invalid.")
                    print(f"Expected time range: {expected_start} to {expected_end}")
                    print(f"File time range: {file_times.min()} to {file_times.max()}")
                    return False
            else:
                print(f"No 'time' coordinate found in {file_path}.")
                return False
    except Exception as e:
        print(f"Error checking time range in file {file_path}: {e}")
        return False


def check_data_quality(file_path, variables):
    """
    Function to check the quality of the data downloaded.
    """
    try:
        with xr.open_dataset(file_path) as ds:
            if variables in ds:
                data = ds[variables].values
                if data.size == 0:
                    print(f"Variable {var} has no data in {file_path}.")
                    return False
                else:
                    pass

                if np.ndim(data) == 4:
                    for dt in range(data[:,0,0,0].size):
                        if np.isnan(data[dt]).all():
                            print(f"All values for variable {var} at timestep {dt} are NaN in {file_path}.")
                            return False
                        else:
                            pass

                elif np.ndim(data) == 3:
                    for dt in range(data[:,0,0].size):
                        if np.isnan(data[dt]).all():
                            print(f"All values for variable {var} at timestep {dt} are NaN in {file_path}.")
                            return False
                        else:
                            pass

                else:
                    print(f"Variable {var} at timestep {dt} has incorrect dimension shape in {file_path}.")
                    return False

            else:
                print(f"Variable {var} not found in {file_path}.")
                return False

        return True

    except Exception as e:
        print(f"Error checking data quality in file {file_path}: {e}")

def validate_download(file_path, expected_vars, expected_start, expected_end):
    """
    Function that calls the validation routines to check the intergrity of the files downloaded.
    """
    if not is_file_valid(file_path):
        return False
    if not is_netcdf_valid(file_path):
        return False
    if not check_variables(file_path, expected_vars):
        return False
    if not check_time_range(file_path, expected_start, expected_end):
        return False
    if not check_data_quality(file_path, expected_vars):
        return False
    return True

def update_var_list(var_list):
    var_metadata = {
        'salinity': {
            "vars": ["salinity"],
            "url": "http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_s3z/FMRC_ESPC-D-V02_s3z_best.ncd",
        },
        'water_temp': {
            "vars": ["water_temp"],
            "url": "http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_t3z/FMRC_ESPC-D-V02_t3z_best.ncd",
        },
        'surf_el': {
            "vars": ["surf_el"],
            "url": "http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_ssh/FMRC_ESPC-D-V02_ssh_best.ncd",
        },
        'water_u': {
            "vars": ["water_u"],
            "url": "http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_u3z/FMRC_ESPC-D-V02_u3z_best.ncd",
        },
        'water_v': {
            "vars": ["water_v"],
            "url": "http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_v3z/FMRC_ESPC-D-V02_v3z_best.ncd",
        }
    }
    return {var: var_metadata[var] for var in var_list if var in var_metadata}

def decode_time_units(time_var):
    try:
        units = time_var.units
        calendar = getattr(time_var, 'calendar', 'standard')
        time = cftime.num2date(
            time_var[:], units=units, calendar=calendar,
            only_use_cftime_datetimes=False, only_use_python_datetimes=True
        )
        return pd.DatetimeIndex(time)
    except AttributeError as e:
        raise ValueError(f"Missing expected attributes in time_var: {e}")
    except Exception as e:
        raise RuntimeError(f"Error decoding time units: {e}")

def download_var(var, metadata, domain, depths, save_dir, run_date, hdays, fdays):
    vars_to_drop = ['salinity_bottom', 'water_temp_bottom', 'water_u_bottom', 'water_v_bottom', 'tau', 'time_offset',
                    'time_run', 'time1_offset', 'sst', 'sss', 'ssu', 'ssv', 'sic', 'sih', 'siu', 'siv', 'surtx',
                    'surty', 'steric_ssh']
    
    lon_range, lat_range, depth_range = slice(domain[0], domain[1]), slice(domain[2], domain[3]), slice(depths[0], depths[1])

    # Calculate time range for subsetting
    start_date = pd.Timestamp(run_date) - timedelta(days=hdays)
    end_date = pd.Timestamp(run_date) + timedelta(days=fdays)
    time_range = slice(start_date, end_date)
    
    variable=None
    
    try:
        print('')
        print(f'Connecting to {metadata[var]["url"]} to subset and download {metadata[var]["vars"][0]}.')
        ds = xr.open_dataset(metadata[var]["url"],
                             drop_variables=vars_to_drop,
                             decode_times=False,
                             engine="netcdf4").sel(lat=lat_range,
                                                   lon=lon_range)
                                                   
        if 'time' in ds: ds['time'] = decode_time_units(ds['time'])
        
        variable = ds[metadata[var]["vars"][0]]
        
        if variable.ndim == 4: variable = variable.sel(depth=depth_range)
        
        if run_date.hour == 0: 
            variable = variable.resample(time='1D').mean()
        elif run_date.hour == 12: 
            variable = variable.resample(time='1D',offset='12h').mean()
        else: 
            print(f'Invalid run date: {run_date.hour}')
        
        variable = variable.sel(time=time_range)
        
        save_path = os.path.join(save_dir, f"hycom_{metadata[var]['vars'][0]}.nc")
        variable.to_netcdf(save_path, 'w')
        ds.close()
        
        if validate_download(save_path, metadata[var]["vars"][0], start_date, end_date):
            print(f'File written to {save_path} and validation was successful.')        
        else:
            print(f"File {save_path} validation failed and retrying the download")

    except Exception as e:
        print(f"Error: {e}")

    finally:
        # Explicitly delete large variables to free up memory
        del variable

def pad_time_step(ds: xr.Dataset) -> xr.Dataset:
    """
    Pad all time-dependent variables in an xarray.Dataset by one timestep at the end,
    duplicating the last time step's data and inferring the time delta.
    """
    # Infer time and padded time
    time = ds['time']
    last_time = pd.to_datetime(time[-1].values)
    time_diff = pd.to_timedelta(time.diff('time').mean().values)
    padded_time = np.datetime64(last_time + time_diff, 'ns')
    new_time = np.append(time.values, padded_time)

    # Pad all time-dependent variables
    padded_vars = {}
    for var_name, da in ds.data_vars.items():
        if 'time' in da.dims:
            new_slice = da.isel(time=-1).expand_dims(time=[padded_time])
            new_slice[:] = da.isel(time=-1)
            da_padded = xr.concat([da, new_slice], dim='time')
            da_padded = da_padded.assign_coords(time=('time', new_time))
            padded_vars[var_name] = da_padded
        else:
            padded_vars[var_name] = da

    # Rebuild dataset excluding old time coord
    ds_padded = xr.Dataset(padded_vars, coords={k: v for k, v in ds.coords.items() if k != 'time'})
    
    # Assign the new time coordinate
    ds_padded = ds_padded.assign_coords(time=('time', new_time))
    
    return ds_padded

def download_hycom(variables, domain, depths, run_date, hdays, fdays, save_dir,pad=False):
    """
    Downloads the HYCOM analysis in daily outputs using xarrray opendap.
    This function does check the integrity of the file/s. 

    INPUTS:
    variables: List of variables to download (e.g. ['salinity', 'water_temp', 'surface_el', 'water_u', 'water_v'])
    domain   : List of geographical coordinates to subset the data and download (e.g. [lon_min,lon_max,lat_min,lat_max]).
    depths   : List of minimum and maximum depths to download. Values must be positive (e.g. [0,5000]).
    run_date : Todays datetime to download (e.g. datetime.datetime(YYYY,MM,DD)).
    hdays    : Days to hindcast (e.g. hdays=5).
    fdays    : Days to forecast (e.g. fdays=5).
    save_dir : Directory to save the downloaded data (eg. save_dir='/path/and/directory/to/save/').
    pad      : Pad all time-dependent variables in the dataset by one timestep at the start and end. 
               At the start, we download and extra day and at the end we copy the last timestep (Default is False).
               This is used operationally for our forecast models. 

    OUTPUT:
    NetCDF file containing the most recent HYCOM forcast run.
    """
    # We add an additional day to ensure that it exceeds the model run time. 
    # We also pad the dataset at the end, but instead of downloading it, we copy the last timestep.
    if pad: hdays = hdays + 1
    start_date = pd.Timestamp(run_date) - timedelta(days=hdays)
    end_date = pd.Timestamp(run_date) + timedelta(days=fdays)
    
    # This block of code is possibly not necessary, but its a nice to have feature.
    # It looks at the requested variables and the files in the save directory. 
    # If some of the variables are already there, they are removed from the list of the requested variables to download. 
    # This avoids having to redownload files that already exists. 
    # The integrity of these files are also checked to ensure they are good, otherwise they are redownloaded.
    files = glob(os.path.join(save_dir, 'hycom_*.nc'))
    # Extract variable names from the filenames
    existing_vars = [os.path.basename(f).replace('hycom_', '').replace('.nc', '') for f in files]
    print('')
    print(f'Requested variables that has already been downloaded: {existing_vars}')
    # Find missing variables
    missing_vars = [var for var in variables if var not in existing_vars]
    for file_path in files:
        var_name = os.path.basename(file_path).replace('hycom_', '').replace('.nc', '')
        if var_name in variables:
            try:
                # Run validation
                if not validate_download(file_path, var_name, start_date, end_date):
                    print('')
                    print(f"Validation failed for {file_path}. Marking as missing.") 
                    os.remove(file_path)
                    missing_vars.append(var_name)
            except Exception as e:
                print(f"Error validating {file_path}: {e}. Marking as missing.")
                missing_vars.append(var_name)
    # Print missing variables list
    print(f'Variables that still needs to be downloaded: {missing_vars}')
    
    # This function creates a metadata dictionary which comtains information about the variables. 
    # we are intersted in downloading
    var_metadata = update_var_list(missing_vars)
    
    # Testing the connection to the Thredds server where the HYCOM analysis is stored.
    server_url = "http://tds.hycom.org/thredds/dodsC/"
    if is_server_reachable(server_url):
        try:
            # We loop through the variable list and download using the download function
            for var in missing_vars:
                download_var(var, var_metadata, domain, depths, save_dir, run_date, hdays, fdays)
            # We count the number of files in the save directory and if it mathces the number of 
            # variables that was specified, we constructy the final HYCOM file. 
            nvar_files = len([os.path.basename(f).replace('hycom_', '').replace('.nc', '') for f in glob(os.path.join(save_dir, 'hycom_*.nc'))])
            if nvar_files==len(variables):
                # Run validation again for all files to ensure they are the correct.
                files = glob(os.path.join(save_dir, 'hycom_*.nc'))
                files_success = True
                for file_path in files:
                    var_name = os.path.basename(file_path).replace('hycom_', '').replace('.nc', '')
                    if var_name in variables:
                        try:
                            if not validate_download(file_path, var_name, start_date, end_date):
                                print('')
                                print(f"Validation failed for {file_path}. Removing file.")
                                files_success = False
                                os.remove(file_path)
                        except Exception as e:
                            print(f"Error validating {file_path}: {e}.")
                
                # If all the files in the variable list passed validation then 
                # we construct the final merged HYCOM file
                if files_success:
                    print('')
                    print("All variables are present and passed the validation test.")
                    print("Creating the combined netCDF file...")
                    ds = xr.open_mfdataset(os.path.join(save_dir, 'hycom_*.nc'))
                    if pad: ds = pad_time_step(ds)
                    outfile = os.path.abspath(os.path.join(save_dir, f"HYCOM_{run_date.strftime('%Y%m%d_%H')}.nc"))
                    if os.path.exists(outfile): os.remove(outfile)
                    ds.to_netcdf(outfile, 'w')
                    os.chmod(outfile, 0o775)
                    print('')
                    print('created: ', outfile)

        except Exception as e:
            print('')
            print(f"Failed to download dataset: {e}")
    else:
        print('')
        print(f"Server {server_url} is not reachable.")

if __name__ == '__main__':
    run_date = pd.to_datetime('2025-05-06 00:00:00')
    hdays = 5
    fdays = 5
    variables = ['salinity','water_temp','surf_el','water_u','water_v']
    domain = [23,24,-37,-36]
    depths = [0,5]
    save_dir = '/home/g.rautenbach/Projects/somisana-croco/DATASETS_CROCOTOOLS/HYCOM/'
    download_hycom(variables, domain, depths, run_date, hdays, fdays, save_dir, pad=True)
