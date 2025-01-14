import xarray as xr
import cftime
import pandas as pd
import os,sys
from datetime import datetime, timedelta
import subprocess
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

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
    
    # Download section
    MAX_RETRIES, RETRY_WAIT = 3, 10
    variable=None
    for attempt in range(MAX_RETRIES):
        print('')
        print(f'Attempt {attempt+1} out of {MAX_RETRIES} tries for {metadata[var]["vars"][0]}.')
        try:
            print(f'Connecting to {metadata[var]["url"]} to subset and download {metadata[var]["vars"][0]}.')
            ds = xr.open_dataset(metadata[var]["url"],
                                 drop_variables=vars_to_drop,
                                 decode_times=False,
                                 engine="netcdf4").sel(lat=lat_range,
                                                       lon=lon_range)


            
            if 'time' in ds:
                ds['time'] = decode_time_units(ds['time'])
                tmax=pd.to_datetime(ds.time.max().values)
                timesteps_to_add=0
                if tmax<=end_date:
                    # Number of time steps requested
                    requested_num_timesteps = (end_date - start_date).days
                    # Number of timesteps available in the dataset
                    available_timesteps = (tmax - start_date).days
                    # Number of timesteps that needs to be added to have a "complete array"
                    timesteps_to_add = requested_num_timesteps - available_timesteps
                    # Subset the dataset temporally based on what is available
                    ds = ds.sel(time=slice(start_date, tmax.replace(hour=0, minute=0, second=0, microsecond=0)))
                else:
                    ds = ds.sel(time=time_range)

            variable = ds[metadata[var]["vars"][0]]
            
            if variable.ndim == 4:
                variable = variable.sel(depth=depth_range)
            
            if run_date.hour == 0:
                variable = variable.resample(time='1D').mean()
            elif run_date.hour == 12:
                variable = variable.resample(time='1D',offset='12h').mean()
            else:
                print(f'Strange run date: {run_date.hour}')
                print(f'Must be either 0 or 12.')
            
            # if the dataset has less timesteps than what was requested, the script makes duplicates of the last time step to fill the data array.
            if timesteps_to_add > 0:
                if run_date.hour == 12:
                    timesteps_to_add+=1
                else:
                    pass
                # Create additional time values
                last_time = variable.coords["time"][-1]
                time_diff = variable.coords["time"][-1] - variable.coords["time"][-2]
                new_times = [last_time + (i + 1) * time_diff for i in range(timesteps_to_add)]

                # Duplicate the last timestep data for each new time
                last_timestep = variable.isel(time=-1)
                new_data = xr.concat([last_timestep.expand_dims("time") for _ in range(timesteps_to_add)], dim="time")

                # Assign new time coordinates to the duplicated data
                timestamps = np.array([da.values for da in new_times]) 
                new_data = new_data.assign_coords(time=timestamps)

                # Combine the original data with the new data
                extended_data = xr.concat([variable, new_data], dim="time")

                # Clean up the extra arrays
                variable = extended_data
                
                del extended_data, last_time, time_diff, new_times, last_timestep, new_data, timestamps

            save_path = os.path.join(save_dir, f"hycom_{metadata[var]['vars'][0]}.nc")
            variable.to_netcdf(save_path, 'w')
            ds.close()
            
            if validate_download(save_path, metadata[var]["vars"][0], start_date, end_date):
                print('')
                print(f'File written to {save_path} and validation was successful.')
                break
            
            else:
                if attempt < MAX_RETRIES - 1:
                    print('')
                    print(f"File {save_path} validation failed and retrying the download")
                    print(f"Retrying in {RETRY_WAIT} seconds...")
                    time.sleep(RETRY_WAIT)
                else:
                    print(f"Failed to download after the maximum number of attempts: {MAX_RETRIES}.")

        except Exception as e:
            print(f"Error: {e}")
            if attempt < MAX_RETRIES - 1:
                print(f"Retrying in {RETRY_WAIT} seconds...")
                time.sleep(RETRY_WAIT)
            else:
                print(f"Failed to download after the maximum number of attempts: {MAX_RETRIES}.")

        finally:
            # Explicitly delete large variables to free up memory
            del variable

def download_vars_parallel(variables, domain, depths, run_date, hdays, fdays, workers, save_dir):
    var_metadata = update_var_list(variables)

    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_var = {
            executor.submit(download_var, var, var_metadata, domain, depths, save_dir, run_date, hdays, fdays): var
            for var, metadata in var_metadata.items()
        }

        for future in as_completed(future_to_var):
            var = future_to_var[future]
            try:
                future.result()
                return True
            except Exception as e:
                print(f"Download failed for {var}: {e}")
                return False
                
def download_hycom(variables, domain, depths, run_date, hdays, fdays, save_dir, workers=None):
    """
    Downloads the HYCOM analysis in daily outputs using xarrray opendap.
    This function does check the integrity of the file. If the file is corrupt the download retries - max retries is 3.

    INPUTS:
    variables: List of variables to download (e.g. ['salinity', 'water_temp', 'surface_el', 'water_u', 'water_v'])
    domain   : List of geographical coordinates to subset the data and download (e.g. [lon_min,lon_max,lat_min,lat_max]).
    depths   : List of minimum and maximum depths to download. Values must be positive (e.g. [0,5000]).
    run_date : Todays datetime to download (e.g. datetime.datetime(YYYY,MM,DD)).
    hdays    : Days to hindcast (e.g. hdays=5).
    fdays    : Days to forecast (e.g. fdays=5).
    save_dir : Directory to save the downloaded data (eg. save_dir='/path/and/directory/to/save/').
    workers  : It is the number of variables to download in parallel. Default is None, in which cases it downloads all the variables in paralell.
               To note, the number of workers cannot exceed the number of variables.

    OUTPUT:
    NetCDF file containing the most recent HYCOM forcast run.
    """
    # We add an additional day to ensure that it exceeds the model run time. 
    hdays, fdays = hdays+1, fdays+1
    
    # I am opting to download in series here. Its been hardcoded and will be changed. Here I am just trying it out eleviate stress on the computre server.
    # it should be workers=len(params)
    if workers is None:
        workers=np.size(variables)
    else:
        pass

    if download_vars_parallel(variables, domain, depths, run_date, hdays, fdays, workers, save_dir):
        ds = xr.open_mfdataset(os.path.join(save_dir, 'hycom_*.nc'))
        outfile = os.path.abspath(os.path.join(save_dir, f"HYCOM_{run_date.strftime('%Y%m%d_%H')}.nc"))
        if os.path.exists(outfile):
            os.remove(outfile)
        ds.to_netcdf(outfile, 'w')
        subprocess.call(["chmod", "-R", "775", outfile])
        print('')
        print('created: ', outfile)
        print('')
    else:
        print('HYCOM download failed.')

if __name__ == '__main__':
    run_date = pd.to_datetime('2025-01-14 12:00:00')
    hdays = 2
    fdays = 7
    variables = ['salinity','water_temp','surf_el','water_u','water_v']
    domain = [23,24,-37,-36]
    depths = [0,5]
    save_dir = '/home/g.rautenbach/Projects/somisana-croco/DATASETS_CROCOTOOLS/HYCOM/'
    download_hycom(variables, domain, depths, run_date, hdays, fdays, save_dir)
