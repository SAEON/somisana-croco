import xarray as xr
import cftime
import pandas as pd
import os
from datetime import datetime, timedelta
import subprocess
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

def check_time_range(start_date, end_date, time_coords):
    if (time_coords[0] <= start_date) & (end_date <= time_coords[-1]):
        print("Time range is within bounds.")
    else:
        print("Warning: Specified time range exceeds dataset range.")

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

def download_var(var, metadata, domain, depths, save_dir):
    vars_to_drop = ['salinity_bottom', 'water_temp_bottom', 'water_u_bottom', 'water_v_bottom', 'tau', 'time_offset',
                    'time_run', 'time1_offset', 'sst', 'sss', 'ssu', 'ssv', 'sic', 'sih', 'siu', 'siv', 'surtx', 
                    'surty', 'steric_ssh']
    lon_range, lat_range, depth_range = slice(domain[0], domain[1]), slice(domain[2], domain[3]), slice(depths[0], depths[1])
    
    MAX_RETRIES, RETRY_WAIT = 3, 20
    variable=None
    for attempt in range(MAX_RETRIES):
        print('')
        print(f'Attempt {attempt} out of {MAX_RETRIES} tries for {metadata["vars"][0]}.')
        try:
            print(f'Connecting to {metadata["url"]} to subset and download {metadata["vars"][0]}.')
            ds = xr.open_dataset(metadata["url"], 
                                 drop_variables=vars_to_drop, 
                                 decode_times=False, 
                                 engine="netcdf4").sel(lat=lat_range, 
                                                       lon=lon_range)
            if 'time' in ds:
                ds['time'] = decode_time_units(ds['time'])                
            
            variable = ds[metadata["vars"][0]]
            if variable.ndim == 4:
                variable = variable.sel(depth=depth_range)
            
            #variable = variable.resample(time='1D').mean()
            variable = variable.resample(time='1D', offset='12h').mean()
            save_path = os.path.join(save_dir, f"hycom_{metadata['vars'][0]}.nc")
            variable.to_netcdf(save_path, 'w')
            print(f'File written to {save_path}')
            ds.close()
            break
            
        except Exception as e:
            print(f"Error: {e}")
            if attempt < MAX_RETRIES - 1:
                print(f"Retrying in {RETRY_WAIT} seconds...")
                time.sleep(RETRY_WAIT)
            else:
                print("Failed to download after multiple attempts.")
                
        finally:
            # Explicitly delete large variables to free up memory
            del variable

def download_vars_parallel(variables, domain, depths, workers, save_dir):
    var_metadata = update_var_list(variables)
    
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_var = {
            executor.submit(download_var, var, metadata, domain, depths, save_dir): var 
            for var, metadata in var_metadata.items()
        }
        
        for future in as_completed(future_to_var):
            var = future_to_var[future]
            try:
                future.result()
            except Exception as e:
                print(f"Download failed for {var}: {e}")

def download_hycom(variables, domain, depths, run_date, hdays, fdays, save_dir,workers=None):
    """
    Downloads the most recent HYCOM analysis in daily outputs using xarrray opendap.

    INPUTS:
    variables: List of variables to download (e.g. ['salinity', 'water_temp', 'surface_el', 'water_u', 'water_v'])
    domain   : List of geographical coordinates to subset the data and download (e.g. [lon_min,lon_max,lat_min,lat_max]).
    depths   : List of minimum and maximum depths to download. Values must be positive (e.g. [0,5000]).
    run_date : Todays datetime to check that the data fall within the range (e.g. datetime.datetime(YYYY,MM,DD)).
    hdays    : Days to hindcast (e.g. hdays=5).
    fdays    : Days to forecast (e.g. fdays=5).
    save_dir : Directory to save the downloaded data (eg. save_dir='/path/and/directory/to/save/').
    workers  : It is the number of variables to download in parallel. Default is None, in which cases it downloads all the variables in paralell. To note, the number of workers cannot exceed the number of variables.

    OUTPUT:
    NetCDF file containing the most recent HYCOM forcast run.
    """    
    if workers is None:
        workers=len(variables)
    else:
        pass
    download_vars_parallel(variables, domain, depths, workers, save_dir) 
    start_date = pd.Timestamp(run_date) - timedelta(days=hdays)
    end_date = pd.Timestamp(run_date) + timedelta(days=fdays)
    ds = xr.open_mfdataset(os.path.join(save_dir, 'hycom_*.nc'))
    check_time_range(start_date, end_date, ds.coords['time'].values)
    outfile = os.path.abspath(os.path.join(save_dir, f"HYCOM_{run_date.strftime('%Y%m%d_%H')}.nc"))
    if os.path.exists(outfile):
        os.remove(outfile)
    ds.to_netcdf(outfile, 'w')
    subprocess.call(["chmod", "-R", "775", outfile])
    print('')
    print('created: ', outfile)
    print('')
