import xarray as xr
import cftime
import pandas as pd
import os
from datetime import datetime, timedelta
from glob import glob
import subprocess
import numpy as np
from threading import Thread
import time

def checkTimeRange(startDate,endDate,timeCoords):
    try:
        # Attempt to divide by zero
       if (timeCoords[0] <= startDate) & (endDate <= timeCoords[-1]):
           print("")
    except FileExceedTimeRange as e:
        # Handle the specific exception
        print(f"Caught an exception: {e}")

def update_varList(varList):
    varNum = len(varList)
    updated_varList = []
    for i in range(varNum):
        if varList[i]=='salinity':
            so_dict = {
                    "name": "salt",
                    "units": "Salinity in practical salinity units (PSU)",
                    "vars": ["salinity"],
                    "url":"http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_s3z/FMRC_ESPC-D-V02_s3z_best.ncd",
                    }
            updated_varList.append(so_dict)
        elif varList[i]=='water_temp':
            theta_dict = {
                    "name": "temp",
                    "units": "Temperature in degrees C",
                    "vars": ["water_temp"],
                    "url":"http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_t3z/FMRC_ESPC-D-V02_t3z_best.ncd",
                    }
            updated_varList.append(theta_dict)
        elif varList[i]=='surf_el':
            zos_dict = {
                    "name": "zeta",
                    "units": "Sea surface height in meters (m)",
                    "vars": ["surf_el"],
                    "url":"http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_ssh/FMRC_ESPC-D-V02_ssh_best.ncd",
                    }
            updated_varList.append(zos_dict)
        elif varList[i]=='water_u':
            uo_dict = {
                    "name": "u",
                    "units": "Eastward velocity component in meters per second (m.s-1)",
                    "vars": ["water_u"],
                    "url":"http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_u3z/FMRC_ESPC-D-V02_u3z_best.ncd",
                    }
            updated_varList.append(uo_dict)
        elif varList[i]=='water_v':
            vo_dict = {
                    "name": "v",
                    "units": "Northward velocity component in meters per second (m.s-1)",
                    "vars": ["water_v"],
                    "url":"http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_v3z/FMRC_ESPC-D-V02_v3z_best.ncd",
                    }
            updated_varList.append(vo_dict)
        else:
            print('Oops, check variables names in varList')
            break

    return updated_varList


def decode_time_units(time_var):
    """
    Convert a time variable from a NetCDF file to a pandas DatetimeIndex.

    Parameters:
    time_var: netCDF4.Variable
        A NetCDF variable with 'units' and optional 'calendar' attributes.

    Returns:
    pandas.DatetimeIndex
        A pandas DatetimeIndex corresponding to the time variable.
    """
    try:
        units = time_var.units
        calendar = getattr(time_var, 'calendar', 'standard')  # 'standard' as a common default
        time = cftime.num2date(
            time_var[:],
            units=units,
            calendar=calendar,
            only_use_cftime_datetimes=False,
            only_use_python_datetimes=True
        )
        return pd.DatetimeIndex(time)
    except AttributeError as e:
        raise ValueError(f"Missing expected attributes in time_var: {e}")
    except Exception as e:
        raise RuntimeError(f"Error decoding time units: {e}")

def download_var(varlist,domain,depths,savedir):
    """
    Main download function that gets called in download_hycom. 
    """
    # List of variables to remove from dataset
    varsToDrop = ['salinity_bottom','water_temp_bottom','water_u_bottom','water_v_bottom',
                  'tau','time_offset','time_run','time1_offset','sst','sss','ssu','ssv',
                  'sic','sih','siu','siv','surtx','surty','steric_ssh']
    
    # lon, lat and depth limits 
    lon_range = slice(domain[0],domain[1])
    lat_range = slice(domain[2],domain[3])
    depth_range = slice(depths[0],depths[1])
    
    ds = None
    MAX_RETRIES = 3
    RETRY_WAIT = 20
    engines = ["pydap","netcdf4"]
    i = 0
    while i < MAX_RETRIES:
        print('')
        print(f"Attempt {i+1} of {MAX_RETRIES}")
        print('')
        for engine in engines:
            try:
                ds = xr.open_dataset(varlist["url"], drop_variables=varsToDrop, decode_times = False, engine=engine)
                print(f"Successfully opened {varlist['vars'][0]} with engine: {engine}")
                # decode the time
                if 'time' in ds:
                    ds['time'] = decode_time_units(ds['time'])
                else:
                    pass
                    # Assign the correct variable
                if varlist["vars"][0]=='surf_el':
                    variable = ds.surf_el
                elif varlist["vars"][0]=='salinity':
                    variable = ds.salinity
                elif varlist["vars"][0]=='water_temp':
                    variable = ds.water_temp
                elif varlist["vars"][0]=='water_u':
                    variable = ds.water_u
                elif varlist["vars"][0]=='water_v':
                    variable = ds.water_v
                else:
                    print(varlist["vars"][0])
                    print('')
                    print(f'Invalid variable name: {varlist["vars"][0]}')
                    print('Valid variable names are: salinity, water_temp, surf_el, water_u, water_v')
                    print('')
                    break
                # Subset the variable
                if np.ndim(variable)==3:
                    var_subset = variable.sel({'lat' : lat_range, 'lon' : lon_range})
                else:
                    var_subset = variable.sel({'lat' : lat_range, 'lon' : lon_range, 'depth' : depth_range})
                # Compute daily averages (this can take a while)
                var_resampled = var_subset.resample(time='1D').mean()
                # Save and download the dataset
                sname = savedir + 'HYCOM_' + varlist["name"] + '.nc'
                var_resampled.to_netcdf(sname,'w')
                ds.close()
                i = MAX_RETRIES
                break
            except Exception as e:
                print(f"Failed to open {varlist['vars'][0]} with engine '{engine}': {e}")
        
        time.sleep(RETRY_WAIT)
        i+=1

    print('File saved: ', sname)

def download_hycom(variables,domain,depths,run_date,hdays,fdays,savedir):
    """
    Download HYCOM analysis using xarrray opendap.

    INPUTS:
    domain   : List of geographical coordinates to subset the data and download (e.g. [lon_min,lon_max,lat_min,lat_max]).
    depths   : List of minimum and maximum depths to download. Values must be positive (e.g. [0,5000]).
    variables: List of variables to download (default = ['salinity', 'water_temp', 'surf_el', 'water_u', 'water_v'])
    run_date : Todays datetime to ensure downloaded data corosponds (e.g. datetime.datetime(YYYY,MM,DD)).
    hdays    : Days to hindcast. Default is 5).
    fdays    : Days to forecast. Default is 5).

    OUTPUT:
    NetCDF file containing the most recent HYCOM forcast run.
    """
    varlist = update_varList(variables)
    def download_worker(var):
        download_var(var,domain,depths,savedir)
    
    threads = []
    for var in varlist:                
        p = Thread(target=download_worker, args=(var,))
        threads.append(p)
        p.start()
    
    for p in threads:
        p.join()
    
    start_date = pd.Timestamp(run_date) + timedelta(days=-hdays)
    end_date = pd.Timestamp(run_date) + timedelta(days=fdays)
    
    ds = xr.open_mfdataset(savedir + '*.nc')
    
    time_coords = ds.coords['time'].values
    
    checkTimeRange(start_date,end_date,time_coords)
    
    outfile = os.path.abspath(os.path.join(savedir, f"HYCOM_{run_date.strftime('%Y%m%d_00')}.nc"))
    
    if os.path.exists(outfile):
            os.remove(outfile)
    
    ds.to_netcdf(outfile,'w')
    
    subprocess.call(["chmod", "-R", "775", outfile])
    
    print('')
    print('created: ', outfile)
    print('')

