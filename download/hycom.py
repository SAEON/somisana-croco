import xarray as xr
import cftime
import pandas as pd
import os
from datetime import datetime, timedelta
import subprocess
import numpy as np

def checkTimeRange(startDate,endDate,timeCoords):
    try:
        # Attempt to divide by zero
       if (timeCoords[0] <= startDate) & (endDate <= timeCoords[-1]):
           print("")
    except FileExceedTimeRange as e:
        # Handle the specific exception
        print(f"Caught an exception: {e}")

def update_varList(varList):
    # Define the mapping for variables to metadata
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

    # Build the dictionary for the output
    updated_varDict = {}

    for var in varList:
        if var in var_metadata:
            updated_varDict[var] = var_metadata[var]
        else:
            print(f"Oops, '{var}' is an invalid variable name. Please check variable names in varList.")

    return updated_varDict

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

def download_vars(variables,domain,depths,savedir):    
    for var in variables:
        varlist = update_varList([var])
        print(f'Connecting to {varlist[var]["url"]} to subset and download {varlist[var]["vars"][0]}.') 
        
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
        download_engine = "netcdf4"
        i = 0
        while i < MAX_RETRIES:
            print('')
            print(f"Attempt {i+1} of {MAX_RETRIES}")
            print('')
        
            ds = xr.open_dataset(varlist[var]["url"], drop_variables=varsToDrop, decode_times = False, engine="netcdf4")
            if 'time' in ds:
                ds['time'] = decode_time_units(ds['time'])
            else:
                pass
        
            # Assign the correct variable
            if varlist[var]["vars"][0]=='surf_el':
                variable = ds.surf_el
            elif varlist[var]["vars"][0]=='salinity':
                variable = ds.salinity    
            elif varlist[var]["vars"][0]=='water_temp':
                variable = ds.water_temp    
            elif varlist[var]["vars"][0]=='water_u':
                variable = ds.water_u
            elif varlist[var]["vars"][0]=='water_v':
                variable = ds.water_v
            else:
                print(f'Invalid variable name: {varlist[var]["vars"]}')
                print(f"Valid variable names are: salinity, water_temp, surface_el, water_u and water_v")
                print('')

            # Subset the variable
            if np.ndim(variable)==3:
                var_subset = variable.sel({'lat' : lat_range, 'lon' : lon_range})
            else:
                var_subset = variable.sel({'lat' : lat_range, 'lon' : lon_range, 'depth' : depth_range})
            
            # Compute daily averages (this can take a while)
            print(f'Subsetting {varlist[var]["vars"][0]}.')
            print('')
            var_resampled = var_subset.resample(time='1D').mean()
            
            # Save and download the dataset
            print(f'Downloading {varlist[var]["vars"][0]}.')
            print('')
            sname = savedir + 'hycom_' + varlist[var]["vars"][0] + '.nc'
            var_resampled.to_netcdf(sname,'w')
            
            if ds is not None:
                ds.close()    
                print(f'File written to {sname}')
                print('')
                i=MAX_RETRIES
            else:
                i+=1

def download_hycom(variables,domain,depths,run_date,hdays,fdays,savedir):
    """
    Download HYCOM analysis using xarrray opendap.

    INPUTS:
    domain   : List of geographical coordinates to subset the data and download (e.g. [lon_min,lon_max,lat_min,lat_max]).
    depths   : List of minimum and maximum depths to download. Values must be positive (e.g. [0,5000]).
    variables: List of variables to download (e.g. ['salinity', 'water_temp', 'surface_el', 'water_u', 'water_v'])
    run_date : Todays datetime to ensure downloaded data corosponds (e.g. datetime.datetime(YYYY,MM,DD)).
    hdays    : Days to hindcast (e.g. hdays = 5).
    fdays    : Days to forecast (e.g. fdays = 5).

    OUTPUT:
    NetCDF file containing the most recent HYCOM forcast run.
    """
    
    download_vars(variables,domain,depths,savedir)

    hdays=hdays+1
    fdays=fdays+1
    
    start_date = pd.Timestamp(run_date) + timedelta(days=-hdays)
    end_date = pd.Timestamp(run_date) + timedelta(days=fdays)
    
    ds = xr.open_mfdataset(savedir + 'hycom_*.nc')
    
    time_coords = ds.coords['time'].values
    
    checkTimeRange(start_date,end_date,time_coords)
    
    outfile = os.path.abspath(os.path.join(savedir, f"HYCOM_{run_date.strftime('%Y%m%d_%H')}.nc"))
    
    if os.path.exists(outfile):
            os.remove(outfile)
    
    ds.to_netcdf(outfile,'w')
    
    subprocess.call(["chmod", "-R", "775", outfile])
    
    print('created: ', outfile)
    print('')
