import xarray as xr
import cftime
import pandas as pd
import os
from datetime import datetime, timedelta
from glob import glob
import subprocess
import numpy as np
from threading import Thread

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
                    }
            updated_varList.append(so_dict)
        elif varList[i]=='water_temp':
            theta_dict = {
                    "name": "temp",
                    "units": "Temperature in degrees C",
                    "vars": ["water_temp"],
                    }
            updated_varList.append(theta_dict)
        elif varList[i]=='surf_el':
            zos_dict = {
                    "name": "zeta",
                    "units": "Sea surface height in meters (m)",
                    "vars": ["surf_el"],
                    }
            updated_varList.append(zos_dict)
        elif varList[i]=='water_u':
            uo_dict = {
                    "name": "u",
                    "units": "Eastward velocity component in meters per second (m.s-1)",
                    "vars": ["water_u"],
                    }
            updated_varList.append(uo_dict)
        elif varList[i]=='water_v':
            vo_dict = {
                    "name": "v",
                    "units": "Northward velocity component in meters per second (m.s-1)",
                    "vars": ["water_v"],
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

def download_hycom_parallel(ds,domain,depths,var,varname,outDir):
    fname = outDir + 'HYCOM_' + varname + '.nc'
    if os.path.exists(fname) is not True:      
        variable = ds[var]
        lon_range = slice(domain[0],domain[1])
        lat_range = slice(domain[2],domain[3])
        for var_name,data_array in variable.data_vars.items():
            ndims = data_array.ndim
        if ndims == 4:
            for n in range(len(depths)):
                if depths[n] < 10:
                    str_depth = '000' + str(depths[n])
                elif (depths[n] >= 10) & (depths[n] < 100):
                    str_depth = '00' + str(depths[n])
                elif (depths[n] >= 100) & (depths[n] < 1000):
                    str_depth = '0' + str(depths[n])
                else:
                    str_depth = str(depths[n])
                sname = outDir + 'HYCOM_' + varname + '_' + str_depth + '.nc'
                if os.path.exists(sname):
                    continue
                else:
                    var_subset = variable.sel({'lat' : lat_range, 'lon' : lon_range, 'depth' : depths[n]})
                    var_resampled = var_subset.resample(time='1D').mean()
                    var_resampled.to_netcdf(sname,'w')
                    var_resampled.close()
            
            filepath = outDir + 'HYCOM_' + varname + '_' + '*.nc'
            file_paths = sorted(glob(filepath))
            combined_ds = xr.open_mfdataset(file_paths, concat_dim='depth', combine='nested')
            ds_reordered = combined_ds.transpose('time','depth','lat','lon')
            ds_reordered.to_netcdf(fname,'w')
            ds_reordered.close()
            
            delFiles = glob(outDir + 'HYCOM_' + varname + '_' + '*.nc')
            
            for file in range(len(delFiles)):
                subprocess.call(["rm", "-rf", delFiles[file]])
        
        elif ndims == 3:
            var_subset = variable.sel({'lat' : lat_range, 'lon' : lon_range})
            var_resampled = var_subset.resample(time='1D').mean()
            var_resampled.to_netcdf(fname,'w')
            var_resampled.close()
        else:
            print('')
            print('Error: number of dimensions are all wrong!')
            print('')
    else:
        print('')
        print('The file already exists: ', varname + '.nc')
        print('')

def download_hycom_series(ds,domain,depths,varList,outDir):
    for vn in range(len(varList)):
        fname = outDir + 'HYCOM_' + varList[vn] + '.nc'

        variable = ds[varList[vn]]

        lon_range = slice(domain[0],domain[1])
        
        lat_range = slice(domain[2],domain[3])

        if variable.ndim == 4:
            for n in range(len(depths)):
                if depths[n] < 10:
                    str_depth = '000' + str(depths[n])
                elif (depths[n] >= 10) & (depths[n] < 100):
                    str_depth = '00' + str(depths[n])
                elif (depths[n] >= 100) & (depths[n] < 1000):
                    str_depth = '0' + str(depths[n])
                else:
                    str_depth = str(depths[n])
                
                sname = outDir + 'HYCOM_' + varList[vn]+ '_' + str_depth + '.nc'
                var_subset = variable.sel({'lat' : lat_range, 'lon' : lon_range, 'depth' : depths[n]})
                var_resampled = var_subset.resample(time='1D').mean()
                var_resampled.to_netcdf(sname,'w')
                var_resampled.close()

            filepath = sorted(glob(outDir + 'HYCOM_' + varList[vn] + '_' + '*.nc'))
            combined_ds = xr.open_mfdataset(filepath, concat_dim='depth', combine='nested')
            ds_reordered = combined_ds.transpose('time','depth','lat','lon')
            ds_reordered.to_netcdf(fname,'w')
            ds_reordered.close()
            delFiles = glob(outDir + 'HYCOM_' + varList[vn] + '_' + '*.nc')
            for file in range(len(delFiles)):
                subprocess.call(["rm", "-rf", delFiles[file]])

        elif variable.ndim == 3:
            var_subset = variable.sel({'lat' : lat_range, 'lon' : lon_range})
            var_resampled = var_subset.resample(time='1D').mean()
            var_resampled.to_netcdf(fname,'w')
            var_resampled.close()
        
        else:
            print('')
            print('Error: number of dimensions are all wrong!')
            print('')
   
def download_hycom(outDir,domain=None,depths=None,varList=None,run_date=None,hdays=None,fdays=None,parallel=True):
    """
    Download HYCOM analysis using xarrray opendap.

    INPUTS:
    outDir   : Directory to save files in (e.g. "/path/to/save/files/in/").
    domain   : List of geographical coordinates to subset the data and download (e.g. []).
    depths   : List of depths to download (e.g. [lomin,lonmax,latmin,latmax]).
    varList  : List of variables to download (e.g. [0,2,4,...,3000,4000,5000])
    run_date : Todays datetime to ensure downloaded data corosponts (e.g. ).
    hdays    : Days to hindcast (e.g. hdays = 5).
    fdays    : Days to forecast (e.g. fdays = 5).
    parallel : Boolean to download in parallel or series (True = download the files in parallel | False = download the files in series).

    OUTPUT:
    NetCDF file containing the most recent HYCOM forcast run.
    """
    print('')
    if domain is None:
        domain = [23,34,-37,-31]
        print('No domain specified. Using dafault domain:')
        print(domain)
        print('')

    if depths is None:
        depths = [0,2,4,6,8,10,12,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,200,250,300,350,400,500,600,700,800,900,1000,1250,1500,2000,2500,3000,4000,5000]
        print('No depths specified. Downloading all depths: ')
        print(depths)
        print('')

    if varList is None:
        varList = ['surf_el','salinity','water_temp','water_v','water_u']
        print('No variables specified. Downloading all variables: ')
        print(varList)
        print('')

    if run_date is None:
        run_date = pd.Timestamp.now()
        print('No run_date specified. Using current run_date: ')
        print(run_date)
        print('')

    if hdays is None:
        hdays = 5
        print('No hdays specified. Using hdays: ')
        print(hdays)
        print('')


    if fdays is None:
        fdays = 5
        print('No fdays specified. Using fdays: ')
        print(fdays)
        print('')

    start_date = pd.Timestamp(run_date) + timedelta(days=-hdays)

    end_date = pd.Timestamp(run_date) + timedelta(days=fdays)

    url = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/FMRC/GLBy0.08_930_FMRC_best.ncd'

    varsToDrop = ['salinity_bottom','water_temp_bottom','water_u_bottom','water_v_bottom','tau','time_offset','time_run']
    
    ds = xr.open_dataset(url, drop_variables=varsToDrop, decode_times = False)

    if 'time' in ds:
        ds['time'] = decode_time_units(ds['time'])

    if parallel:
        varList = update_varList(varList)

        def download_worker(var):
            download_hycom_parallel(ds,domain,depths,var['vars'],var['name'],outDir)

        threads = []
        for var in varList:                
            p = Thread(target=download_worker, args=(var,))
            threads.append(p)
            p.start()

        for p in threads:
            p.join()

    else:
        download_hycom_series(ds,domain,depths,varList,outDir)

    ds.close()

    ds = xr.open_mfdataset(outDir + '*.nc')

    time_coords = ds.coords['time'].values

    checkTimeRange(start_date,end_date,time_coords)

    outfile = os.path.abspath(os.path.join(outDir, f"HYCOM_{run_date.strftime('%Y%m%d_%H')}.nc"))

    if os.path.exists(outfile):
            os.remove(outfile)

    ds.to_netcdf(outfile,'w')

    subprocess.call(["chmod", "-R", "775", outfile])

    print('')
    print('created: ' + outDir + 'HYCOM.nc')
    print('')

