import xarray as xr
from datetime import datetime, timedelta
import calendar
import numpy as np
import pandas as pd
import os, sys, glob
# import fnmatch
from scipy.interpolate import griddata
import crocotools_py.postprocess as post
# functions from croco_pytools submodule
# (I'm sure there's a better way of importing these functions, but this works)
sys.path.append(os.path.dirname(__file__) + "/croco_pytools/prepro/Modules/")
sys.path.append(os.path.dirname(__file__) + "/croco_pytools/prepro/Readers/")
import Cgrid_transformation_tools as grd_tools
import interp_tools
import sigmagrid_tools as sig_tools
import croco_class as Croco
import ibc_class as Inp
import pylab as plt
import netCDF4 as netcdf

def fill_blk(croco_grd,croco_blk_file_in,croco_blk_file_out):
    '''
    This is a little hack function to help us getting around the specific problem
    of having random blocks of missing data in the WASA3 zarr files 
    So we're just using this function to interpolate over these blocks.
    
    I'm specifically not doing this automatically in make_WASA3_from_blk()
    as the user should really know exactly where the missing data are and 
    first see if this hack is appropriate
    
    Parameters
    ----------
    croco_grd      : /your_dir/croco_grd.nc
    croco_blk_file_in : /your_dir/croco_blk_file_needing_filling
    croco_blk_file_out : /your_dir/filled_croco_blk_file

    '''
    # get the croco grid variables
    ds_croco_grd=xr.open_dataset(croco_grd)
    # handling uwnd and vwnd on their own grids
    lon_u=ds_croco_grd.lon_u.values
    lat_u=ds_croco_grd.lat_u.values
    lon_v=ds_croco_grd.lon_v.values
    lat_v=ds_croco_grd.lat_v.values
    ds_croco_grd.close()
    # for use in griddata later:
    source_points_u = np.column_stack((lon_u.flatten(), lat_u.flatten()))
    source_points_v = np.column_stack((lon_v.flatten(), lat_v.flatten()))
    
    ds_blk = xr.open_dataset(croco_blk_file_in, decode_times=False)
    u_blk_filled = ds_blk.uwnd.values
    v_blk_filled = ds_blk.vwnd.values
    for t in range(len(ds_blk.bulk_time)):
        
        # Extract wasa values at this specific time-step
        # flattening for use in griddata
        u_blk_now = u_blk_filled[t,:,:].flatten()
        v_blk_now = v_blk_filled[t,:,:].flatten()
        
        # check for missing values and fill them where needed
        # handling u and v separately as they're on their own grids with different indices
        if np.any(np.isnan(u_blk_now)):
            # get the indices for which we have valid data
            idx = np.logical_not(np.isnan(u_blk_now))
            # Perform interpolation using valid data
            u_blk_filled[t,:,:] = griddata(source_points_u[idx,:], u_blk_now[idx], (lon_u, lat_u), method='nearest')
        if np.any(np.isnan(v_blk_now)):
            # get the indices for which we have valid data
            idx = np.logical_not(np.isnan(v_blk_now))
            # Perform interpolation using valid data
            v_blk_filled[t,:,:] = griddata(source_points_v[idx,:], v_blk_now[idx], (lon_v, lat_v), method='nearest')
        
    # update wspd (on the rho grid)
    spd_blk_filled = np.hypot(post.u2rho(u_blk_filled), post.v2rho(v_blk_filled))
    
    # Update uwnd in ds_blk
    ds_blk['uwnd'] = (('bulk_time', 'eta_u', 'xi_u'), u_blk_filled)
    ds_blk['vwnd'] = (('bulk_time', 'eta_v', 'xi_v'), v_blk_filled)
    ds_blk['wspd'] = (('bulk_time', 'eta_rho', 'xi_rho'), spd_blk_filled)
    
    # write float32 instead of float64 in an attempt to save space
    encoding = {
        "tair": {"dtype": "float32"},
        "rhum": {"dtype": "float32"},
        "prate": {"dtype": "float32"},
        "wspd": {"dtype": "float32"},
        "radlw_in": {"dtype": "float32"},
        "radsw": {"dtype": "float32"},
        "uwnd": {"dtype": "float32"},
        "vwnd": {"dtype": "float32"},
        "bulk_time": {"dtype": "float32"},
    }
    
    # write the new blk file
    ds_blk.to_netcdf(croco_blk_file_out,
                     mode='w', 
                     encoding=encoding)
    
    ds_blk.close()


def subset_WASA3(ds_wasa,ds_wasa_grd,extents,dl=0.1):
    '''
    Parameters
    ----------
    ds_wasa : xarray dataset extracted from wasa3
    ds_wasa_grd : xarray dataset extracted from wasa3 grid file
    extents : [lon0,lon1,lat0,lat1]
    dl      : buffer to add around the extents (degrees)

    Returns
    -------
    spatially subsetted dataset based on the input extents

    '''
    
    def find_nearest_wasa_point(lon_wasa, lat_wasa, lon_in, lat_in):
        # subfunction to find the nearest indices of the wasa grid to a specific lon, lat coordinate
        # based on postprocess.find_nearest_point() in this repo
        
        # Calculate the distance between (lon_in, lat_in) and all grid points
        distance = ((lon_wasa - lon_in) ** 2 +
                    (lat_wasa - lat_in) ** 2) ** 0.5

        # Find the indices of the minimum distance
        # unravel_index method Converts a flat index or array of flat indices into a tuple of coordinate 
        # arrays: https://numpy.org/doc/stable/reference/generated/numpy.unravel_index.html
        min_index = np.unravel_index(distance.argmin(), distance.shape)

        j, i = min_index

        return j, i
    
    lon_wasa = ds_wasa_grd.XLONG_M.squeeze()
    lat_wasa = ds_wasa_grd.XLAT_M.squeeze()
    
    # get the nearest wasa grid indices corresponding to the corners of the croco grid   
    # bl = bottom left, br = bottom right, tr = top right, tl = top_left
    j_bl, i_bl = find_nearest_wasa_point(lon_wasa, lat_wasa, extents[0] - dl, extents[2] - dl)
    j_br, i_br = find_nearest_wasa_point(lon_wasa, lat_wasa, extents[1] + dl, extents[2] - dl)
    j_tr, i_tr = find_nearest_wasa_point(lon_wasa, lat_wasa, extents[1] + dl, extents[3] + dl)
    j_tl, i_tl = find_nearest_wasa_point(lon_wasa, lat_wasa, extents[0] - dl, extents[3] + dl)
    # use the corners to get the extreme indices
    j_min=min(j_bl,j_br)
    j_max=max(j_tl,j_tr)
    i_min=min(i_bl,i_tl)
    i_max=max(i_br,i_tr)
    #
    # finally, subset wasa datasets based on these indices
    ds_wasa=ds_wasa.isel(south_north=slice(j_min,j_max+1),
                         west_east=slice(i_min,i_max+1))
    ds_wasa_grd=ds_wasa_grd.isel(south_north=slice(j_min,j_max+1),
                         west_east=slice(i_min,i_max+1))
    
    return ds_wasa,ds_wasa_grd

def make_WASA3_from_blk(wasa_grid, 
                 wasa_zarr_dir, 
                 croco_grd,
                 croco_blk_dir,
                 out_wasa_dir,
                 ref_date,
                 start_date,
                 end_date,
                 croco_atmos='ERA5',
                 interp_method='linear'
                 ):
    '''
    Interpolate WASA3 u10, v10 data onto existing blk files
    Recommended to use ERA5 for the existing blk files, 
    as WASA3 is a downscaling of ERA5
    So we're just increasing the resolution of u10, v10
    
    Parameters
    ----------
    wasa_grid     : /your_dir/geo_em.d03.nc
    wasa_zarr_dir : /your_dir/uv_ds_10m
    croco_grd     : /your_dir/croco_grd.nc
    croco_blk_dir : the directory containing your existing blk files
    out_wasa_dir  : the directory to save blk files with the WASA data
    ref_date      : datetime object corresponding to the reference date
    croco_atmos   : the string defining the existing atmospheric forcing
                    default is 'ERA5'.
    interp_method : The default is 'linear'. Can use 'nearest' to speed things up?

    Returns
    -------
    Writes blk netcdf files in out_wasa_dir

    '''
        
    # get the relevant wasa data
    ds_wasa=xr.open_zarr(wasa_zarr_dir)
    ds_wasa_grd=xr.open_dataset(wasa_grid)
    
    # get the croco grid variables
    ds_croco_grd=xr.open_dataset(croco_grd)
    lon_rho=ds_croco_grd.lon_rho
    lat_rho=ds_croco_grd.lat_rho
    angle=ds_croco_grd.angle
    ds_croco_grd.close()
    
    # Subset WASA spatially to speed up interpolation later
    # this helps a lot so worth the mess!
    extents=[min(lon_rho[0,0],lon_rho[-1,0]), # western extent
             max(lon_rho[0,-1],lon_rho[-1,-1]), # eastern extent
             min(lat_rho[0,0],lat_rho[0,-1]), # southern extent
             max(lat_rho[-1,0],lat_rho[-1,-1])] # northern extent
    ds_wasa,ds_wasa_grd=subset_WASA3(ds_wasa,ds_wasa_grd,extents)
        
    # extract the wasa grid variables we need from this spatial subset
    lon_wasa = ds_wasa_grd.XLONG_M.squeeze()
    lat_wasa = ds_wasa_grd.XLAT_M.squeeze() 
    cosa = ds_wasa_grd.COSALPHA.squeeze()
    sina = ds_wasa_grd.SINALPHA.squeeze()
    # flatten lon,lat arays for use in interpolation later
    source_points = np.column_stack((np.array(lon_wasa).flatten(), np.array(lat_wasa).flatten()))

    date_now=start_date
    while date_now <= end_date:
        croco_blk_file = 'croco_blk_'+croco_atmos+str(date_now.strftime('_Y%YM%m'))+'.nc'
   
        # get the corresponding WASA blk filename
        croco_blk_file_splt=croco_blk_file.split(croco_atmos)
        croco_blk_file_WASA=croco_blk_file_splt[0]+'WASA3'+croco_blk_file_splt[1]
        
        if os.path.isfile(out_wasa_dir+'/'+croco_blk_file_WASA):
            print(croco_blk_file_WASA + ' already exists - skipping')
        else:
            print('\nworking on '+croco_blk_file)
            
            # get a dataset for the croco blk file
            # NB to have decode_times=False
            ds_blk = xr.open_dataset(croco_blk_dir+'/'+croco_blk_file, decode_times=False)
            
            # get an array of datetimes for this file
            blk_time = np.array([ref_date + timedelta(days=day) for day in ds_blk.bulk_time.values])
            
            # subset the WASA data based on the blk times
            # let's just get the closest WASA data to our blk time-steps
            # WASA data is actually half hourly! 
            # So with this method we're only taking
            # every second WASA record if our blk file is hourly!
            ds_wasa_now = ds_wasa.sel(Time=blk_time, method="nearest")
            # 'Time' in ds_wasa_now is now the same length as 'bulk_time' in ds_blk
            
            # rotate the wasa vectors to be east,north components
            u = ds_wasa_now.U10
            v = ds_wasa_now.V10
            u_wasa = u * cosa - v * sina
            v_wasa = v * cosa + u * sina
            
            # check if there are any missing data and print a warning
            # This is needed due to some random blocks of missing values we've found
            # I'm intentionally not automatically interpolating over these data
            # as we should first check to make sure we are happy to do so
            # ultimately, we need to get a un-corrupted copy of the data!
            if np.any(np.isnan(u_wasa.values.flatten())):
                print('wasa contains missing data for '+croco_blk_file_WASA+'. Check this out!')
            elif np.any(np.isnan(v_wasa.values.flatten())):
                print('wasa contains missing data for '+croco_blk_file_WASA+'. Check this out!')
            
            # interpolate the wasa u,v wind vectors onto the croco grid   
            # I tried to find an in-built xarray interpolation method to do the job
            # efficiently, but it only works for regular grids
            # So we've got to extract the data and use scipy's griddata function
            # Note I'm puposefully interpolating onto the croco rho grid (not the croco u,v grids)
            # because we have to rotate the u,v vectors (on the rho grid)
            # before writing the data to file, so would have to use u2rho() and v2rho() 
            # which would introduce additional unnecessary interpolation
            u_wasa_interp = np.zeros_like(ds_blk.wspd) # using 'wspd' to get the rho grid
            v_wasa_interp = np.zeros_like(ds_blk.wspd)
            for t in range(len(blk_time)):
                if t % 50 == 0:
                    percentage_complete = (t/len(blk_time)) * 100
                    print(f"{percentage_complete:.0f}% complete")
                
                # Extract wasa values at this specific time-step
                # flattening for use in griddata
                u_wasa_now = u_wasa.isel(Time=t).values.flatten()
                v_wasa_now = v_wasa.isel(Time=t).values.flatten()
                
                # Perform interpolation 
                # 'nearest' method speeds it up a bit
                # but I guess 'linear' is better
                u_wasa_now_interp = griddata(source_points, u_wasa_now, (lon_rho, lat_rho), method=interp_method) 
                v_wasa_now_interp = griddata(source_points, v_wasa_now, (lon_rho, lat_rho), method=interp_method)             
                
                # hmmm, with 'linear' interpolation we're left with some nan's in the case where our domain 
                # extends a little past the WASA data, so default to nearest wasa data here
                # this will slow things down a bit, but not sure what else to do, other than regenerate the croco grid!
                if np.any(np.isnan(u_wasa_now_interp)):
                    u_wasa_now_interp_nearest = griddata(source_points, u_wasa_now, (lon_rho, lat_rho), method='nearest') 
                    v_wasa_now_interp_nearest = griddata(source_points, v_wasa_now, (lon_rho, lat_rho), method='nearest')
                    # now replace missing data with nearest available data
                    u_wasa_now_interp[np.isnan(u_wasa_now_interp)] = u_wasa_now_interp_nearest[np.isnan(u_wasa_now_interp)]
                    v_wasa_now_interp[np.isnan(v_wasa_now_interp)] = v_wasa_now_interp_nearest[np.isnan(v_wasa_now_interp)]
                
                u_wasa_interp[t,:,:] = u_wasa_now_interp
                v_wasa_interp[t,:,:] = v_wasa_now_interp
                
            # calculate the wspd from the wasa components (on the rho grid)
            spd_wasa_interp = np.hypot(u_wasa_interp, v_wasa_interp)
            
            # uwnd and vwnd are grid aligned, so we need to rotate vectors
            # (based on official croco_tools interp_ERA5.m)
            cos_a = np.cos(angle.values)
            sin_a = np.sin(angle.values)
            # have to use the data on the rho grid for the rotation calc 
            u_out = u_wasa_interp*cos_a + v_wasa_interp*sin_a
            v_out = v_wasa_interp*cos_a - u_wasa_interp*sin_a
            # but get onto the u,v grids to write to file
            u_out = post.rho2u(u_out)
            v_out = post.rho2v(v_out)
            
            # Update uwnd in ds_blk
            ds_blk['uwnd'] = (('bulk_time', 'eta_u', 'xi_u'), u_out)
            ds_blk['vwnd'] = (('bulk_time', 'eta_v', 'xi_v'), v_out)
            ds_blk['wspd'] = (('bulk_time', 'eta_rho', 'xi_rho'), spd_wasa_interp)
            
            # I'm pretty sure we don't need sustr and svstr in blk files as this is computed online
            # So I'm removing these variables to be sure
            # leaving them hanging around seems a little messy since we're replacing uwnd and vwnd
            # while we're at it, radlw is also not used (radlw_in is), so remove that as well
            ds_blk = ds_blk.drop_vars(["sustr", "svstr", "radlw"])
            
            # write float32 instead of float64 in an attempt to save space
            encoding = {
                "tair": {"dtype": "float32"},
                "rhum": {"dtype": "float32"},
                "prate": {"dtype": "float32"},
                "wspd": {"dtype": "float32"},
                "radlw_in": {"dtype": "float32"},
                "radsw": {"dtype": "float32"},
                "uwnd": {"dtype": "float32"},
                "vwnd": {"dtype": "float32"},
                "bulk_time": {"dtype": "float32"},
            }
            
            # write the new blk file
            ds_blk.to_netcdf(out_wasa_dir+'/'+croco_blk_file_WASA,
                             encoding=encoding)
            
            ds_blk.close()
            
        date_now=date_now+timedelta(days=32) # 32 days ensures we get to the next month
        date_now=datetime(date_now.year, date_now.month, 1) # set the first day of the month 
        
    ds_wasa.close()
    ds_wasa_grd.close()

def WASA3_2_nc(wasa_grid, 
                 wasa_zarr_dir, 
                 out_dir,
                 start_date,
                 end_date,
                 extents=None
                 ):
    '''
    Create monthly nc files from WASA3 u10, v10 data
    To be used as OpenDrift forcing
    
    Parameters
    ----------
    wasa_grid     : /your_dir/geo_em.d03.nc
    wasa_zarr_dir : /your_dir/uv_ds_10m
    out_wasa_dir  : the directory to save blk files with the WASA data
    start_date    : start datetime for processing (to the nearest month)
    end_date      : end datetime for processing (to the nearest month)
    extents       : option to only get a spatial subset [lon0,lon1,lat0,lat1]
    
    Returns
    -------
    Writes monthly netcdf files in out_dir

    '''
        
    # get the relevant wasa data
    ds_wasa=xr.open_zarr(wasa_zarr_dir)
    ds_wasa_grd=xr.open_dataset(wasa_grid)
    
    if extents is not None:
        ds_wasa,ds_wasa_grd=subset_WASA3(ds_wasa,ds_wasa_grd,extents)
    
    # extract the wasa grid variables we need from this spatial subset
    lon_wasa = ds_wasa_grd.XLONG_M.squeeze()
    lat_wasa = ds_wasa_grd.XLAT_M.squeeze() 
    cosa = ds_wasa_grd.COSALPHA.squeeze()
    sina = ds_wasa_grd.SINALPHA.squeeze()
    
    # prefer to keep rename Time to time
    ds_wasa = ds_wasa.rename({'Time': 'time'})
    
    # subset time so that we only keep every second time-step
    # The raw data are half hourly so this makes it hourly
    # plenty for our purposes and we half the disk space
    ds_wasa = ds_wasa.isel(time=slice(None, None, 2))
    
    # loop through months
    date_now=start_date
    while date_now <= end_date:
        file_out = 'WASA3'+date_now.strftime('_Y%YM%m')+'.nc'
   
        if os.path.isfile(out_dir+'/'+file_out):
            print(file_out + ' already exists - skipping')
        else:
            print('\nworking on '+file_out)
            
            # start and end days of this month
            start_extract=datetime(date_now.year,date_now.month,1,0,0,0)
            day_end = calendar.monthrange(date_now.year, date_now.month)[1]
            end_extract=datetime(date_now.year,date_now.month,day_end,23,59,59)
            
            # subset the WASA data for this month
            ds_wasa_now = ds_wasa.sel(time=slice(start_extract,end_extract))
            time_now = ds_wasa_now.time
            
            # rotate the wasa vectors to be east,north components
            u = ds_wasa_now.U10
            v = ds_wasa_now.V10
            u_wasa = u * cosa - v * sina
            v_wasa = v * cosa + u * sina
            
            # check if there are any missing data
            # This is needed due to some random blocks of missing values we've found
            # I'm automatically interpolating over these data so we have workable files
            # THIS IS A QUICK FIX AND NEEDS TO BE SORTED OUT WITH CSAG!
            # Then we can remove this part of the code
            if np.any(u_wasa.isnull().data.flatten()): # this check works for both u and v due to the rotation calc above - any missing data in one translates to missing in the other
                print('filling missing values in '+file_out)
                # Rechunk the data to avoid issues with Dask and interpolation below
                u_wasa = u_wasa.chunk({'south_north': -1, 'west_east': -1})
                u_wasa = u_wasa.interpolate_na(dim='south_north', method='linear', fill_value='extrapolate')
            
                # Rechunk the data to avoid issues with Dask and interpolation below
                v_wasa = v_wasa.chunk({'south_north': -1, 'west_east': -1})
                v_wasa = v_wasa.interpolate_na(dim='south_north', method='linear', fill_value='extrapolate')
            
            # add a couple of attributes
            u_wasa.attrs["standard_name"] = 'x_wind'
            u_wasa.attrs["units"] = 'm s**-1'
            v_wasa.attrs["standard_name"] = 'y_wind'
            v_wasa.attrs["units"] = 'm s**-1'
            
            # Create the dataset
            ds_out = xr.Dataset(
                {
                    'time': time_now,
                    'longitude': lon_wasa,
                    'latitude': lat_wasa,
                    'x_wind': u_wasa,
                    'y_wind': v_wasa,
                }
            )
            ds_out.to_netcdf(out_dir+'/'+file_out)
            
        date_now=date_now+timedelta(days=32) # 32 days ensures we get to the next month
        date_now=datetime(date_now.year, date_now.month, 1) # set the first day of the month 
        
    ds_wasa.close()
    ds_wasa_grd.close()  

    
def make_SAWS_from_blk(saws_dir, 
                 croco_grd,
                 croco_blk_in,
                 croco_blk_out,
                 ref_date,
                 croco_atmos='GFS',
                 interp_method='linear'
                 ):
    '''
    Interpolate operational SAWS UM data onto existing blk files
    
    Parameters
    ----------
    ...
    ref_date      : datetime object corresponding to the reference date
    interp_method : The default is 'linear'. Can use 'nearest' to speed things up?

    Returns
    -------
    Writes blk netcdf file with SAWS data

    '''
    
    # get the croco grid variables...
    ds_croco_grd=xr.open_dataset(croco_grd)
    lon_rho=ds_croco_grd.lon_rho
    lat_rho=ds_croco_grd.lat_rho
    angle=ds_croco_grd.angle
    ds_croco_grd.close()
    # ... and extract the spatial limits to be used in subsetting the SAWS data
    # (useful to speed up the code)
    dl=0.5
    lon_min=min(lon_rho[0,0],lon_rho[-1,0])-dl # western extent
    lon_max=max(lon_rho[0,-1],lon_rho[-1,-1])+dl # eastern extent
    lat_min=min(lat_rho[0,0],lat_rho[0,-1])-dl # southern extent
    lat_max=max(lat_rho[-1,0],lat_rho[-1,-1])+dl # northern extent
    
    # Get a list of all .nc files in the saws directory...
    files = sorted(glob.glob(f"{saws_dir}/*.nc"))
    # ... and reverse the order to give priority to the latest files
    # (could be a weak point in the code if the file names change?)
    files.reverse()
    
    # Open datasets and combine them, giving priority to the latest files
    # this is needed since the times in the files overlap
    print('extracting SAWS UM data')
    ds_saws = None
    for file in files:
        ds = xr.open_dataset(file)
        # do the spatial subset (rather here than later to save time in the concatenation step, which is a bit slow)
        ds = ds.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))
        if ds_saws is None:
            ds_saws = ds
        else:
            # extract the non-overlapping data
            ds = ds.sel(time=slice(ds.time[0],ds_saws.time[0]-1))
            # and combine with the full dataset
            ds_saws = xr.concat([ds, ds_saws], dim='time')
        ds.close()   
    
    print('extracting CROCO blk data')
    # get a dataset for the croco blk file
    # NB to have decode_times=False (since CROCO files don't come with the ref date in the attributes)
    ds_blk = xr.open_dataset(croco_blk_in, decode_times=False)
    # get an array of datetimes for this file
    blk_time = np.array([ref_date + timedelta(days=day) for day in ds_blk.bulk_time.values])
    
    # now get the blk file onto the SAWS time steps
    # we are getting the blkfile times onto the SAWS times rather than the other way around
    # because the saws data only contains 3 days of forecasts, so it will be a subset of the normal blk file which would have 5 days of forecast data
    # so the output blk file will have fewer time-steps than the original
    #
    # start by cutting off saws times which are outside the blk file (we may have read in more hindcast data than we need)
    ds_saws = ds_saws.sel(time=slice(blk_time[0],blk_time[-1]))
    # and then subset the bulk times to the overlapping SAWS data
    # (to do this we need to get the saws times in days since ref_date, as per the ds_blk.bulk_time)
    # (we intentionally don't change ds_blk.bulk_time into real datetimes because the output blk must be the same format as the input blk)
    saws_days_since_ref_date = (pd.to_datetime(ds_saws.time.values) - ref_date).total_seconds()/86400
    ds_blk = ds_blk.sel(bulk_time=saws_days_since_ref_date, method="nearest")
    # now ds_blk and ds_wasa have the exact same number of timesteps
    
    # extract the saws grid variables for doing the spatial interpolation
    lon_saws = ds_saws.lon.squeeze()
    lat_saws = ds_saws.lat.squeeze() 
    lon_saws, lat_saws = np.meshgrid(lon_saws,lat_saws)
    # flatten lon,lat arays for use in interpolation later
    source_points = np.column_stack((np.array(lon_saws).flatten(), np.array(lat_saws).flatten()))
    
    # interpolate the saws data onto the croco grid   
    # I tried to find an in-built xarray interpolation method to do the job
    # efficiently, but it only works for regular grids
    # So we're extracting the data and using scipy's griddata function
    # Note I'm puposefully interpolating saws wind components onto the croco rho grid (not the croco u,v grids)
    # because we have to rotate the u,v vectors (on the rho grid)
    # before writing the data to file, so would have to use u2rho() and v2rho() 
    # which would introduce additional unnecessary interpolation
    print('interpolating SAWS data onto CROCO grid...')
    u_saws_interp = np.zeros_like(ds_blk.wspd) # using 'wspd' to get the rho grid
    v_saws_interp = np.zeros_like(ds_blk.wspd)
    t2_saws_interp = np.zeros_like(ds_blk.wspd)
    rhum_saws_interp = np.zeros_like(ds_blk.wspd)
    for t in range(len(ds_blk.bulk_time)):
        if t % 50 == 0:
            percentage_complete = (t/len(ds_blk.bulk_time)) * 100
            print(f"{percentage_complete:.0f}% complete")
        
        # Extract saws values at this specific time-step
        # flattening for use in griddata
        # (we can use .isel for this since we know the blk and saws times are aligned)
        ds_saws_now = ds_saws.isel(time=t)
        u_saws_now = ds_saws_now['10u'].values.flatten()
        v_saws_now = ds_saws_now['10v'].values.flatten()
        t2_saws_now = ds_saws_now['2t'].values.flatten() - 273.15 # converting K to degrees celcius
        rhum_saws_now = ds_saws_now['r'].values.flatten() / 100 # converting percent to fraction
        
        # Perform interpolation
        # 'nearest' method speeds it up a bit
        # but I guess 'linear' is better
        u_saws_interp[t,:,:] = griddata(source_points, u_saws_now, (lon_rho, lat_rho), method=interp_method) 
        v_saws_interp[t,:,:] = griddata(source_points, v_saws_now, (lon_rho, lat_rho), method=interp_method) 
        t2_saws_interp[t,:,:] = griddata(source_points, t2_saws_now, (lon_rho, lat_rho), method=interp_method)
        rhum_saws_interp[t,:,:] = griddata(source_points, rhum_saws_now, (lon_rho, lat_rho), method=interp_method)              
        
    # calculate the wspd from the wasa components (on the rho grid)
    spd_saws_interp = np.hypot(u_saws_interp, v_saws_interp)
    
    # uwnd and vwnd are grid aligned, so we need to rotate vectors
    # (based on official croco_tools interp_ERA5.m)
    cos_a = np.cos(angle.values)
    sin_a = np.sin(angle.values)
    # have to use the data on the rho grid for the rotation calc 
    u_out = u_saws_interp*cos_a + v_saws_interp*sin_a
    v_out = v_saws_interp*cos_a - u_saws_interp*sin_a
    # but get onto the u,v grids to write to file
    u_out = post.rho2u(u_out)
    v_out = post.rho2v(v_out)
    
    # Update the interpolated saws data within ds_blk
    ds_blk['uwnd'] = (('bulk_time', 'eta_u', 'xi_u'), u_out)
    ds_blk['vwnd'] = (('bulk_time', 'eta_v', 'xi_v'), v_out)
    ds_blk['wspd'] = (('bulk_time', 'eta_rho', 'xi_rho'), spd_saws_interp)
    ds_blk['tair'] = (('bulk_time', 'eta_rho', 'xi_rho'), t2_saws_interp)
    ds_blk['rhum'] = (('bulk_time', 'eta_rho', 'xi_rho'), rhum_saws_interp)
    
    # I'm pretty sure we don't need sustr and svstr in blk files as this is computed online
    # So I'm removing these variables to be sure
    # leaving them hanging around seems a little messy since we're replacing uwnd and vwnd
    # while we're at it, radlw is also not used (radlw_in is), so remove that as well
    ds_blk = ds_blk.drop_vars(["sustr", "svstr", "radlw"])
    
    # write float32 instead of float64 in an attempt to save space
    encoding = {
        "tair": {"dtype": "float32"},
        "rhum": {"dtype": "float32"},
        "prate": {"dtype": "float32"},
        "wspd": {"dtype": "float32"},
        "radlw_in": {"dtype": "float32"},
        "radsw": {"dtype": "float32"},
        "uwnd": {"dtype": "float32"},
        "vwnd": {"dtype": "float32"},
        "bulk_time": {"dtype": "float32"},
    }
    
    # write the new blk file
    ds_blk.to_netcdf(croco_blk_out,
                      encoding=encoding)
    
    ds_blk.close()
    
    ds_saws.close()

def reformat_saws_atm(saws_dir,backup_dir,out_dir,run_date,hdays,Yorig):
    '''
    Convert the SAWS UM files into a format which can be ingested by CROCO
    using the ONLINE cpp key for online interpolation of the surface forcing
    (we will use the default 'CFSR' file format)
    
    Parameters
    ----------
    saws_dir      : the directory where the saws um files are located
    backup_dir    : directory containing already reformatted data, used for variables not provided by SAWS
    output_dir    : the directory where the croco-friendly files will be saved
    run_date      : a datetime.datetime object of the forecast initialisation time
    hdays         : integer of the hindcast days component of the croco run (used to identify relevatn saws files)
    Yorig         : the origin year used in setting up the croco times
    '''
    
    # Get a list of all .nc files in the saws directory...
    files = sorted(glob.glob(f"{saws_dir}/SAWSUM_SA4-*.nc"))
    # ... and reverse the order to give priority to the latest files
    files.reverse()
    
    # check if latest saws file isn't older than 12 hours since run_date
    saws_latest_datetime = datetime.strptime(files[0][-13:-4], "%Y%m%d%H")
    saws_latest_datetime_allowed = run_date - timedelta(hours=12)
    if saws_latest_datetime < saws_latest_datetime_allowed:
        # end with error
        raise ValueError("The latest SAWS file is more than 12 hours older than the run_date - aborting SAWS forced run")

    # Filter the files list to include only files after the CROCO start datetime
    cutoff_datetime = run_date - timedelta(days=(hdays+1)) # extend by a day to make sure we cover our CROCO run
    filtered_files = [
        file for file in files
        if datetime.strptime(file[-13:-4], "%Y%m%d%H") >= cutoff_datetime
    ]

    # get the backup (currently gfs) grid using a template file (could be any variable)
    # this is used to subset to SAWS data, since it covers a much bigger domain than we need
    grid_file = os.path.join(backup_dir, "Temperature_height_above_ground_Y9999M1.nc")
    ds_grid = xr.open_dataset(grid_file)
    lon_grid = ds_grid.lon.values
    lon_lims = slice(lon_grid.min(),lon_grid.max())
    lat_grid = ds_grid.lat.values
    lat_lims = slice(lat_grid.min(),lat_grid.max())
    ds_grid.close()

    # only some variables are provided to us by SAWS
    # so we have a "existsSAWS" flag to indicate which variables we have
    # the ones we don't must be interpolated from another source (e.g. GFS or other) onto the SAWS grid for the CROCO model to be able to run
    variables = {
            "Temperature_height_above_ground": {
                "shortName": "2t",
                "existsSAWS": "true",
            },
            "U-component_of_wind": {
                "shortName": "10u",
                "existsSAWS": "true",
            },
            "V-component_of_wind": {
                "shortName": "10v",
                "existsSAWS": "true",
            },
            "Specific_humidity": {
                "shortName": "2sh",
                "existsSAWS": "false",
            },
            "Precipitation_rate": {
                "shortName": "prate",
                "existsSAWS": "false",
            },
            "Downward_Short-Wave_Rad_Flux_surface": {
                "shortName": "dswrf",
                "existsSAWS": "false",
            },
            "Upward_Short-Wave_Rad_Flux_surface": {
                "shortName": "uswrf",
                "existsSAWS": "false",
            },
            "Downward_Long-Wave_Rad_Flux": {
                "shortName": "dlwrf",
                "existsSAWS": "false",
            },
            "Upward_Long-Wave_Rad_Flux_surface": {
                "shortName": "ulwrf",
                "existsSAWS": "false",
            },
        }
    
    for var in variables:
        
        # get an xarray dataset for this variable
        print('working on '+var)
        var_dict = variables[var]
        
        if var_dict['existsSAWS'] == 'true': # we can use saws data for this variable

            # get an xarray dataset for the SAWS data for this variable
            # giving preference to more recent files where data overlap
            da = None
            for file in filtered_files:
                ds_file = xr.open_dataset(file)
                # subset the data spatially
                ds_file = ds_file.sel(lon=lon_lims, lat=lat_lims)
                # extract the variable
                da_file = ds_file[var_dict['shortName']].squeeze()
                # the concatenation step below seems to take quite long (it's a big domain!) 
                # but pre-chunking the lon and lat dimensions seems to speed it up a bit
                da_file = da_file.chunk({'time': 1, 'lat': 100, 'lon': 100}) 
                if da is None:
                    da = da_file
                else:
                    # extract the non-overlapping data
                    da_file = da_file.sel(time=slice(da_file.time[0],da.time[0]-1))
                    # and combine with the full dataset
                    da = xr.concat([da_file, da], dim='time')
                da_file.close()
                ds_file.close()
            
            # make a dataset from the dataarray
            ds = xr.Dataset({var: da})
            
            # we need to convert time to days since Yorig!
            # Reference date
            reference_date = np.datetime64(str(Yorig)+'-01-01T00:00:00')
            # time_in_ns = ds['time'].astype('datetime64[ns]')
            # Convert the time dimension to days since the reference date
            ds['time'] = (ds['time'].astype('datetime64[ns]') - reference_date) / np.timedelta64(1, 'D')
            # Set the units attribute for the time coordinate
            ds['time'].attrs['units'] = 'days since 1-Jan-'+str(Yorig)+' 00:00:00'
            
        else: # we have to interpolate previously formatted data onto the SAWS grid to have the full compliment of variables

            # get the backup data which we need to interpolate onto the saws grid
            backup_file = os.path.join(backup_dir, var+"_Y9999M1.nc")
            ds_backup = xr.open_dataset(backup_file)
            
            # start by getting a template saws file for interpolating 
            # we're using the Temperature_height_above_ground file, which would have already been processed
            template_file = os.path.join(out_dir, "Temperature_height_above_ground_Y9999M1.nc")
            ds_template = xr.open_dataset(template_file)
            
            # do the interpolation
            # this succinct approach works because our files have the exact same dimensions
            # ds = ds_backup.interp_like(ds_template, method='linear')
            # Or we can specify the interpolation explicitly (not sure if it makes a difference)
            ds = ds_backup.interp(
                lon=ds_template.lon, 
                lat=ds_template.lat, 
                time=ds_template.time, 
                method="linear"
            )
            
            ds_backup.close()
            ds_template.close()
            
        # write the nc file
        # the ONLINE cppkey is designed for use with monly interannual simulations
        # where the year and month of the file name is appended to the end of the file
        # Since we are using this option with forecasts we'll just put dummy values
        # for the year and month, and put these values in the *.in file making it 
        # something we don't have to handle separately
        print('writing the file...')
        fname_out = os.path.join(out_dir,var+"_Y9999M1.nc")
        if os.path.exists(fname_out):
            os.remove(fname_out)
        ds.astype('float32').to_netcdf(fname_out) # writing as float32 halves the file size!
        ds.close()

def reformat_gfs_atm(gfs_dir,out_dir,Yorig):
    '''
    Convert the GFS atmospheric forecast grb files downloaded by the cli.py function download_gfs_atm
    and convert the data into nc files in a format which can be ingested by CROCO
    using the ONLINE cpp key for online interpolation of the surface forcing
    (we will use the default 'CFSR' file format)
    '''
    
    # Path to the directory containing your GRIB files
    gfs_files = os.path.join(gfs_dir,"*.grb")
    
    # List all GRIB files
    file_paths = sorted(glob.glob(gfs_files))
    
    def open_grib_file(file_path,var_dict):
        return xr.open_dataarray(
            file_path,
            engine='cfgrib',
            filter_by_keys={'shortName': var_dict['shortName'],
                            'stepType': var_dict['stepType'],
                            })

    # (you can see the vars in the files using e.g.
    # grib_ls -P shortName,typeOfLevel,level 2024080106_f001.grb)
    
    variables = {
            "Temperature_height_above_ground": {
                "shortName": "2t",
                "stepType": "instant",
            },
            "Specific_humidity": {
                "shortName": "2sh",
                "stepType": "instant",
            },
            "Precipitation_rate": {
                "shortName": "prate",
                "stepType": "instant",
            },
            "Downward_Short-Wave_Rad_Flux_surface": {
                "shortName": "dswrf",
                "stepType": "avg",
            },
            "Upward_Short-Wave_Rad_Flux_surface": {
                "shortName": "uswrf",
                "stepType": "avg",
            },
            "Downward_Long-Wave_Rad_Flux": {
                "shortName": "dlwrf",
                "stepType": "avg",
            },
            "Upward_Long-Wave_Rad_Flux_surface": {
                "shortName": "ulwrf",
                "stepType": "avg",
            },
            "U-component_of_wind": {
                "shortName": "10u",
                "stepType": "instant",
            },
            "V-component_of_wind": {
                "shortName": "10v",
                "stepType": "instant",
            },
            # will need to add patm here
        }
    
    for var in variables:
        
        # get an xarray dataset for this variable
        print('working on '+var)
        var_dict = variables[var]
        datasets = [open_grib_file(fp,var_dict) for fp in file_paths]
        da = xr.concat(datasets, dim='valid_time')
        
        if var_dict['stepType']=='avg':
            # A bunch of variables are (rather annoyingly) written out as the
            # accumulated average over each 6 hour forecast period.
            # We want to convert these accumulated averages into individual one-hour averages
            #
            # see FAQ "How can the individual one-hour averages be computed" from https://rda.ucar.edu/datasets/ds093.0/#docs/FAQs_6hrly.html
            # Excerpt from there:
            # You can compute the one-hour average (X) ending at hour N by using the N-hour average (a) and the (N-1)-hour average (b) as follows:
            # X = N*a - (N-1)*b
            # So if you want the 1-hour Average for the period initial+3 to initial+4 (X), you would use the 4-hour Average (initial+0 to initial+4) as (a) and the 3-hour Average (initial+0 to initial+3) as (b) as follows:
            # X = 4*a - 3*b
                
            # start by extracting the forecast hour for each time-step (can be derived from the 'step' variable)
            # frcst = xr.DataArray(da.step / np.timedelta64(1, 'h'), dims='valid_time')
            frcst = da.step.values / np.timedelta64(1, 'h')
            
            # the averaging period for these variables is (again rather annoyingly) reset every 6 hours, 
            # even if we are looking at forecast hours greater than 6
            # so here we compute the hour within in the 6 hour forecast cycle 
            # (frcst_ave will range from 1-6 by definition)
            frcst_ave = np.mod(frcst - 1, 6) + 1
            
            # get arrays for doing the calcs and updating the output array
            data = da.values
            data_output = data.copy()
            
            for i in range(len(frcst)):
                if frcst_ave[i]>1:
                    # only do the conversion for forecast hours greater than 1
                    if frcst[i]<=120:
                        data_output[i,::] = data[i,::]*frcst_ave[i]-data[i-1,::]*frcst_ave[i-1]
                    else:
                        # after 120 hrs the forecasts are provided at three hourly intervals
                        # so each of these represents the accumulated hourly averages over 3 and 6 hours
                        # so we have to handle this separately to get back to hourly averages
                        # 
                        if frcst_ave[i]==3:
                            data_output[i,::] = data[i,::]/3
                        else: # frcst_ave[i]==6
                            data_output[i,::] = (data[i,::]*frcst_ave[i]-data[i-1,::]*frcst_ave[i-1])/3
            # Create a new DataArray with the same dimensions and coordinates as the original
            # but with the updated output
            da = xr.DataArray(
                data_output,
                dims=da.dims,
                coords=da.coords,
                attrs=da.attrs
            )
        
        # rename the dimensions to match what is expected by CROCO ONLINE option
        da = da.drop_vars(['time','step'])
        da = da.rename({'longitude': 'lon', 'latitude': 'lat', 'valid_time': 'time'})
        
        # handle any 3 hourly time-steps at the end of the data 
        # by interpolating onto an hourly time axis (I'm not totally sure this is needed but no harm done)
        time_equidistant = pd.date_range(start=da.time.min().values, end=da.time.max().values, freq='h')
        da = da.interp(time=time_equidistant)
        
        # make a dataset from the dataarray
        ds = xr.Dataset({var: da})
        
        # we need to convert time to days since Yorig!
        # Reference date
        reference_date = np.datetime64(str(Yorig)+'-01-01T00:00:00')
        # time_in_ns = ds['time'].astype('datetime64[ns]')
        # Convert the time dimension to days since the reference date
        ds['time'] = (ds['time'].astype('datetime64[ns]') - reference_date) / np.timedelta64(1, 'D')
        # Set the units attribute for the time coordinate
        ds['time'].attrs['units'] = 'days since 1-Jan-'+str(Yorig)+' 00:00:00'
        
        # write the nc file
        # the ONLINE cppkey is designed for use with monly interannual simulations
        # where the year and month of the file name is appended to the end of the file
        # since we are using this option with forecasts we'll just put dummy values
        # for the year and month, and put these values in the *.in file making it 
        # something we don't have to handle separately
        fname_out = os.path.join(out_dir,var+"_Y9999M1.nc")
        ds.to_netcdf(fname_out)
        ds.close()
        
        # Delete all temporary files created during reading the .grb files (not sure what these are but they're not needed)
        for gbx9_file in glob.glob(os.path.join(gfs_dir, '*.gbx9')):
            os.remove(gbx9_file)
        for ncx_file in glob.glob(os.path.join(gfs_dir, '*.ncx')):
            os.remove(ncx_file)
        for idx_file in glob.glob(os.path.join(gfs_dir, '*.idx')):
            os.remove(idx_file)
        
def make_ini_fcst(input_file,output_dir,run_date,hdays):
    '''
    Make CROCO initial file for SAOMISANA
    
    output_dir - the directory where the ini file will be saved
                 NB - there needs to be a crocotools_param.py file in this directory which includes all the configurable parameters
    run_date   - the time when the operational run was initialised, as a datetime.datetime object
    hdays      - the number of hindcast days used in the operational run (time of the ini file will be run_date - hdays)
    
    '''  
    # --------
    # imports
    # --------
    #
    sys.path.append(output_dir)
    import crocotools_param as params
    
    # edit ini_filename to add starting date
    ini_filename = params.ini_filename.replace('.nc', run_date.strftime('_%Y%m%d_%H.nc'))

    # Load croco_grd
    crocogrd = Croco.CROCO_grd(params.croco_grd, params.sigma_params)
    #--- Load input (restricted to croco_grd) ----------------------------  
    
    multi_files=False
    print(params.inputdata)
    print(input_file)
    print(crocogrd)
    inpdat=Inp.getdata(params.inputdata,input_file,crocogrd,multi_files,params.tracers)
    #inpdat=Inp.getdata(params.inputdata,input_file,crocogrd,multi_files,params.tracers,bdy=[params.obc_dict,params.cycle_bry])
    print(inpdat)
    # --- Create the initial file -----------------------------------------

    
    Croco.CROCO.create_ini_nc(None,''.join((params.croco_dir + ini_filename)),crocogrd,
                              tracers=params.tracers)

    # --- Handle initial time ---------------------------------------------

    ini_date_num = run_date - timedelta(days=hdays)
    ini_date_num = ini_date_num.timestamp()
    
    day_zero_num = datetime(int(params.Yorig), int(params.Morig), int(params.Dorig))
    
    day_zero_num = day_zero_num.timestamp()
    
    tstart=0
    
    scrumt = ini_date_num

    oceant = ini_date_num
    
    tend=0.

    # time index to use in the file
    
    tndx = 0

   #  --- Compute and save variables on CROCO grid ---------------

    for vars in ['ssh','tracers','velocity']:
        print('\nProcessing *%s*' %vars)
        nc=netcdf.Dataset(params.croco_dir+ini_filename, 'a')
        if vars == 'ssh' :
            (zeta,NzGood) = interp_tools.interp_tracers(inpdat,vars,-1,crocogrd,tndx,tndx)
            nc.variables['zeta'][0,:,:] = zeta*crocogrd.maskr
            nc.Input_data_type=params.inputdata
            nc.variables['ocean_time'][:] = oceant
            nc.variables['scrum_time'][:] = scrumt
            nc.variables['scrum_time'].units='seconds since %s-%s-%s 00:00:00' %(params.Yorig,params.Morig,params.Dorig)
            nc.variables['tstart'][:] = tstart
            nc.variables['tend'][:] = tend
            z_rho = crocogrd.scoord2z_r(zeta=zeta)
            z_w   = crocogrd.scoord2z_w(zeta=zeta)
        elif vars == 'tracers':
            for tra in params.tracers:
                print(f'\nIn tracers processing {tra}')
                trac_3d= interp_tools.interp(inpdat,tra,params.Nzgoodmin,z_rho,crocogrd,tndx,tndx)
                nc.variables[tra][0,:,:,:] = trac_3d*crocogrd.mask3d()

        elif vars == 'velocity':

            cosa=np.cos(crocogrd.angle)
            sina=np.sin(crocogrd.angle)

            [u,v,ubar,vbar]=interp_tools.interp_uv(inpdat,params.Nzgoodmin,z_rho,cosa,sina,crocogrd,tndx,tndx)
              
            conserv=1  # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport          
            if conserv == 1:
                (ubar_croco,h0)=sig_tools.vintegr(u,grd_tools.rho2u(z_w),grd_tools.rho2u(z_rho),np.nan,np.nan)/grd_tools.rho2u(crocogrd.h)
                (vbar_croco,h0)=sig_tools.vintegr(v,grd_tools.rho2v(z_w),grd_tools.rho2v(z_rho),np.nan,np.nan)/grd_tools.rho2v(crocogrd.h)

                u = u - ubar_croco ; u = u + np.tile(ubar,(z_rho.shape[0],1,1))
                v = v - vbar_croco ; v = v + np.tile(vbar,(z_rho.shape[0],1,1))
           
            nc.variables['u'][0,:,:,:] = u *crocogrd.umask3d()
            nc.variables['v'][0,:,:,:] = v * crocogrd.vmask3d()
            nc.variables['ubar'][0,:,:] = ubar *crocogrd.umask
            nc.variables['vbar'][0,:,:] = vbar * crocogrd.vmask

   
    nc.close()
    
    print('')
    print(' Initial file created ')
    print(' Path to file is ', params.croco_dir + ini_filename)
    print('')

def make_bry_fcst(input_file,output_dir,run_date,hdays):
    '''
    
    Make CROCO boundary file for SAOMISANA
    
    '''
    # --------
    # imports
    # --------
    #
    sys.path.append(output_dir)
    import crocotools_param as params
    
    # --- Make filename --------------------------------------------------

    bry_filename = params.bry_filename.replace('.nc', run_date.strftime('_%Y%m%d_%H.nc'))

    # --- Load croco_grd --------------------------------------------------

    crocogrd = Croco.CROCO_grd(params.croco_grd, params.sigma_params)

    # --- Initialize boundary vars ----------------------------------------

    crocogrd.WEST_grid()
    crocogrd.EAST_grid()
    crocogrd.SOUTH_grid()
    crocogrd.NORTH_grid()

    # --- Initialize input data class -------------------------------------
    
    multi_files=False
    inpdat = Inp.getdata(params.inputdata,input_file,crocogrd,multi_files,
                         params.tracers,
                         bdy=[params.obc_dict,params.cycle_bry])

    Croco.CROCO.create_bry_nc(None,params.croco_dir + bry_filename,crocogrd,params.obc_dict,params.cycle_bry,tracers=params.tracers)

    # --- Handle bry_time --------------------------------------------

    nc=netcdf.Dataset(params.croco_dir + bry_filename, 'a')
    
    bry_start_date = plt.date2num(run_date - timedelta(days=hdays))
    bry_end_date = plt.date2num(run_date + timedelta(days=1))

    # Load full time dataset
    time = plt.date2num(inpdat.ncglo['time'].values)
    
    # find index for the time range
    ind = np.where((time>bry_start_date) & (time<=bry_end_date))
    [dtmin,dtmax] = np.min(ind),np.max(ind)
    bry_time = time[int(dtmin):int(dtmax)]
    nc.Input_data_type=params.inputdata
    nc.variables['bry_time'].cycle=params.cycle_bry
    nc.variables['bry_time'][:]=bry_time
    
    if params.cycle_bry==0:
        nc.variables['bry_time'].units='days since %s-01-01 00:00:00' %(params.Yorig)
        
    # --- Loop on boundaries ------------------------------------------
    prev=1
    nxt=1

    if len(params.tracers) == 0:
        var_loop = ['ssh','velocity']
    else:
        var_loop = ['ssh','tracers','velocity']

    for boundary, is_open in zip(params.obc_dict.keys(), params.obc_dict.values()):
        if is_open:
            for vars in var_loop:
                print('\n     Processing *%s* for %sern boundary' %(vars, boundary))
                print('     ------------------------------------------')
                if vars == 'ssh':
                    (zeta,NzGood) = interp_tools.interp_tracers(inpdat,vars,-1,crocogrd,dtmin,dtmax,prev,nxt,boundary[0].upper())
                    z_rho = crocogrd.scoord2z_r(zeta=zeta,bdy="_"+boundary)
                    z_w   = crocogrd.scoord2z_w(zeta=zeta,bdy="_"+boundary)

                elif vars == 'tracers':
                    trac_dict = dict()
                    for trc in params.tracers:
                        print(f'\nIn tracers processing {trc}')
                        trac_dict[trc] = interp_tools.interp(inpdat,trc,params.Nzgoodmin,z_rho,crocogrd,dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())

                elif vars == 'velocity':
                    cosa=np.cos(eval(''.join(('crocogrd.angle_',boundary))) )
                    sina=np.sin(eval(''.join(('crocogrd.angle_',boundary))) )

                    [u,v,ubar,vbar]=interp_tools.interp_uv(inpdat,params.Nzgoodmin,z_rho,cosa,sina,crocogrd,dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())

                    conserv=1  # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport

                    if conserv == 1:
                        ubar_croco=sig_tools.vintegr4D(u,grd_tools.rho2u(z_w),grd_tools.rho2u(z_rho),np.nan,np.nan)[0]/grd_tools.rho2u(eval(''.join(('crocogrd.h_'+boundary))))
                        vbar_croco=sig_tools.vintegr4D(v,grd_tools.rho2v(z_w),grd_tools.rho2v(z_rho),np.nan,np.nan)[0]/grd_tools.rho2v(eval(''.join(('crocogrd.h_'+boundary))))

                        u = u - np.tile(ubar_croco[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                        u = u + np.tile(ubar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))

                        v = v - np.tile(vbar_croco[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                        v = v + np.tile(vbar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))

    # --- Saving in netcdf ------------------------------------------------

            print('\nSaving %sern boundary in Netcdf' % boundary)
            print('----------------------------------')

             # handle indices (as 2 points where taken next to bdy)
            if str(boundary) == 'west' and is_open:
                indices3D="[:,:,:,0]" # T,N,J,i=0
                indices2D="[:,:,0]"   # T,J,i=0
            elif str(boundary) == 'east' and is_open:
                indices3D="[:,:,:,-1]" # T,N,J,i=last
                indices2D="[:,:,-1]"   # T,J,i=last
            elif str(boundary) == 'south' and is_open:
                indices3D="[:,:,0,:]" # T,N,j=0,I
                indices2D="[:,0,:]"   # T,j=0,I
            elif str(boundary) == 'north' and is_open:
                indices3D="[:,:,-1,:]" # T,N,j=last,I
                indices2D="[:,-1,:]"   # T,j=last,I

            mask_zet = np.tile(eval(''.join(('crocogrd.maskr_',boundary))),[zeta.shape[0],1,1])
            if "velocity" in var_loop:
                mask_u   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],u.shape[1],1,1])
                mask_v   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[u.shape[0],u.shape[1],1,1])
                mask_ubar   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],1,1])
                mask_vbar   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[v.shape[0],1,1])

            nc.variables['zeta_'+str(boundary)][:]=eval(''.join(('zeta',indices2D)))*eval(''.join(('mask_zet',indices2D)))
            nc.variables['u_'+str(boundary)][:]   =eval(''.join(('u',indices3D)))*eval(''.join(('mask_u',indices3D)))
            nc.variables['v_'+str(boundary)][:]   =eval(''.join(('v',indices3D)))*eval(''.join(('mask_v',indices3D)))
            nc.variables['ubar_'+str(boundary)][:]=eval(''.join(('ubar',indices2D)))*eval(''.join(('mask_ubar',indices2D)))
            nc.variables['vbar_'+str(boundary)][:]=eval(''.join(('vbar',indices2D)))*eval(''.join(('mask_vbar',indices2D)))

            if 'tracers' in var_loop:
                for varname, value in zip(trac_dict.keys(), trac_dict.values()):
                    mask_tra = np.tile(eval(''.join(('crocogrd.maskr_',boundary))),[value.shape[0],value.shape[1],1,1])
                    nc.variables[f"{varname}_{boundary}"][:] = eval(f'value{indices3D}')*eval(''.join(('mask_tra',indices3D)))

    nc.close()

    print('')
    print(' Boundary file created ')
    print(' Path to file is ', params.croco_dir + bry_filename)
    print('')


if __name__ == '__main__':
    #input_file = '/home/rautenbach/SOMISANA/croco/somisana-croco/DATA/MERCATOR/MERCATOR.nc'
    input_file = '/home/rautenbach/SOMISANA/croco/somisana-croco/DATA/HYCOM/HYCOM.nc'
    output_dir='/home/rautenbach/SOMISANA/croco/somisana-croco/configs/sa_west_02/croco_v1.3.1/MERCATOR/'
    run_date = datetime(2024,8,18,0,0,0)
    hdays = 5
    make_bry_fcst(input_file,output_dir,run_date,hdays)
    #make_ini_fcst(input_file,output_dir,run_date,hdays)
        
