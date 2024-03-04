import xarray as xr
from datetime import datetime, timedelta
import numpy as np
import os
import fnmatch
from scipy.interpolate import griddata
import crocotools_py.postprocess as post

def fill_blk(croco_grd,croco_blk_file_in,croco_blk_file_out):
    '''
    This is a little hack function to help us getting around the specific problem
    of having a block of missing data in the WASA3 zarr files for about 13 days in April 2010
    So we're just using this function to interpolate over this block.
    Maybe this function will be usefull in other cases too?
    
    I'm specifically not doing this automatically in make_WASA3_from_blk()
    as the user should really know exactly where the missing data are and 
    first see if this hack is appropriate
    
    If the block of missing data is too big, maybe filling with ERA5 would be more appropriate
    Or maybe we have to make another plan
    
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
            u_blk_filled[t,:,:] = griddata(source_points_u[idx,:], u_blk_now[idx], (lon_u, lat_u), method='linear')
        if np.any(np.isnan(v_blk_now)):
            # get the indices for which we have valid data
            idx = np.logical_not(np.isnan(v_blk_now))
            # Perform interpolation using valid data
            v_blk_filled[t,:,:] = griddata(source_points_v[idx,:], v_blk_now[idx], (lon_v, lat_v), method='linear')
        
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
    
def make_WASA3_from_blk(wasa_grid, 
                 wasa_zarr_dir, 
                 croco_grd,
                 croco_blk_dir,
                 out_wasa_dir,
                 ref_date,
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
    lon_wasa = ds_wasa_grd.XLONG_M.squeeze()
    lat_wasa = ds_wasa_grd.XLAT_M.squeeze()
    
    # get the croco grid variables
    ds_croco_grd=xr.open_dataset(croco_grd)
    lon_rho=ds_croco_grd.lon_rho
    lat_rho=ds_croco_grd.lat_rho
    angle=ds_croco_grd.angle
    ds_croco_grd.close()
    
    # Subset WASA spatially to speed up interpolation later
    # this helps a lot so worth the mess!
    def find_nearest_wasa_point(lon_wasa, lat_wasa, lon_in, lat_in):
        # Sub-function to find the nearest indices of the wasa grid to a specific lon, lat coordinate
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
    # get the nearest wasa grid indices corresponding to the corners of the croco grid
    dl = 0.1 # buffer    
    # bl = bottom left, br = bottom right, tr = top right, tl = top_left
    j_bl, i_bl = find_nearest_wasa_point(lon_wasa, lat_wasa, lon_rho[0,0] - dl, lat_rho[0,0] - dl)
    j_br, i_br = find_nearest_wasa_point(lon_wasa, lat_wasa, lon_rho[0,-1] + dl, lat_rho[0,-1] - dl)
    j_tr, i_tr = find_nearest_wasa_point(lon_wasa, lat_wasa, lon_rho[-1,-1] + dl, lat_rho[-1,-1] + dl)
    j_tl, i_tl = find_nearest_wasa_point(lon_wasa, lat_wasa, lon_rho[-1,0] - dl, lat_rho[-1,0] + dl)
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
    
    # extract the wasa grid variables we need from this spatial subset
    lon_wasa = ds_wasa_grd.XLONG_M.squeeze()
    lat_wasa = ds_wasa_grd.XLAT_M.squeeze() 
    cosa = ds_wasa_grd.COSALPHA.squeeze()
    sina = ds_wasa_grd.SINALPHA.squeeze()
    # flatten lon,lat arays for use in interpolation later
    source_points = np.column_stack((np.array(lon_wasa).flatten(), np.array(lat_wasa).flatten()))

    # Loop through croco blk files
    croco_blk_files = sorted(fnmatch.filter(os.listdir(croco_blk_dir), '*blk*.nc*'))
    for croco_blk_file in croco_blk_files:
            
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
            u = ds_wasa_now.U
            v = ds_wasa_now.V
            u_wasa = u * cosa - v * sina
            v_wasa = v * cosa + u * sina
            
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
            ds_blk.drop_vars(["sustr", "svstr", "radlw"])
            
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
        
    ds_wasa.close()
    ds_wasa_grd.close()

    
    
