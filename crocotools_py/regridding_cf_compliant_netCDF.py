import crocotools_py.postprocess as post
import numpy as np
import xarray as xr
import os,sys
from datetime import datetime,timedelta
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from scipy.spatial import Delaunay

def get_tri_coef(lonold_vals, latold_vals, lonnew_vals, latnew_vals, verbose=0):
    """
    Compute Delaunay triangulation coefficients for xarray DataArrays.
    
    Inputs:
    - lonold_vals, latold_vals: 2D numpy arrays from original xarray grids
    - lonnew_vals, latnew_vals: 2D numpy arrays from target xarray grids
    
    Returns:
    - elem: Indices of triangulation points used for interpolation
    - coef: Interpolation coefficients
    """
    
    # Flatten 2D grids into 1D arrays
    Xp = np.array([lonold_vals.ravel(), latold_vals.ravel()]).T
    Xc = np.array([lonnew_vals.ravel(), latnew_vals.ravel()]).T

    # Compute Delaunay triangulation
    if verbose: print("Computing Delaunay triangulation...")
    tri = Delaunay(Xp)

    # Find which triangle each new point falls into
    simplex_indices = tri.find_simplex(Xc)

    # Create a mask for points inside the triangulation
    valid_points = simplex_indices >= 0

    # Initialize empty outputs with NaNs
    elem = np.full((Xc.shape[0], 3), -1, dtype=int)  # Default to -1 for out-of-bounds points
    coef = np.full((Xc.shape[0], 3), np.nan)

    # Process only valid points
    if np.any(valid_points):
        valid_Xc = Xc[valid_points]
        valid_simplices = tri.simplices[simplex_indices[valid_points]]

        points = tri.points[valid_simplices]
        p = np.zeros((valid_Xc.shape[0], 3))

        # Compute barycentric coordinates
        for i in range(valid_Xc.shape[0]):
            A = np.append(np.ones((3, 1)), points[i], axis=1)
            B = np.append(1., valid_Xc[i])
            p[i, :] = np.linalg.lstsq(A.T, B.T, rcond=None)[0]

        elem[valid_points] = valid_simplices
        coef[valid_points] = p

    # Reshape results to match new grid shape
    elem = elem.reshape(lonnew_vals.shape + (3,))
    coef = coef.reshape(lonnew_vals.shape + (3,))

    return elem, coef

def horiz_interp_delaunay(lonold, latold, varold, lonnew, latnew, elem=None, coef=None):
    """
    Perform horizontal interpolation using Delaunay triangulation for xarray DataArray objects.

    Parameters:
    - lonold, latold: xarray.DataArray (2D) - Original longitude and latitude grids
    - varold: xarray.DataArray (2D) - Original variable to be interpolated
    - lonnew, latnew: xarray.DataArray (2D) - New longitude and latitude grids
    - elem, coef: (Optional) Precomputed triangulation elements and coefficients for efficiency

    Returns:
    - varnew: xarray.DataArray (2D) - Interpolated variable on the new grid
    - (elem, coef): Only returned when computing coefficients
    """

    # Ensure input is xarray DataArray and extract numerical values
    lonold_vals, latold_vals, varold_vals = lonold.values, latold.values, varold.values
    lonnew_vals, latnew_vals = lonnew.values, latnew.values

    if (elem is None) | (coef is None):
        # Compute interpolation coefficients
        elem, coef = get_tri_coef(lonold_vals, latold_vals, lonnew_vals, latnew_vals)

        # Normalize coefficients
        coef /= np.sum(coef, axis=2, keepdims=True)

        # Perform interpolation
        varnew_vals = np.sum(coef * varold.values.ravel()[elem], axis=2)

        # Return coefficients for future reuse
        return elem, coef, xr.DataArray(varnew_vals, coords=lonnew.coords, dims=lonnew.dims, attrs=varold.attrs)

    else:
        # Use precomputed interpolation coefficients
        varnew_vals = np.sum(coef * varold.values.ravel()[elem], axis=2)
        return xr.DataArray(varnew_vals, coords=lonnew.coords, dims=lonnew.dims, attrs=varold.attrs)

def write_attrs(info):
    attrs={
        "title"          : info['Description/Guidance'][3],
        "institution"    : info['Description/Guidance'][6],
        "source"         : "CROCO model",
        "history"        : "Created " + str(datetime.strftime(datetime.now(),"%Y-%m-%d %H:%M:%S")),
        "reference"      : "DOI ###",
        "conventions"    : "CF-1.8",
        "author"         : info['Description/Guidance'][4],
        "contact person" : info['Description/Guidance'][5],
        "summery"        : info['Description/Guidance'][10],
        "keywords"       : info['Description/Guidance'][8],
    }
    return attrs

def write_coords(ds, regrid_ver, ref_date=None):
    if ref_date is None:
        ref_date = datetime(2000, 1, 1)

    if (regrid_ver == 1) | (regrid_ver == 2):
        coords = {
            "xi_rho": ("xi_rho", ds.xi_rho.data, {
                "standard_name": "projection_x_coordinate",
                "long_name": "xi_rho grid index",
                "units": "1",
                "axis": "X"
            }),
            "eta_rho": ("eta_rho", ds.eta_rho.data, {
                "standard_name": "projection_y_coordinate",
                "long_name": "eta_rho grid index",
                "units": "1",
                "axis": "Y"
            }),
        }
    elif regrid_ver == 3:
        coords = {
            "longitude": ("longitude", ds.longitude.data, {
                "standard_name": "longitude",
                "long_name": "Longitude",
                "units": "degrees_east",
                "axis": "X"
            }),
            "latitude": ("latitude", ds.latitude.data, {
                "standard_name": "latitude",
                "long_name": "Latitude",
                "units": "degrees_north",
                "axis": "Y"
            }),
        }
    else:
        print('')
        print(f'Regrid version is {regrid_ver}.')
        print(f'Valid versions are: 1, 2 and 3')
        print('')

    # Conditionally add the vertical coordinate before 'time'
    if regrid_ver == 1:
        coords["s_rho"] = ("s_rho", ds.s_rho.data, {
            "standard_name": "ocean_sigma_coordinate",
            "long_name": "terrain-following vertical coordinate",
            "units": "1",
            "positive": "up",
            "axis": "Z"
        })
    else:
        coords["depth"] = ("depth", ds.depths.data, {
            "standard_name": "depth",
            "long_name": "Fixed depth levels",
            "units": "m",
            "positive": "up",
            "axis": "Z"
        })

    # Ensure 'time' is always last
    coords["time"] = ("time", ((pd.to_datetime(ds.time.values) - ref_date) / pd.Timedelta(days=1)).astype(np.float32), {
        "standard_name": "time",
        "long_name": "Time",
        "units": f"days since {ref_date}",
        "calendar": "standard",
        "axis": "T",
    })

    return coords

def write_vars(ds, regrid_ver, fill_value=None, ref_date=None):
    
    if fill_value is None:
        fill_value = -9999.0
        
    if ref_date is None:
        ref_date = datetime(2000,1,1)

    
    
    valid_zeta = ds.zeta.where(ds.zeta.values!=fill_value)
    valid_temp = ds.temp.where(ds.temp.values!=fill_value)
    valid_salt = ds.salt.where(ds.salt.values!=fill_value)
    valid_u = ds.u.where(ds.u.values!=fill_value)
    valid_v = ds.v.where(ds.v.values!=fill_value)
    
    if regrid_ver == 1:
        valid_mask = ds.mask.where(ds.mask.values!=fill_value)
        valid_depth = ds.depth.where(ds.depth.values!=fill_value)
        valid_h = ds.h.where(ds.h.values!=fill_value)
        data_vars={
            "lat_rho": (["eta_rho", "xi_rho"], ds.lat_rho.data, {
                "standard_name": "latitude",
                "long_name": "latitude",
                "units": "degrees_nort",
                "valid_min": np.min(ds.lat_rho.data),
                "valid_max": np.max(ds.lat_rho.data),
                "_FillValue": fill_value,
                "coordinates": "eta_rho xi_rho",
                }),
            "lon_rho": (["eta_rho", "xi_rho"], ds.lon_rho.data, {
                "standard_name": "longitude",
                "long_name": "longitude",
                "units": "degrees_east",
                "valid_min": np.min(ds.lon_rho.data),
                "valid_max": np.max(ds.lon_rho.data),
                "_FillValue": fill_value,
                "coordinates": "eta_rho xi_rho",
                }),
            "h": (["eta_rho", "xi_rho"], ds.h.data, {
                "standard_name": "sea_floor_depth",
                "long_name": "Depth of the sea floor",
                "units": "m",
                "valid_min": valid_h.min().item(),
                "valid_max": valid_h.max().item(),
                "_FillValue": fill_value,
                "coordinates": "eta_rho xi_rho",
                }),
            "mask": (["eta_rho", "xi_rho"], ds.mask.data, {
                "standard_name": "land_binary_mask",
                "long_name": "Land-sea mask (1=water, 0=land)",
                "units": "1",
                "flag_values":[0, 1],
                "flag_meanings": "land water",
                "valid_min": valid_mask.min().item(),
                "valid_max": valid_mask.max().item(),
                "_FillValue": fill_value,
                "coordinates": "eta_rho xi_rho",
                }),
            "zeta": (["time", "eta_rho", "xi_rho"], ds.zeta.data, {
                "standard_name": "sea_surface_elevation",
                "long_name": "Sea Surface Elevation",
                "units": "m",
                "valid_min": valid_zeta.min().item(),
                "valid_max": valid_zeta.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time eta_rho xi_rho",
                }),        
            "depth": (["time", "s_rho", "eta_rho", "xi_rho"], ds.depth.data, {
                "standard_name": "depth",
                "long_name": "Depth",
                "units": "m",
                "positive":"up",
                "valid_min": valid_depth.min().item(),
                "valid_max": valid_depth.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho eta_rho xi_rho",
                }),
            "temp": (["time", "s_rho", "eta_rho", "xi_rho"], ds.temp.data, {
                "standard_name": "sea_water_temperature",
                "long_name": "Sea Water Temperature",
                "units": "degC",
                "valid_min": valid_temp.min().item(),
                "valid_max": valid_temp.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho eta_rho xi_rho",
                }),
            "salt": (["time", "s_rho", "eta_rho", "xi_rho"], ds.salt.data, {
                "standard_name": "sea_water_salinity",
                "long_name": "Sea Water Salinity",
                "units": "1",
                "valid_min": valid_salt.min().item(),
                "valid_max": valid_salt.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho eta_rho xi_rho",
                }),
            "u": (["time", "s_rho", "eta_rho", "xi_rho"], ds.u.data, {
                "standard_name": "eastward_sea_water_velocity",
                "long_name": "Eastward Sea Water Velocity",
                "units": "m s-1",
                "valid_min": valid_u.min().item(),
                "valid_max": valid_u.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho eta_rho xi_rho",
                }),
            "v": (["time", "s_rho", "eta_rho", "xi_rho"], ds.v.data, {
                "standard_name": "northward_sea_water_velocity",
                "long_name": "Northward Sea Water Velocity",
                "units": "m s-1",
                "valid_min": valid_v.min().item(),
                "valid_max": valid_v.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho eta_rho xi_rho",
                }),
            }
    elif regrid_ver == 2:
        valid_h = ds.h.where(ds.h.values!=fill_value)
        
        data_vars={
            "lat_rho": (["eta_rho", "xi_rho"], ds.lat_rho.data, {
                "standard_name": "latitude",
                "long_name": "latitude",
                "units": "degrees_nort",
                "valid_min": np.min(ds.lat_rho.data),
                "valid_max": np.max(ds.lat_rho.data),
                "_FillValue": fill_value,
                "coordinates": "eta_rho xi_rho",
                }),
            "lon_rho": (["eta_rho", "xi_rho"], ds.lon_rho.data, {
                "standard_name": "longitude",
                "long_name": "longitude",
                "units": "degrees_east",
                "valid_min": np.min(ds.lon_rho.data),
                "valid_max": np.max(ds.lon_rho.data),
                "_FillValue": fill_value,
                "coordinates": "eta_rho xi_rho",
                }),
            "h": (["eta_rho", "xi_rho"], ds.h.data, {
                "standard_name": "sea_floor_depth",
                "long_name": "Depth of the sea floor",
                "units": "m",
                "valid_min": valid_h.min().item(),
                "valid_max": valid_h.max().item(),
                "_FillValue": fill_value,
                "coordinates": "eta_rho xi_rho",
                }),
            "zeta": (["time", "eta_rho", "xi_rho"], ds.zeta.data, {
                "standard_name": "sea_surface_elevation",
                "long_name": "Sea Surface Elevation",
                "units": "m",
                "valid_min": valid_zeta.min().item(),
                "valid_max": valid_zeta.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time eta_rho xi_rho",
                }),
            "temp": (["time", "depth", "eta_rho", "xi_rho"], ds.temp.data, {
                "standard_name": "sea_water_temperature",
                "long_name": "Sea Water Temperature",
                "units": "degC",
                "valid_min": valid_temp.min().item(),
                "valid_max": valid_temp.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho eta_rho xi_rho",
                }),
            "salt": (["time", "depth", "eta_rho", "xi_rho"], ds.salt.data, {
                "standard_name": "sea_water_salinity",
                "long_name": "Sea Water Salinity",
                "units": "1",
                "valid_min": valid_salt.min().item(),
                "valid_max": valid_salt.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho eta_rho xi_rho",
                }),
            "u": (["time", "depth", "eta_rho", "xi_rho"], ds.u.data, {
                "standard_name": "eastward_sea_water_velocity",
                "long_name": "Eastward Sea Water Velocity",
                "units": "m s-1",
                "valid_min": valid_u.min().item(),
                "valid_max": valid_u.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho eta_rho xi_rho",
                }),
            "v": (["time", "depth", "eta_rho", "xi_rho"], ds.v.data, {
                "standard_name": "northward_sea_water_velocity",
                "long_name": "Northward Sea Water Velocity",
                "units": "m s-1",
                "valid_min": valid_v.min().item(),
                "valid_max": valid_v.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho eta_rho xi_rho",
                }),
            }
    elif regrid_ver == 3:
        
        data_vars={
            "zeta": (["time", "latitude", "longitude"], ds.zeta.data, {
                "standard_name": "sea_surface_elevation",
                "long_name": "Sea Surface Elevation",
                "units": "m",
                "valid_min": valid_zeta.min().item(),
                "valid_max": valid_zeta.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time latitude longitude",
                }),
            "temp": (["time", "depth", "latitude", "longitude"], ds.temp.data, {
                "standard_name": "sea_water_temperature",
                "long_name": "Sea Water Temperature",
                "units": "degC",
                "valid_min": valid_temp.min().item(),
                "valid_max": valid_temp.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho latitude longitude",
                }),
            "salt": (["time", "depth", "latitude", "longitude"], ds.salt.data, {
                "standard_name": "sea_water_salinity",
                "long_name": "Sea Water Salinity",
                "units": "1",
                "valid_min": valid_salt.min().item(),
                "valid_max": valid_salt.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho latitude longitude",
                }),
            "u": (["time", "depth", "latitude", "longitude"], ds.u.data, {
                "standard_name": "eastward_sea_water_velocity",
                "long_name": "Eastward Sea Water Velocity",
                "units": "m s-1",
                "valid_min": valid_u.min().item(),
                "valid_max": valid_u.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho latitude longitude",
                }),
            "v": (["time", "depth", "latitude", "longitude"], ds.v.data, {
                "standard_name": "northward_sea_water_velocity",
                "long_name": "Northward Sea Water Velocity",
                "units": "m s-1",
                "valid_min": valid_v.min().item(),
                "valid_max": valid_v.max().item(),
                "_FillValue": fill_value,
                "coordinates": "time s_rho latitude longitude",
                }),
            }
    
    return data_vars



def write_cf_compliant_netCDF(file_out, ds_in, ref_date, info_dir, regrid_ver):

    info = pd.read_excel(info_path)

    # Use a fill value to ensure Compliance
    fill_value = -9999.0

    ds_in = ds_in.fillna(fill_value)

    ds = xr.Dataset(
        coords=write_coords(ds_in,regrid_ver),
        data_vars=write_vars(ds_in,regrid_ver),
        attrs=write_attrs(info)
    )

    ds.to_netcdf(file_out, format="NETCDF4")

    ds.close()

    print('')
    print(f'Created: {file_out}')
    print('')

def regrid1_cf_compliant(fname_in,info_path,out_dir=None,ref_date=None):
    
    if ref_date is None:
        ref_date = datetime(2000,1,1)
    
    if type(fname_in) == str:
        if fname_in.find('*') < 0: 
            fname_in = [fname_in]
        elif fname_in.find('*') > 0:
            fname_in = glob(fname_in)
        else:
            print('Error: unkown input format. Input variable fname_in needs to be str or list.')
            sys.exit()
    elif type(fname_in) == list:
        fname_in=fname_in
    else:
        print('Error: unkown input format. Input variable fname_in needs to be str or list.')
        sys.exit()

    
    for file in fname_in:
        print('')
        print('Opening: ', file)
        print('')
        
        print("Extracting the model output variables we need")
    
        # using the get_var() function so masking is included
        # returns datasets
        ds_temp = post.get_var(file, "temp", ref_date=ref_date)
        ds_salt = post.get_var(file, "salt", ref_date=ref_date)
        ds_u    = post.get_var(file, "u"   , ref_date=ref_date)
        ds_v    = post.get_var(file, "v"   , ref_date=ref_date)

        ds_out = xr.Dataset(
            {
                "lon_rho" : (["eta_rho", "xi_rho"], ds_temp.lon_rho.data),
                "lat_rho" : (["eta_rho", "xi_rho"], ds_temp.lat_rho.data),
                "h"       : (["eta_rho", "xi_rho"], ds_temp.h.data),
                "mask"    : (["eta_rho", "xi_rho"], ds_temp.mask.data),
                "zeta"    : (["time", "eta_rho", "xi_rho"], ds_temp.zeta.data),
                "depth"   : (["time", "s_rho", "eta_rho", "xi_rho"], ds_temp.depth.data),
                "temp"    : (["time", "s_rho", "eta_rho", "xi_rho"], ds_temp.temp.data),
                "salt"    : (["time", "s_rho", "eta_rho", "xi_rho"], ds_salt.salt.data),
                "u"       : (["time", "s_rho", "eta_rho", "xi_rho"], ds_u.u.data),
                "v"       : (["time", "s_rho", "eta_rho", "xi_rho"], ds_v.v.data)
                
                },
            
            coords={
                
                "time"    : ds_temp.time,
                "s_rho"   : ds_temp.s_rho,
                "eta_rho" : ds_temp.eta_rho,
                "xi_rho"  : ds_temp.xi_rho,
                
            })
        

        # Ensure the directory for the specified output exists
        if out_dir is None:
            out_dir = os.path.join(os.path.dirname(file), 'regrid_tier_1/',)
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        else:
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        
        fname_out = os.path.abspath(os.path.join(os.path.dirname(out_dir), os.path.basename(file[:-3]) + '_regrid1_CF_COMPLIANT.nc'))

        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass
        
        try:
            write_cf_compliant_netCDF(fname_out, ds_out, ref_date, info_path, regrid_ver=1)
        except Exception as e:
            print(f"Error: {e}")

def regrid2_cf_compliant(fname_in,info_dir,out_dir=None,depths=None,ref_date=None):
    if type(fname_in) == str:
        if fname_in.find('*') < 0:
            fname_in = [fname_in]
        elif fname_in.find('*') > 0:
            fname_in = glob(fname_in)
        else:
            print('Error: unkown input format. Input variable fname_in needs to be str or list.')
            sys.exit()
    elif type(fname_in) == list:
        fname_in=fname_in
    else:
        print('Error: unkown input format. Input variable fname_in needs to be str or list.')
        sys.exit()

    if depths is None:
        depths = [0,-5,-10,-20,-50,-75,-100,-200,-500,-1000]

    if ref_date is None:
        ref_date = datetime(2000,1,1)
    
    for file in fname_in:
        print('')
        print(f'Opening: {file}')
        print('')
    
        ds = xr.open_dataset(file)

        # Create an empty dataset with NaNs
        ds_out = xr.Dataset(
            
            {
                "lon_rho" : (["eta_rho", "xi_rho"], ds.lon_rho.data),
                "lat_rho" : (["eta_rho", "xi_rho"], ds.lat_rho.data),
                "h"       : (["eta_rho", "xi_rho"], ds.h.data),
                "mask"    : (["eta_rho", "xi_rho"], ds.lat_rho.data),
                "zeta"    : (["time", "eta_rho", "xi_rho"], ds.zeta.data),
                "temp"    : (["time", "depth", "eta_rho", "xi_rho"], np.full((np.size(ds.time), np.size(depths), np.size(ds.eta_rho),np.size(ds.xi_rho)), np.nan)),
                "salt"    : (["time", "depth", "eta_rho", "xi_rho"], np.full((np.size(ds.time), np.size(depths), np.size(ds.eta_rho),np.size(ds.xi_rho)), np.nan)),
                "u"       : (["time", "depth", "eta_rho", "xi_rho"], np.full((np.size(ds.time), np.size(depths), np.size(ds.eta_rho),np.size(ds.xi_rho)), np.nan)),
                "v"       : (["time", "depth", "eta_rho", "xi_rho"], np.full((np.size(ds.time), np.size(depths), np.size(ds.eta_rho),np.size(ds.xi_rho)), np.nan))
                
                },
            
            coords={
                
                "time"    : ds.time,
                "depths"  : np.array(depths),
                "eta_rho" : ds.eta_rho,
                "xi_rho"  : ds.xi_rho,
                
            })
        
        print("Doing the vertical interpolation to constant depth levels")
        print('')
        for i,d in enumerate(depths):
            print(f"Depth = {d} m")
            ds_out["temp"][:, i, :, :] = post.hlev_xarray(ds.temp, ds.depth, d)
            ds_out["salt"][:, i, :, :] = post.hlev_xarray(ds.salt, ds.depth, d)
            ds_out["u"][:, i, :, :]    = post.hlev_xarray(ds.u   , ds.depth, d)
            ds_out["v"][:, i, :, :]    = post.hlev_xarray(ds.v   , ds.depth, d)


        # Ensure the directory for the specified output exists
        if out_dir is None:
            out_dir = os.path.join(os.path.dirname(os.path.dirname(file)), 'regrid_tier_2/',)
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        else:
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        
        fname_out = os.path.abspath(os.path.join(os.path.dirname(out_dir), os.path.basename(file)[:18] + '_regrid2_CF_COMPLIANT.nc'))
            
        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass
        
        try:
            write_cf_compliant_netCDF(fname_out, ds_out, ref_date, info_path, regrid_ver=2)
        except Exception as e:
            print(f"Error: {e}")


def regrid3_cf_compliant(fname_in,out_dir=None,spacing=None,ref_date=None):
    if type(fname_in) == str:
        if fname_in.find('*') < 0:
            fname_in = [fname_in]
        elif fname_in.find('*') > 0:
            fname_in = glob(fname_in)
        else:
            print('Error: unkown input format. Input variable fname_in needs to be str or list.')
            sys.exit()
    elif type(fname_in) == list:
        fname_in=fname_in
    else:
        print('Error: unkown input format. Input variable fname_in needs to be str or list.')
        sys.exit()

    if spacing is None:
        spacing = 0.01

    if ref_date is None:
        ref_date = datetime(2000,1,1)

    for file in fname_in:
        print('')
        print(f'Opening: {file}')
        print('')

        ds = xr.open_dataset(file)
        Nt, Nz = np.shape(ds.temp[:,:,0,0].values)
        lon_rho_1d = np.ravel(ds.lon_rho.values)
        lat_rho_1d = np.ravel(ds.lat_rho.values)

        # input for griddata function later
        lonlat_input = np.array([lon_rho_1d, lat_rho_1d]).T

        # get the model boundary polygon
        lon_boundary = np.hstack(
            (ds.lon_rho.values[:, 0], ds.lon_rho.values[-1, 1:], ds.lon_rho.values[-1::-1, -1], ds.lon_rho.values[0, -2::-1])
        )

        lat_boundary = np.hstack(
            (ds.lat_rho.values[:, 0], ds.lat_rho.values[-1, 1:], ds.lat_rho.values[-1::-1, -1], ds.lat_rho.values[0, -2::-1])
        )

        # find the corners of the output regular grid (just big enough to cover the model domain)
        lon_min = np.floor(np.min(lon_boundary) / spacing) * spacing
        lon_max = np.ceil(np.max(lon_boundary) / spacing) * spacing
        lat_min = np.floor(np.min(lat_boundary) / spacing) * spacing
        lat_max = np.ceil(np.max(lat_boundary) / spacing) * spacing

        # generate the regular grid
        Nlon = int(np.rint((lon_max - lon_min) / spacing)) + 1
        Nlat = int(np.rint((lat_max - lat_min) / spacing)) + 1
        lon_out = np.linspace(lon_min, lon_max, num=Nlon, endpoint=True)
        lat_out = np.linspace(lat_min, lat_max, num=Nlat, endpoint=True)
        lon_out_grd, lat_out_grd = np.meshgrid(lon_out, lat_out)

        lon_out_grd = xr.DataArray(lon_out_grd, dims=("latitude", "longitude"))
        lat_out_grd = xr.DataArray(lat_out_grd, dims=("latitude", "longitude"))

        Ny, Nx = np.shape(lon_out_grd)

        # get a mask for the output grid which tells us which points are inside the CROCO model grid
        poly_boundary = mplPath.Path(np.array([lon_boundary, lat_boundary]).T)
        mask_out = np.zeros_like(lon_out_grd)
        for y in np.arange(Nlat):
            for x in np.arange(Nlon):
                if poly_boundary.contains_point((lon_out_grd[y, x], lat_out_grd[y, x])):
                    mask_out[y, x] = 1
                else:
                    mask_out[y, x] = 0

        # Here we do the horizontal interpolation
        # Even though the horizantal interpolation is made to work for an Xarray dataset, it does not
        # work in a loop when filling the xarray datasets because xarray is lazy!
        # When we do the interpolation for the variables in a loop, we have to write to numpy arrays,
        # otherwise xarray does not write the outputs properly and you end up with the incorrect values.
        # Therefore, we write to numpy arrays, then construct the xarray dataset!

        # make empty arrays
        zeta = np.zeros((Nt, Ny, Nx)) + np.nan
        temp = np.zeros((Nt, Nz, Ny, Nx)) + np.nan
        salt = np.zeros((Nt, Nz, Ny, Nx)) + np.nan
        u = np.zeros((Nt, Nz, Ny, Nx)) + np.nan
        v = np.zeros((Nt, Nz, Ny, Nx)) + np.nan

        # here we are just interested in calculating the elem and coef arrays.
        elem, coef, _ = horiz_interp_delaunay(ds.lon_rho, ds.lat_rho, ds.zeta[0], lon_out_grd,lat_out_grd)


        # do the main loop to do the horizontal interpolation
        for t in range(Nt):
            zeta[t,:,:]=horiz_interp_delaunay(ds.lon_rho,ds.lat_rho,ds.zeta[t],lon_out_grd,lat_out_grd,elem,coef).values
            for z in range(Nz):
                temp[t,z,:,:] = horiz_interp_delaunay(ds.lon_rho,ds.lat_rho,ds.temp[t,z,:,:],lon_out_grd,lat_out_grd,elem,coef).values

                salt[t,z,:,:] = horiz_interp_delaunay(ds.lon_rho[:],ds.lat_rho[:],ds.salt[t,z,:,:],lon_out_grd,lat_out_grd,elem,coef).values

                u[t,z,:,:] = horiz_interp_delaunay(ds.lon_rho,ds.lat_rho,ds.u[t,z,:,:],lon_out_grd,lat_out_grd,elem,coef).values

                v[t,z,:,:] = horiz_interp_delaunay(ds.lon_rho,ds.lat_rho,ds.v[t,z,:,:],lon_out_grd,lat_out_grd,elem,coef).values

        # Construvt the dataset which contains the regridded variables
        ds_out = xr.Dataset(
            coords={
                "time" : ds.time.values,
                "depths": ds.depth.values,
                "latitude" : lat_out,
                "longitude" : lon_out
            },
            data_vars={
                "zeta": (["time", "latitude", "longitude"], zeta),
                "temp": (["time", "depths", "latitude", "longitude"], temp),
                "salt": (["time", "depths", "latitude", "longitude"], salt),
                "u": (["time", "depths", "latitude", "longitude"], u),
                "v": (["time", "depths", "latitude", "longitude"], v)
            }
        )

        # Ensure the directory for the specified output exists
        if out_dir is None:
            out_dir = os.path.join(os.path.dirname(os.path.dirname(file)), 'regrid_tier_3/',)
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        else:
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)

        fname_out = os.path.abspath(os.path.join(os.path.dirname(out_dir), os.path.basename(file)[:18] + '_regrid3_CF_COMPLIANT.nc'))

        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass

        try:
            write_cf_compliant_netCDF(fname_out, ds_out, ref_date, info_path, regrid_ver=3)
        except Exception as e:
            print(f"Error: {e}")

        ds.close()
        ds_out.close()

if __name__ == "__main__":
    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg_Y2009M09.nc'
    info_path = '/home/g.rautenbach/Scripts/metadata/SAEON ODP metadata sheet - ISO19115 (data provider).xlsx'
    regrid1_cf_compliant(fname_in,info_path)

    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/regrid_tier_1/croco_avg_Y2009M09_regrid1_CF_COMPLIANT.nc'
    info_path = '/home/g.rautenbach/Scripts/metadata/SAEON ODP metadata sheet - ISO19115 (data provider).xlsx'
    regrid2_cf_compliant(fname_in,info_path)

    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/regrid_tier_2/croco_avg_Y2009M09_regrid2_CF_COMPLIANT.nc'
    ds_out = regrid3_cf_compliant(fname_in)
