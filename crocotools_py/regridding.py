import crocotools_py.postprocess as post
import numpy as np
import xarray as xr
import os,sys
from datetime import datetime
import pandas as pd
import matplotlib.path as mplPath
from dask import delayed
import dask.array as da
from scipy.interpolate import griddata
from glob import glob

def write_attrs(info=None,doi_link=None):
    if info is not None:
        attrs={
            "title"          : info['Description/Guidance'][3],
            "institution"    : info['Description/Guidance'][6],
            "source"         : "CROCO model",
            "history"        : "Created " + str(datetime.strftime(datetime.now(),"%Y-%m-%d %H:%M:%S")),
            "conventions"    : "CF-1.8",
            "author"         : info['Description/Guidance'][4],
            "contact person" : info['Description/Guidance'][5],
            "summery"        : info['Description/Guidance'][10],
            "keywords"       : info['Description/Guidance'][8],
        }

    else:
        attrs={
            "source"         : "CROCO model",
            "history"        : "Created " + str(datetime.strftime(datetime.now(),"%Y-%m-%d %H:%M:%S")),
            "conventions"    : "CF-1.8",
        }

    if doi is not None: 
        attrs["reference"] = f"doi: {doi}"
        
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

def write_netCDF(ds_in, file_out, regrid_ver, doi_link=None, ref_date=None, info_dir=None):
    
    if info_dir is not None:
        info_file = pd.read_excel(info_dir)
    else:
        info_file = None

    # Use a fill value to ensure Compliance
    fill_value = -9999.0

    ds_in = ds_in.fillna(fill_value)
    
    ds = xr.Dataset(
        coords=write_coords(ds_in,regrid_ver),
        data_vars=write_vars(ds_in,regrid_ver),
        attrs=write_attrs(info_file,doi_link)
    )
    
    ds.to_netcdf(file_out, format="NETCDF4")

    ds.close()

    print('')
    print(f'Created: {file_out}')
    print('')


def regrid_tier1(fname_in,fname_out,info_file=None,doi_link=None,ref_date=None):
    '''
    tier 1 regridding of a raw CROCO output file:
        -> regrids u/v to the density (rho) grid so all parameters are on the same horizontal grid
        -> rotates u/v to be east/north components instead of grid-aligned components
        -> adds a 'depth' variable providing the depths of each sigma level at each time-step

    Parameters
    ----------
    fname_in  : path to input CROCO file(s) (required = True)
    fname_out : path to output tier 1 netcdf file (required = True)
    info_file : path to CF-Compliance spreadsheet (required = False)
    doi_link  : doi link in string (required = False)
    ref_date  : reference datetime used in croco runs (must be a datetime.datetime object, required = False, standard = 2000,1,1)
    
    '''

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
        
        fname_out = os.path.abspath(os.path.join(os.path.dirname(fname_out), os.path.basename(file[:-3]) + '_t1.nc'))
        
        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass

        try:
            write_netCDF(ds_out, fname_out, regrid_ver=1, doi_link=doi_link, ref_date=ref_date, info_dir=info_file)
        except Exception as e:
            print(f"Error: {e}")

def regrid_tier2(fname_in,fname_out,info_file=None,doi_link=None,ref_date=None, depths=None):
    '''
    tier 2 regridding of a CROCO output:
      -> as per tier1 regridding but we regrid vertically to constant z levels
      -> output variables are the same as tier 1, only 'depth' is now a dimension with the user specified values

    Parameters
    ----------
    fname_in  : path to input tier 2 netcdf file (required = True).
    fname_out : path to output tier 2 netcdf file (required = True).
    info_file : path to CF-Compliance spreadsheet (required = False).
    doi_link  : doi link in string (required = False).
    ref_date  : reference datetime used in croco runs (must be a datetime.datetime object, required = False, standard = 2000,1,1).
    depths    : list of depths to extract (in metres, negative down, required = False).
                If not specified depth = [0,-5,-10,-20,-50,-75,-100,-200,-500,-1000].
                A value of 0 denotes the surface and a value of -99999 denotes the bottom layer.

    '''


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


        fname_out = os.path.abspath(os.path.join(os.path.dirname(fname_out), os.path.basename(file)[:18] + '_t2.nc'))
            
        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass
        
        try:
            write_netCDF(ds_out, fname_out, regrid_ver=2, doi_link=doi_link, ref_date=ref_date, info_dir=info_file)
        except Exception as e:
            print(f"Error: {e}")

def regrid_tier3(fname_in,fname_out,info_file=None,doi_link=None,ref_date=None,spacing=None):
    '''
    tier 3 regridding of a CROCO output:
      -> takes the output of regrid-tier2 as input and
      -> regrids the horizontal grid to be regular with a specified grid spacing
      -> output variables are the same as tier 1 and 2, 
         only horizontal grid is now rectilinear with hz dimensions of longitude,latitude
         i.e. horizontal grid is no longer curvilinear
         the extents of the rectilinear grid are automatically determined using the curvilinear grid extents
         
    Parameters
    ----------
    fname_in  : path to input tier 2 netcdf file (required = True)
    fname_out : path to output tier 3 netcdf file (required = True)
    info_file : path to CF-Compliance spreadsheet (required = False)
    doi_link  : doi link in string (required = False)
    ref_date  : reference datetime used in croco runs (must be a datetime.datetime object, required = False, standard = 2000,1,1)
    spacing   : constant horizontal grid spacing (in degrees) to be used for the horizontal interpolation of the output
    

    CAREFUL! tier3 output is useful for website visualisation (it's intended use), 
    but don't use it for reasearch/analysis as it's interpolated, so can be at a 
    totally different resolution to the native CROCO grid 

    '''

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

        @delayed
        def compute_2d_chunk(t, variable, method="nearest"):
            return (
                griddata(
                    lonlat_input,
                    np.ravel(variable[t, ::]),
                    (lon_out_grd, lat_out_grd),
                    method,
                )
                * mask_out / mask_out
            )
    
        @delayed
        def compute_3d_chunk(t, variable, n, method="nearest"):
            return (
                griddata(
                    lonlat_input,
                    np.ravel(variable[t, n, ::]),
                    (lon_out_grd, lat_out_grd),
                    method,
                )
                * mask_out / mask_out
            )



        # Separate lists for each time step
        zeta_out = []
        temp_out_time = []
        salt_out_time = []
        u_out_time = []
        v_out_time = []

        print('Interpolating the model output onto the regular horizontal output grid...')

        for t in np.arange(Nt):
            # Lists for each depth level
            temp_out_depth = []
            salt_out_depth = []
            u_out_depth = []
            v_out_depth = []

            zeta_out.append(
                da.from_delayed(
                    compute_2d_chunk(t, ds.zeta.values),
                    shape=(Nlat, Nlon),
                    dtype=float,
                )
            )

            for n in np.arange(Nz):
                temp_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, ds.temp.values, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
                salt_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, ds.salt.values, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
                u_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, ds.u.values, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
                v_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, ds.v.values, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
    
    
            # Stack the depth dimension and append to the time list
            temp_out_time.append(da.stack(temp_out_depth, axis=0))
            salt_out_time.append(da.stack(salt_out_depth, axis=0))
            u_out_time.append(da.stack(u_out_depth, axis=0))
            v_out_time.append(da.stack(v_out_depth, axis=0))
    
        # Stack the time dimension
        zeta_out = da.stack(zeta_out, axis=0)
        temp_out = da.stack(temp_out_time, axis=0)
        salt_out = da.stack(salt_out_time, axis=0)
        u_out = da.stack(u_out_time, axis=0)
        v_out = da.stack(v_out_time, axis=0)

        # Construct the dataset which contains the regridded variables       
        ds_out = xr.Dataset(
            coords={
                "time" : ds.time.values,
                "depths": ds.depth.values,
                "latitude" : lat_out,
                "longitude" : lon_out
            },
            data_vars={
                "zeta": (["time", "latitude", "longitude"], zeta_out.compute()),
                "temp": (["time", "depths", "latitude", "longitude"], temp_out.compute()),
                "salt": (["time", "depths", "latitude", "longitude"], salt_out.compute()),
                "u": (["time", "depths", "latitude", "longitude"], u_out.compute()),
                "v": (["time", "depths", "latitude", "longitude"], v_out.compute())
            }
        )

        fname_out = os.path.abspath(os.path.join(os.path.dirname(fname_out), os.path.basename(file)[:18] + '_t3.nc'))
        
        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass
            
        try:
            write_netCDF(ds_out, fname_out, regrid_ver=3, doi_link=doi_link, ref_date=ref_date, info_dir=info_file)
        except Exception as e:
            print(f"Error: {e}")

        ds.close()
        ds_out.close()


if __name__ == "__main__":
    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg_Y2009M09.nc'
    fname_out = '/home/g.rautenbach/Data/models/sa_southeast/'
    info_dir = '/home/g.rautenbach/Scripts/metadata/SAEON ODP metadata sheet - ISO19115 (data provider).xlsx'
    doi = '10.15493/SOMISANA.26032025'
    regrid_tier1(fname_in,fname_out,info_dir,doi)

    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg_Y2009M09_t1.nc'
    regrid_tier2(fname_in,fname_out,info_dir,doi)

    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg_Y2009M09_t2.nc'
    regrid_tier3(fname_in,fname_out,info_dir,doi)
