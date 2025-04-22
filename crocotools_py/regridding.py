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

def write_attrs(info=None, doi_link=None):
    attrs = {
        "title": "Model Outputs Regridded Tier 1",
        "source": "CROCO model",
        "history": "Created " + datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S"),
        "conventions": "CF-1.8",
    }

    if info is not None: 
        df = pd.read_excel(info)
        
        desc = df["Description/Guidance"].tolist()
        
        if isinstance(desc, list):
            attrs.update({
                "title": desc[3],
                "institution": desc[6],
                "author": desc[4],
                "contact person": desc[5],
                "summery": desc[10],
                "keywords": desc[8],
            })
            
    if doi_link is not None:
        attrs.update({
            "reference" :f"{doi_link}"
        })

    return attrs

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

    os.makedirs(os.path.dirname(fname_out), exist_ok=True)

    if ref_date is None:
        ref_date = datetime(2000,1,1,0,0)
    
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
    
        ds_temp = post.get_var(file, "temp", ref_date=ref_date)
        ds_salt = post.get_var(file, "salt", ref_date=ref_date)
        ds_uv   = post.get_uv(file,ref_date=ref_date)
        
        ds_out = xr.Dataset(
            {
                "lon_rho" : ds_temp.lon_rho,
                "lat_rho" : ds_temp.lat_rho,
                "h"       : ds_temp.h,
                "mask"    : ds_temp.mask,
                "zeta"    : ds_temp.zeta,
                "depth"   : ds_temp.depth,
                "temp"    : ds_temp.temp,
                "salt"    : ds_salt.salt,
                "u"       : ds_uv.u,
                "v"       : ds_uv.v
                
                },
            
            coords={
                   
                "time"    : ds_temp.time,
                "s_rho"   : ds_temp.s_rho,
                "eta_rho" : ds_temp.eta_rho,
                "xi_rho"  : ds_temp.xi_rho,
                
            }, 

            attrs = write_attrs(info=info_file, doi_link=doi_link)

            )
        
        ds_out["time"].attrs = {"long_name": "Time","standard_name":"time"}

        fname_out = os.path.abspath(os.path.join(os.path.dirname(fname_out), os.path.basename(file)[:-3] + '_t1.nc'))
        
        try:
            ds_out.to_netcdf(fname_out,format="NETCDF4",encoding={"time":{
                                                                            "units": f"hours since {ref_date.strftime('%Y-%m-%d %H:%M:%S')}",
                                                                            "calendar": "standard",
                                                                            "dtype": "float32"
                                                                            }
                                                                  }
                             )
            ds_out.close()
            print('')
            print(f'Created: {fname_out}')
            print('')
        except Exception as e:
            print(f"Error: {e}")

def regrid_tier2(fname_in,fname_out,info_file=None, doi_link=None,ref_date=None,depths=None):
    '''
    tier 2 regridding of a CROCO output:
      -> as per tier1 regridding but we regrid vertically to constant z levels
      -> output variables are the same as tier 1, only 'depth' is now a dimension with the user specified values

    Parameters
    ----------
    fname_in  : path to input tier 2 netcdf file (required = True).
    fname_out : path to output tier 2 netcdf file (required = True).
    ref_date  : reference datetime used in croco runs (must be a datetime.datetime object, required = False, standard = 2000,1,1).
    depths    : list of depths to extract (in metres, negative down, required = False).
                If not specified depth = [0,-5,-10,-20,-50,-75,-100,-200,-500,-1000].
                A value of 0 denotes the surface and a value of -99999 denotes the bottom layer.

    '''

    os.makedirs(os.path.dirname(fname_out), exist_ok=True)

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
                "lon_rho" : ds.lon_rho,
                "lat_rho" : ds.lat_rho,
                "h"       : ds.h,
                "mask"    : ds.mask,
                "zeta"    : ds.zeta,
                "temp"    : (["time", "depth", "eta_rho", "xi_rho"], np.full((np.size(ds.time), np.size(depths), np.size(ds.eta_rho),np.size(ds.xi_rho)), np.nan)),
                "salt"    : (["time", "depth", "eta_rho", "xi_rho"], np.full((np.size(ds.time), np.size(depths), np.size(ds.eta_rho),np.size(ds.xi_rho)), np.nan)),
                "u"       : (["time", "depth", "eta_rho", "xi_rho"], np.full((np.size(ds.time), np.size(depths), np.size(ds.eta_rho),np.size(ds.xi_rho)), np.nan)),
                "v"       : (["time", "depth", "eta_rho", "xi_rho"], np.full((np.size(ds.time), np.size(depths), np.size(ds.eta_rho),np.size(ds.xi_rho)), np.nan))
                
                },
            
            coords={
                
                "time"    : ds.time,
                "depth"  : np.array(depths),
                "eta_rho" : ds.eta_rho,
                "xi_rho"  : ds.xi_rho,
                
            },

            attrs=write_attrs(info=info_file, doi_link=doi_link)

            )
        
        ds_out["time"].attrs = ds.time.attrs


        print("Doing the vertical interpolation to constant depth levels")
        print('')
        for i,d in enumerate(depths):
            print(f"Depth = {d} m")
            ds_out["temp"][:, i, :, :] = post.hlev_xarray(ds.temp, ds.depth, d)
            ds_out["salt"][:, i, :, :] = post.hlev_xarray(ds.salt, ds.depth, d)
            ds_out["u"][:, i, :, :]    = post.hlev_xarray(ds.u   , ds.depth, d)
            ds_out["v"][:, i, :, :]    = post.hlev_xarray(ds.v   , ds.depth, d)
        
        ds_out["temp"].attrs = ds.temp.attrs
        ds_out["salt"].attrs = ds.salt.attrs
        ds_out["u"].attrs = ds.u.attrs
        ds_out["v"].attrs = ds.v.attrs
        ds_out["depth"].attrs = {"long_name": "Fixed depth levels","units": "m","standard_name": "depth"}

        fname_out = os.path.abspath(os.path.join(os.path.dirname(fname_out), os.path.basename(file)[:-6] + '_t2.nc'))

        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass

        try:
            ds_out.to_netcdf(fname_out, format="NETCDF4",encoding={"time":{
                                                                            "units": f"hours since {ref_date.strftime('%Y-%m-%d %H:%M:%S')}",
                                                                            "calendar": "standard",
                                                                            "dtype": "float32"
                                                                            }
                                                                  }
                             )
            ds_out.close()
            print('')
            print(f'Created: {fname_out}')
            print('')
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
                "time" : ds.time,
                "depth": ds.depth,
                "latitude" : lat_out,
                "longitude" : lon_out            
                },
            data_vars={
                "zeta": (["time", "latitude", "longitude"], zeta_out.compute()),
                "temp": (["time", "depths", "latitude", "longitude"], temp_out.compute()),
                "salt": (["time", "depths", "latitude", "longitude"], salt_out.compute()),
                "u": (["time", "depths", "latitude", "longitude"], u_out.compute()),
                "v": (["time", "depths", "latitude", "longitude"], v_out.compute())
                },
            attrs=write_attrs(
                info=info_file, doi_link=doi_link
                )
            )
        
        # We just reassign the old attributes because they are CF-Compliant
        ds_out["longitude"].attrs = ds.lon_rho.attrs
        ds_out["latitude"].attrs = ds.lat_rho.attrs
        ds_out["zeta"].attrs = ds.zeta.attrs
        ds_out["temp"].attrs = ds.temp.attrs
        ds_out["salt"].attrs = ds.salt.attrs
        ds_out["u"].attrs = ds.u.attrs
        ds_out["v"].attrs = ds.v.attrs
        
        fname_out = os.path.abspath(os.path.join(os.path.dirname(fname_out), os.path.basename(file)[:-6] + '_t3.nc'))
 
        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass

        try:
            ds_out.to_netcdf(fname_out, format="NETCDF4",encoding={"time":{
                                                                "units": f"hours since {ref_date.strftime('%Y-%m-%d %H:%M:%S')}",
                                                                "calendar": "standard",
                                                                "dtype": "float32"
                                                                }
                                                      }
                             )
            ds_out.close()
            print('')
            print(f'Created: {fname_out}')
            print('')
        except Exception as e:
            print(f"Error: {e}")

if __name__ == "__main__":
    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg.nc'
    fname_out = '/home/g.rautenbach/Data/models/sa_southeast/'
    info_dir = '/home/g.rautenbach/Scripts/metadata/SAEON ODP metadata sheet - ISO19115 (data provider).xlsx'
    doi = '10.15493/SOMISANA.26032025'
    ref_date=datetime(2000,1,1)

    #regrid_tier1(fname_in,fname_out,info_dir,doi,ref_date=ref_date)
    #fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg_t1.nc'
    #regrid_tier2(fname_in,fname_out,info_dir,doi,ref_date=ref_date,depths=[0,10,20])
    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg_t2.nc'
    regrid_tier3(fname_in,fname_out,info_dir,doi,ref_date=ref_date,spacing=0.02)
