"""
Some scripts for regridding CROCO output with the hope that the files 
are easier to use for non-CROCO folk.
These functions are used in the regridding of SOMISANA's operational CROCO output
I don't think the output netcdf files are totally CF compliant - 
something to come back to when we start archiving on MIMS!
"""

import crocotools_py.postprocess as post
import numpy as np
import xarray as xr
import os
from datetime import datetime
import subprocess
from dask import delayed
import dask
import dask.array as da
from scipy import interpolate
import matplotlib.path as mplPath

def regrid_tier1(fname_in, fname_out, ref_date):
    '''
    tier 1 regridding of a raw CROCO output file:
        -> regrids u/v to the density (rho) grid so all parameters are on the same horizontal grid
        -> rotates u/v to be east/north components instead of grid-aligned components
        -> adds a 'depth' variable providing the depths of each sigma level at each time-step

    Parameters
    ----------
    fname_in : path to input CROCO file(s)
    fname_out : path to output tier 1 netcdf file
    ref_date : reference datetime used in croco runs
               (must be a datetime.datetime object)
    
    '''

    # Ensure the directory for the specified output exists
    os.makedirs(os.path.dirname(fname_out), exist_ok=True)

    print("Extracting the model output variables we need")
    
    ds_temp = post.get_var(fname_in, "temp", ref_date=ref_date)
    ds_salt = post.get_var(fname_in, "salt", ref_date=ref_date)
    ds_uv = post.get_uv(fname_in, ref_date=ref_date)
    
    ds_all = xr.merge([ds_temp,ds_salt,ds_uv])
    
    encoding = {
        "zeta": {"dtype": "float32"},
        "temp": {"dtype": "float32"},
        "salt": {"dtype": "float32"},
        "u": {"dtype": "float32"},
        "v": {"dtype": "float32"},
        "depth": {"dtype": "float32"},
        "h": {"dtype": "float32"},
        "lon_rho": {"dtype": "float32"},
        "lat_rho": {"dtype": "float32"},
        "s_rho": {"dtype": "float32"},
        "time": {"dtype": "i4"},
    }
    
    print('writing the netcdf file')
    ds_all.to_netcdf(fname_out, encoding=encoding, mode="w")
    
    subprocess.call(["chmod", "-R", "775", fname_out])

    print("Done!")
    
def regrid_tier2(fname_in, fname_out, ref_date,
                 depths=[0,-5,-10,-20,-50,-75,-100,-200,-500,-1000]):
    '''
    tier 2 regridding of a CROCO output:
      -> as per tier1 regridding but we regrid vertically to constant z levels
      -> output variables are the same as tier 1, only 'depth' is now a dimension with the user specified values

    Parameters
    ----------
    fname_in : path to input CROCO file(s)
    fname_out : path to output tier 2 netcdf file
    depths : list of depths to extract (in metres, negative down). 
            A value of 0 denotes the surface and a value of -99999 denotes the bottom layer)

    '''
    
    print("Extracting the model output variables we need")
    
    ds_temp = post.get_var(fname_in, "temp", ref_date=ref_date, level=depths)
    ds_salt = post.get_var(fname_in, "salt", ref_date=ref_date, level=depths)
    ds_uv = post.get_uv(fname_in, ref_date=ref_date, level=depths)
    
    ds_all = xr.merge([ds_temp,ds_salt,ds_uv])
    
    encoding = {
        "zeta": {"dtype": "float32"},
        "temp": {"dtype": "float32"},
        "salt": {"dtype": "float32"},
        "u": {"dtype": "float32"},
        "v": {"dtype": "float32"},
        "depth": {"dtype": "float32"},
        "h": {"dtype": "float32"},
        "lon_rho": {"dtype": "float32"},
        "lat_rho": {"dtype": "float32"},
        "time": {"dtype": "i4"},
        
    }

    print("Writing NetCDF file")
    ds_all.to_netcdf(fname_out, encoding=encoding, mode="w")

    subprocess.call(["chmod", "-R", "775", fname_out])
    
    ds_all.close()
    
    print("Done!")

def regrid_tier3(fname_in, fname_out, spacing='0.01'):
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
    fname_in : path to input tier 2 netcdf file
    fname_out : path to output tier 3 netcdf file
    spacing : constant horizontal grid spacing (in degrees) to be used for the 
              horizontal interpolation of the output
    

    CAREFUL! tier3 output is useful for website visualisation (it's intended use), 
    but don't use it for reasearch/analysis as it's interpolated, so can be at a 
    totally different resolution to the native CROCO grid 

    '''
        
    spacing = float(spacing)
    
    print("Extracting the tier 2 re-gridded model output")
    with xr.open_dataset(fname_in) as ds:
        zeta_in = ds.zeta.values
        temp_in = ds.temp.values
        salt_in = ds.salt.values
        u_in = ds.u.values
        v_in = ds.v.values
        depth = ds.depth.values
        lon_rho = ds.lon_rho.values
        lat_rho = ds.lat_rho.values
        Nt, Nz, Ny, Nx = np.shape(temp_in)
        lon_rho_1d = np.ravel(lon_rho)
        lat_rho_1d = np.ravel(lat_rho)

        # input for griddata function later
        lonlat_input = np.array([lon_rho_1d, lat_rho_1d]).T

        print("Generating the regular horizontal output grid")
        # get the model boundary polygon
        lon_boundary = np.hstack(
            (lon_rho[0:, 0], lon_rho[-1, 1:-1], lon_rho[-1::-1, -1], lon_rho[0, -2::-1])
        )
        lat_boundary = np.hstack(
            (lat_rho[0:, 0], lat_rho[-1, 1:-1], lat_rho[-1::-1, -1], lat_rho[0, -2::-1])
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

        # get a mask for the output grid which tells us which points are inside the CROCO model grid
        poly_boundary = mplPath.Path(np.array([lon_boundary, lat_boundary]).T)
        mask_out = np.zeros_like(lon_out_grd)
        for y in np.arange(Nlat):
            for x in np.arange(Nlon):
                if poly_boundary.contains_point((lon_out_grd[y, x], lat_out_grd[y, x])):
                    mask_out[y, x] = 1
                else:
                    mask_out[y, x] = np.nan

        print("Interpolating the model output onto the regular horizontal output grid")

        @delayed
        def compute_2d_chunk(t, variable, method="nearest"):
            return (
                interpolate.griddata(
                    lonlat_input,
                    np.ravel(variable[t, ::]),
                    (lon_out_grd, lat_out_grd),
                    method,
                )
                * mask_out
            )

        @delayed
        def compute_3d_chunk(t, variable, n, method="nearest"):
            return (
                interpolate.griddata(
                    lonlat_input,
                    np.ravel(variable[t, n, ::]),
                    (lon_out_grd, lat_out_grd),
                    method,
                )
                * mask_out
            )

        # Separate lists for each time step
        zeta_out = []
        temp_out_time = []
        salt_out_time = []
        u_out_time = []
        v_out_time = []

        for t in np.arange(Nt):
            # Lists for each depth level
            temp_out_depth = []
            salt_out_depth = []
            u_out_depth = []
            v_out_depth = []

            zeta_out.append(
                da.from_delayed(
                    compute_2d_chunk(t, zeta_in),
                    shape=(Nlat, Nlon),
                    dtype=float,
                )
            )

            for n in np.arange(Nz):
                print(
                    f"Timestep {str(t+1).zfill(3)}/{Nt}. Depth {str(np.round(depth[n]).astype(int)).zfill(5)}m (lvl {str(n+1).zfill(2)}/{Nz})."
                )
                temp_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, temp_in, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
                salt_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, salt_in, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
                u_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, u_in, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
                v_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, v_in, n),
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

        # Create new xarray dataset with selected variables
        print("Generating dataset")
        data_out = xr.Dataset(
            attrs={
                "description": "tier 3 regridded CROCO output: as per tier2, but data interpolated to a constant horizontal grid.",
            },
            data_vars={
                "zeta": xr.Variable(
                    ["time", "latitude", "longitude"],
                    zeta_out,
                    {
                        "long_name": "averaged free-surface",
                        "units": "meter",
                        "standard_name": "sea_surface_height",
                    },
                ),
                "temp": xr.Variable(
                    ["time", "depth", "latitude", "longitude"],
                    temp_out,
                    {
                        "long_name": "averaged potential temperature",
                        "units": "Celsius",
                        "standard_name": "sea_water_potential_temperature",
                    },
                ),
                "salt": xr.Variable(
                    ["time", "depth", "latitude", "longitude"],
                    salt_out,
                    {
                        "long_name": "averaged salinity",
                        "units": "PSU",
                        "standard_name": "sea_water_salinity",
                    },
                ),
                "u": xr.Variable(
                    ["time", "depth", "latitude", "longitude"],
                    u_out,
                    {
                        "long_name": "Eastward velocity",
                        "units": "meters per second",
                        "standard_name": "eastward_sea_water_velocity",
                    },
                ),
                "v": xr.Variable(
                    ["time", "depth", "latitude", "longitude"],
                    v_out,
                    {
                        "long_name": "Northward velocity",
                        "units": "meters per second",
                        "standard_name": "northward_sea_water_velocity",
                    },
                ),
            },
            coords={
                "longitude": xr.Variable(
                    ["longitude"],
                    lon_out,
                    {
                        "long_name": "Longitude",
                        "units": "degrees_east",
                        "standard_name": "longitude",
                    },
                ),
                "latitude": xr.Variable(
                    ["latitude"],
                    lat_out,
                    {
                        "long_name": "Latitude",
                        "units": "degrees_west",
                        "standard_name": "latitude",
                    },
                ),
                "time": xr.Variable(
                    ["time"],
                    ds.time.values,
                    {"description": "time"},
                ),
                "depth": xr.Variable(
                    ["depth"],
                    ds.depth.values,
                    {
                        "long_name": "water depth from free surface",
                        "units": "meter",
                        "positive": "up",
                        "standard_name": "depth",
                    },
                ),
            },
        )

        # Explicitly set chunk sizes of some dimensions
        chunksizes = {
            "time": 24,
            "depth": 1,
        }

        # For data_vars, set chunk sizes for each dimension
        # This is either the override specified in "chunksizes"
        # or the length of the dimension
        default_chunksizes = {dim: len(data_out[dim]) for dim in data_out.dims}

        encoding = {
            var: {
                "dtype": "float32", 
                "chunksizes": [chunksizes.get(dim, default_chunksizes[dim]) for dim in data_out[var].dims]
            }
            for var in data_out.data_vars
        }

        # Adjust for non-chunked variables
        encoding["time"] = {"dtype": "i4"}
        encoding['latitude'] = {"dtype": "float32"}
        encoding['longitude'] = {"dtype": "float32"}
        encoding['depth'] = {"dtype": "float32"}

        print("Generating NetCDF data")
        write_op = data_out.to_netcdf(
            fname_out,
            encoding=encoding,
            mode="w",
            compute=False,
        )

        print("Writing NetCDF file")
        dask.compute(write_op)

        subprocess.call(["chmod", "-R", "775", fname_out])
        print("Done!")

if __name__ == "__main__":
    
    # Example of how you'd use these functions
    ref_date = datetime(2000, 1, 1, 0, 0, 0)
    fname = '/home/gfearon/tmp/20240122_hourly_avg.nc'
    fname_t1 = '/home/gfearon/tmp/20240122_hourly_avg_tier1.nc'
    fname_t2 = '/home/gfearon/tmp/20240122_hourly_avg_tier2.nc'
    fname_t3 = '/home/gfearon/tmp/20240122_hourly_avg_tier3.nc'
    # regrid_tier1(fname, fname_t1, ref_date)
    regrid_tier2(fname_t1, fname_t2, depths=[0,-10,-100,-1000,-99999])
    # regrid_tier3(fname_t2, fname_t3, spacing='0.05')
