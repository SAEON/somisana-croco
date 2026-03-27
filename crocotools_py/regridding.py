import crocotools_py.postprocess as post
import numpy as np
import xarray as xr
import os,sys
from datetime import datetime
#import pandas as pd
import matplotlib.path as mplPath
from dask import delayed
import dask.array as da
from scipy.interpolate import griddata
from glob import glob
import subprocess
import dask
#import dask.array as da
import re

# Known vector pairs in CROCO output (u-component, v-component)
# When both components are in varList, get_uv() is used to extract and rotate them
# When only one component is present, get_var() extracts it grid-aligned (no rotation)
VECTOR_PAIRS = [
    ('u', 'v'),
    ('ubar', 'vbar'),
    ('sustr', 'svstr'),
    ('bustr', 'bvstr'),
]

def _partition_vars(varList):
    """
    Split varList into scalar variables and vector pairs.
    Returns (scalars, pairs) where:
      - scalars: list of variable names to extract with get_var()
      - pairs: list of (var_u, var_v) tuples to extract with get_uv()
    """
    remaining = set(varList)
    pairs = []
    for var_u, var_v in VECTOR_PAIRS:
        if var_u in remaining and var_v in remaining:
            pairs.append((var_u, var_v))
            remaining.discard(var_u)
            remaining.discard(var_v)
    # preserve original ordering for scalars
    scalars = [v for v in varList if v in remaining]
    return scalars, pairs

def regrid_tier1(fname_in, dir_out, grdname=None, Yorig=2000, doi_link=None,
                 varList=['temp','salt','u','v']):
    '''
    tier 1 regridding of a raw CROCO output file(s):
        -> regrids u/v to the density (rho) grid so all parameters are on the same horizontal grid
        -> rotates u/v to be east/north components instead of grid-aligned components
        -> adds a 'depth' variable providing the depths of each sigma level at each time-step

    Parameters
    ----------
    fname_in  : path to input CROCO file(s). Can include wildcards *. (required = True)
    dir_out   : path to output directory (required = True)
    grdname   : optional name of your croco grid file (only needed if the grid info is not in fname)
    Yorig     : Origin year used in setting up CROCO time i.e. seconds since Yorig-01-01
    doi_link  : doi link in string (required = False)
    varList   : list of variable names to extract (default=['temp','salt','u','v']).
                If both components of a vector pair (e.g. u/v, ubar/vbar, sustr/svstr, bustr/bvstr)
                are present, they are extracted together via get_uv() and rotated to east/north components.
                If only one component is present, it is extracted via get_var() and remains grid-aligned.
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

    
    for file in fname_in:
        print('')
        print('Opening: ', file)
        print("Extracting the model output variables we need")
    
        scalars, pairs = _partition_vars(varList)
        datasets = []
        for var in scalars:
            datasets.append(post.get_var(file, var, grdname=grdname, Yorig=Yorig))
        for var_u, var_v in pairs:
            datasets.append(post.get_uv(file, grdname=grdname, Yorig=Yorig, var_u=var_u, var_v=var_v))
        ds_all = xr.merge(datasets, compat='override')
        
        ds_all.attrs["title"] = "Regridded CROCO output created by the regrid_tier1 function"
        ds_all.attrs["source"] = file
        ds_all.attrs["history"] = "Created on " + datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S")
        ds_all.attrs["conventions"] = "CF-1.8"
        ds_all.attrs['references'] = 'Project: Sustainable Ocean Modelling Initiative: a South AfricaN Approach (SOMISANA; https://somisana.ac.za/); Tools: Regridding Code (https://github.com/SAEON/somisana-croco)'
        if doi_link is not None: ds_all.attrs.update({"doi" :f"https://doi.org/{doi_link}"})

        # build encoding dynamically based on variables present in ds_all
        encoding = {}
        for var in ds_all.data_vars:
            encoding[var] = {"dtype": "float32"}
        for var in ['h', 'mask', 'lon_rho', 'lat_rho']:
            if var in ds_all:
                encoding[var] = {"dtype": "float32"}
        if 'depth' in ds_all.variables:  # allow for both cases where depth is or isn't in the dataset (e.g. if a surface file is used)
            encoding['depth'] = {"dtype": "float32"}
        encoding["time"] = {"units": f"seconds since {Yorig}-01-01 00:00:00",
                            "calendar": "standard",
                            "dtype": "i4"}
      
        # robust way of getting the file extension, including CROCO child domains e.g. *nc.2, *.nc.2 etc
        basename = os.path.basename(file)
        match = re.search(r"(\.nc\.\d+)$|(\.nc)$", basename)
        basename_no_extension = basename[:match.start()]
        extension = match.group(0)
        
        fname_out = os.path.abspath(os.path.join(dir_out, basename_no_extension + '_t1' + extension))

        # If the file already exists, we remove it to avoid permission errors.
        if os.path.exists(fname_out): os.remove(fname_out)        
        ds_all.to_netcdf(fname_out, encoding=encoding, mode="w")
        
        subprocess.call(["chmod", "-R", "775", fname_out])
        
        ds_all.close()

        print(f'Created: {fname_out}')

def regrid_tier2(fname_in, dir_out, grdname=None, Yorig=2000, doi_link=None,
                 depths=[0,-5,-10,-20,-50,-75,-100,-200,-500,-1000],
                 varList=['temp','salt','u','v']):
    '''
    tier 2 regridding of a CROCO output:
      -> as per tier1 regridding but we regrid vertically to constant z levels
      -> output variables are the same as tier 1, only 'depth' is now a dimension with the user specified values

    Parameters
    ----------
    fname_in  : path to input tier 2 netcdf file (required = True).
    dir_out   : path to output directory (required = True).
    grdname   : optional name of your croco grid file (only needed if the grid info is not in fname)
    Yorig     : Origin year used in setting up CROCO time i.e. seconds since Yorig-01-01
    depths    : list of depths to extract (in metres, negative down, required = False).
                If not specified depth = [0,-5,-10,-20,-50,-75,-100,-200,-500,-1000].
                A value of 0 denotes the surface and a value of -99999 denotes the bottom layer.
    doi_link  : doi link in string (required = False)
    varList   : list of variable names to extract (default=['temp','salt','u','v']).
                If both components of a vector pair (e.g. u/v, ubar/vbar, sustr/svstr, bustr/bvstr)
                are present, they are extracted together via get_uv() and rotated to east/north components.
                If only one component is present, it is extracted via get_var() and remains grid-aligned.
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

    for file in fname_in:
        print('')
        print(f'Opening: {file}')
        print("Extracting the model output variables we need")

        scalars, pairs = _partition_vars(varList)
        datasets = []
        for var in scalars:
            datasets.append(post.get_var(file, var, grdname=grdname, Yorig=Yorig, level=depths))
        for var_u, var_v in pairs:
            datasets.append(post.get_uv(file, grdname=grdname, Yorig=Yorig, level=depths, var_u=var_u, var_v=var_v))
        ds_all = xr.merge(datasets, compat='override')

        ds_all.attrs["title"] = "Regridded CROCO output created by the regrid_tier2 function"
        ds_all.attrs["source"] = file
        ds_all.attrs["history"] = "Created on " + datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S")
        ds_all.attrs["conventions"] = "CF-1.8"
        ds_all.attrs['references'] = 'Project: Sustainable Ocean Modelling Initiative: a South AfricaN Approach (SOMISANA; https://somisana.ac.za/); Tools: Regridding Code (https://github.com/SAEON/somisana-croco)'

        if doi_link is not None: ds_all.attrs.update({"doi" :f"https://doi.org/{doi_link}"})

        # build encoding dynamically based on variables present in ds_all
        encoding = {}
        for var in ds_all.data_vars:
            encoding[var] = {"dtype": "float32"}
        for var in ['h', 'mask', 'lon_rho', 'lat_rho']:
            if var in ds_all:
                encoding[var] = {"dtype": "float32"}
        if 'depth' in ds_all.variables:
            encoding['depth'] = {"dtype": "float32"}
        encoding["time"] = {"units": f"seconds since {Yorig}-01-01 00:00:00",
                            "calendar": "standard",
                            "dtype": "i4"}

        # robust way of getting the file extension, including CROCO child domains e.g. *nc.2, *.nc.2 etc
        basename = os.path.basename(file)
        match = re.search(r"(\.nc\.\d+)$|(\.nc)$", basename)
        basename_no_extension = basename[:match.start()]
        extension = match.group(0)

        fname_out = os.path.abspath(os.path.join(dir_out, basename_no_extension + '_t2' + extension))

        # If the file already exists, we remove it to avoid permission errors.
        if os.path.exists(fname_out): os.remove(fname_out)
        ds_all.to_netcdf(fname_out, encoding=encoding, mode="w")

        subprocess.call(["chmod", "-R", "775", fname_out])

        print(f'Created: {fname_out}')

def regrid_tier3(fname_in, dir_out, Yorig=2000, doi_link=None, spacing=0.01,
                 varList=['temp','salt','u','v']):
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
    dir_out   : path to output directory (required = True)
    Yorig     : Origin year used in setting up CROCO time i.e. seconds since Yorig-01-01
    spacing   : constant horizontal grid spacing (in degrees) to be used for the horizontal interpolation of the output (type: str or float).
                If None, the default is 0.01.
    doi_link  : doi link in string (required = False)
    varList   : list of variable names to interpolate (default=['temp','salt','u','v']).
                Variables are read from the tier 2 input file and interpolated onto the regular grid.
                3D variables (with a depth dimension) and 2D variables are handled automatically.

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

    for file in fname_in:
        print('')
        print(f'Opening: {file}')
        print("Extracting the tier 2 re-gridded model output")

        ds = xr.open_dataset(file)
        Nt = ds.time.size
        has_depth = 'depth' in ds.dims
        Nz = ds.depth.size if has_depth else 0
        lon_rho_1d = np.ravel(ds.lon_rho.values)
        lat_rho_1d = np.ravel(ds.lat_rho.values)

        # classify variables as 2D or 3D based on whether they have a depth dimension
        vars_3d = [v for v in varList if v in ds and 'depth' in ds[v].dims]
        vars_2d = [v for v in varList if v in ds and 'depth' not in ds[v].dims]

        print("Generating the regular horizontal output grid")
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

        print("Interpolating the model output onto the regular horizontal output grid")

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

        # zeta is always interpolated as a grid variable
        zeta_out = []
        # dicts to accumulate dask arrays per variable
        out_2d_time = {v: [] for v in vars_2d}
        out_3d_time = {v: [] for v in vars_3d}

        for t in np.arange(Nt):
            zeta_out.append(
                da.from_delayed(
                    compute_2d_chunk(t, ds.zeta.values),
                    shape=(Nlat, Nlon),
                    dtype=float,
                )
            )

            for v in vars_2d:
                out_2d_time[v].append(
                    da.from_delayed(
                        compute_2d_chunk(t, ds[v].values),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )

            if vars_3d:
                out_3d_depth = {v: [] for v in vars_3d}
                for n in np.arange(Nz):
                    for v in vars_3d:
                        out_3d_depth[v].append(
                            da.from_delayed(
                                compute_3d_chunk(t, ds[v].values, n),
                                shape=(Nlat, Nlon),
                                dtype=float,
                            )
                        )
                # Stack the depth dimension and append to the time list
                for v in vars_3d:
                    out_3d_time[v].append(da.stack(out_3d_depth[v], axis=0))

        # Stack the time dimension
        zeta_out = da.stack(zeta_out, axis=0)
        for v in vars_2d:
            out_2d_time[v] = da.stack(out_2d_time[v], axis=0)
        for v in vars_3d:
            out_3d_time[v] = da.stack(out_3d_time[v], axis=0)

        # Create new xarray dataset
        print("Generating dataset")
        data_vars = {
            "zeta": xr.Variable(
                ["time", "latitude", "longitude"],
                zeta_out,
                ds.zeta.attrs,
            ),
        }
        for v in vars_2d:
            data_vars[v] = xr.Variable(
                ["time", "latitude", "longitude"],
                out_2d_time[v],
                ds[v].attrs,
            )
        for v in vars_3d:
            data_vars[v] = xr.Variable(
                ["time", "depth", "latitude", "longitude"],
                out_3d_time[v],
                ds[v].attrs,
            )

        coords = {
            "longitude": xr.Variable(
                ["longitude"],
                lon_out,
                ds.lon_rho.attrs,
            ),
            "latitude": xr.Variable(
                ["latitude"],
                lat_out,
                ds.lat_rho.attrs,
            ),
            "time": xr.Variable(
                ["time"],
                ds.time.values,
                ds.time.attrs,
            ),
        }
        if has_depth:
            coords["depth"] = xr.Variable(
                ["depth"],
                ds.depth.values,
                ds.depth.attrs,
            )

        data_out = xr.Dataset(data_vars=data_vars, coords=coords)

        data_out.attrs["title"] = "Regridded CROCO output created by the regrid_tier3 function"
        data_out.attrs["source"] = file
        data_out.attrs["history"] = "Created on " + datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S")
        data_out.attrs["conventions"] = "CF-1.8"
        data_out.attrs['references'] = 'Project: Sustainable Ocean Modelling Initiative: a South AfricaN Approach (SOMISANA; https://somisana.ac.za/); Tools: Regridding Code (https://github.com/SAEON/somisana-croco)'
        if doi_link is not None: data_out.attrs.update({"doi" :f"https://doi.org/{doi_link}"})

        # Explicitly set chunk sizes of some dimensions
        chunksizes = {
            "time": ds.time.size,
        }
        if has_depth:
            chunksizes["depth"] = 1

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
        encoding["time"] = {"units": f"seconds since {Yorig}-01-01 00:00:00",
                            "calendar": "standard",
                            "dtype": "i4"}
        encoding['latitude'] = {"dtype": "float32"}
        encoding['longitude'] = {"dtype": "float32"}
        if has_depth:
            encoding['depth'] = {"dtype": "float32"}

        # robust way of getting the file extension, including CROCO child domains e.g. *nc.2, *.nc.2 etc
        basename = os.path.basename(file)
        match = re.search(r"(\.nc\.\d+)$|(\.nc)$", basename)
        basename_no_extension = basename[:match.start()]
        extension = match.group(0)
        basename_no_extension = basename_no_extension.replace('_t2', '_t3')

        fname_out = os.path.abspath(os.path.join(dir_out, basename_no_extension + extension))

        # If the file already exists, we remove it to avoid permission errors.
        if os.path.exists(fname_out): os.remove(fname_out)

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
        print(f'Created: {fname_out}')

if __name__ == "__main__":
    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg.nc'
    dir_out = '/home/g.rautenbach/Data/models/sa_southeast/'
    doi = '10.15493/SOMISANA.26032025'
    Yorig=2000
    
    # fname_in = '/home/gfearon/code/somisana-croco/configs/test_02/croco_v1.3.1/C04_I01_GLORYS_ERA5/output/croco_avg_Y2019M05.nc'
    # grdname = '/home/gfearon/croco_grd.nc'
    # dir_out = '/home/gfearon/test_cfcompliance/'

    # regrid_tier1(fname_in, dir_out, grdname=grdname,Yorig=Yorig, doi_link=doi)
    regrid_tier2(fname_in, dir_out, Yorig=Yorig, doi_link=doi, depths=[0,-5,-10])
    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg_t2.nc'
    # regrid_tier3(fname_in, dir_out, Yorig=Yorig, doi_link=doi, spacing=0.05)
