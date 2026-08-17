import crocotools_py.postprocess as post
import numpy as np
import xarray as xr
import os,sys,shutil
from datetime import datetime
#import pandas as pd
import matplotlib.path as mplPath
from dask import delayed
import dask.array as da
from scipy.interpolate import griddata
from scipy.spatial import Delaunay, cKDTree
from glob import glob
import subprocess
import dask
#import dask.array as da
import re

import crocotools_py.define_attrs as cf_attrs

# Discrete/categorical variables that must use nearest-neighbor vertical
# interpolation instead of linear - linear interpolation of a category flag
# (e.g. MHW/MCS category) produces meaningless fractional values and can
# blend fill values into valid data near mask boundaries.
CATEGORICAL_VARS = {'category'}


# Global attributes of the regridded output
#
# All three tiers describe themselves with define_attrs.global_attrs(), the same
# builder the products file uses, so a tier file carries the same descriptive
# header whether it was regridded from raw CROCO output or from a products file.
#
# The one thing that has to travel from the input is 'history': a products file
# records the climatology it was compared against, the MHW/MCS persistence
# criterion, the stratification depths and so on, and without chaining that
# forward a tier file cannot say where its variables came from. Coverage is not
# carried - global_attrs() derives it from the dataset being written, which is
# what keeps tier 3 honest after it builds its own regular grid.
#
# Anything specific to a single variable stays on that variable. get_var() and
# get_uv() carry those attributes through the interpolation, so the titles below
# are deliberately about the file, not about its contents.
TIER_TITLES = {
    1: 'SOMISANA CROCO output regridded to the rho grid (tier 1)',
    2: 'SOMISANA CROCO output on constant depth levels (tier 2)',
    3: 'SOMISANA CROCO output on a regular horizontal grid (tier 3)',
}

def source_history(file):
    """
    The 'history' of an input file, for global_attrs() to chain this step onto.

    Reads the header only, and tolerates anything unreadable by returning an
    empty string - the provenance chain is worth having but never worth failing
    a regrid over. Raw CROCO output has no 'history', so regridding it simply
    starts one.
    """
    try:
        with xr.open_dataset(file, decode_times=False) as ds_src:
            return str(ds_src.attrs.get('history', ''))
    except Exception as err:
        print(f'  could not read the global attributes of {file} ({err}) - '
              'the output will not carry its provenance')
        return ''

def regrid_tier1(fname_in, dir_out, grdname=None, Yorig=2000, doi_link=None,
                 varList=['temp','salt','u','v'], compress=False):
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
    compress  : deflate the data variables (default False, i.e. as written up to now).
                The file stays NETCDF4 either way - this only turns on the internal
                HDF5 filter, so the values are unchanged and any reader, THREDDS
                included, decompresses it transparently. Roughly halves the file at
                a cost of a few tens of seconds per file, which is worth it for a
                long hindcast being written across a network mount, and not worth
                bothering with for a single operational forecast.
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
    
        scalars, pairs = post.partition_vars(varList)
        datasets = []
        for var in scalars:
            datasets.append(post.get_var(file, var, grdname=grdname, Yorig=Yorig))
        for var_u, var_v in pairs:
            datasets.append(post.get_uv(file, grdname=grdname, Yorig=Yorig, var_u=var_u, var_v=var_v))
        ds_all = xr.merge(datasets, compat='override')
        
        ds_all.attrs = cf_attrs.global_attrs(
            title=TIER_TITLES[1], source=file,
            ds=ds_all, src_history=source_history(file), doi_link=doi_link,
            action=f'regrid_tier1: u/v regridded to the rho grid and rotated to '
                   f'east/north, sigma level depths added, from {file}')

        # build encoding dynamically based on variables present in ds_all
        encoding = {}
        for var in ds_all.data_vars:
            encoding[var] = {"dtype": "float32"}
            # the encoding sets the on-disk type, so flag_values and friends have
            # to be re-typed to match it or a float32 variable ends up describing
            # itself with attributes of whatever type it happened to have in
            # memory - see postprocess.retype_attrs()
            ds_all[var].attrs = post.retype_attrs(ds_all[var].attrs, 'float32')
        for var in ['h', 'mask', 'lon_rho', 'lat_rho']:
            if var in ds_all:
                encoding[var] = {"dtype": "float32"}
        if 'depth' in ds_all.variables:  # allow for both cases where depth is or isn't in the dataset (e.g. if a surface file is used)
            encoding['depth'] = {"dtype": "float32"}
        if compress:
            # deflate everything already in the encoding dict - the time axis is
            # added below and is far too small to be worth compressing
            for enc in encoding.values():
                enc.update({"zlib": True, "complevel": 2})
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

        scalars, pairs = post.partition_vars(varList)
        datasets = []
        for var in scalars:
            interp_method = 'nearest' if var in CATEGORICAL_VARS else 'linear'
            datasets.append(post.get_var(file, var, grdname=grdname, Yorig=Yorig, level=depths, interp_method=interp_method))
        for var_u, var_v in pairs:
            datasets.append(post.get_uv(file, grdname=grdname, Yorig=Yorig, level=depths, var_u=var_u, var_v=var_v))
        ds_all = xr.merge(datasets, compat='override')

        ds_all.attrs = cf_attrs.global_attrs(
            title=TIER_TITLES[2], source=file,
            ds=ds_all, src_history=source_history(file), doi_link=doi_link,
            action=f'regrid_tier2: interpolated to depth levels '
                   f'{", ".join(str(d) for d in depths)} m, from {file}')

        # build encoding dynamically based on variables present in ds_all
        encoding = {}
        for var in ds_all.data_vars:
            encoding[var] = {"dtype": "float32"}
            # the encoding sets the on-disk type, so flag_values and friends have
            # to be re-typed to match it or a float32 variable ends up describing
            # itself with attributes of whatever type it happened to have in
            # memory - see postprocess.retype_attrs()
            ds_all[var].attrs = post.retype_attrs(ds_all[var].attrs, 'float32')
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
                 varList=['temp','salt','u','v'], output_format='netcdf',
                 method='nearest'):
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
    output_format : 'netcdf' (default) or 'zarr'.
                When 'zarr', writes a compressed Zarr store (a directory) instead of a NetCDF file.
                Zarr output uses Blosc/Zstd compression and chunking tuned for per-map access
                (time=1, depth=1, full lat/lon plane per chunk) — well suited for web
                visualisation pipelines that read one forecast time step at a time.
                The output path is derived from the input by replacing '.nc' with '.zarr'
                (e.g. 'croco_avg_t2.nc.2' -> 'croco_avg_t3.zarr.2').
    method    : interpolation method passed to scipy.interpolate.griddata
                (default='nearest'). Supported values: 'nearest', 'linear', 'cubic'.

    CAREFUL! tier3 output is useful for website visualisation (it's intended use),
    but don't use it for reasearch/analysis as it's interpolated, so can be at a
    totally different resolution to the native CROCO grid

    '''

    if output_format not in ('netcdf', 'zarr'):
        raise ValueError(f"output_format must be 'netcdf' or 'zarr', got '{output_format}'")

    if method not in ('nearest', 'linear', 'cubic'):
        raise ValueError(f"method must be 'nearest', 'linear' or 'cubic', got '{method}'")

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
        xi_flat = np.column_stack([np.asarray(lon_out_grd).ravel(),
                                   np.asarray(lat_out_grd).ravel()])
        mask_out = poly_boundary.contains_points(xi_flat).reshape(Nlat, Nlon).astype(float)

        # Precompute the spatial interpolation structure ONCE. The source mesh
        # (lonlat_input) and the output points (xi_flat) are identical for every
        # (time, depth, variable) tile, so the expensive triangulation/KDTree
        # build and the per-output-point simplex/neighbor lookup are hoisted
        # out of the per-chunk hot loop. Each compute_*_chunk call collapses
        # to a fancy-indexed gather + (for linear) weighted sum.
        if method == 'linear':
            print("Precomputing Delaunay triangulation and barycentric weights")
            tri = Delaunay(lonlat_input)
            simplex = tri.find_simplex(xi_flat)             # (Nout,), -1 outside hull
            valid = simplex >= 0
            safe_simplex = np.where(valid, simplex, 0)
            vtx = tri.simplices[safe_simplex]               # (Nout, 3)
            T = tri.transform[safe_simplex]                 # (Nout, 3, 2)
            delta = xi_flat - T[:, 2, :]                    # (Nout, 2)
            bary = np.einsum('nij,nj->ni', T[:, :2, :], delta)
            weights = np.column_stack([bary, 1.0 - bary.sum(axis=1)])  # (Nout, 3)
            weights[~valid] = 0.0
            valid_2d = valid.reshape(Nlat, Nlon)

            def _interp_field(values_flat):
                gathered = values_flat[vtx]                  # (Nout, 3)
                out = np.einsum('ni,ni->n', gathered, weights).reshape(Nlat, Nlon)
                out = np.where(valid_2d, out, np.nan)
                return out * mask_out / mask_out
        elif method == 'nearest':
            print("Precomputing nearest-neighbour index")
            kd = cKDTree(lonlat_input)
            _, nn_idx = kd.query(xi_flat, k=1)               # (Nout,)

            def _interp_field(values_flat):
                out = values_flat[nn_idx].reshape(Nlat, Nlon)
                return out * mask_out / mask_out
        else:
            # cubic: scipy has no clean precomputed form, fall back to per-call griddata
            def _interp_field(values_flat):
                out = griddata(lonlat_input, values_flat,
                               (lon_out_grd, lat_out_grd), method)
                return out * mask_out / mask_out

        print("Interpolating the model output onto the regular horizontal output grid")

        @delayed
        def compute_2d_chunk(t, variable):
            return _interp_field(np.ravel(variable[t, ::]))

        @delayed
        def compute_3d_chunk(t, variable, n):
            return _interp_field(np.ravel(variable[t, n, ::]))

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

        # source_history() rather than the already-open ds.attrs, so that all
        # three tiers treat their input's header identically
        data_out.attrs = cf_attrs.global_attrs(
            title=TIER_TITLES[3], source=file,
            ds=data_out, src_history=source_history(file), doi_link=doi_link,
            action=f'regrid_tier3: interpolated onto a regular {spacing} degree '
                   f'grid using {method} interpolation, from {file}')

        # Explicitly set chunk sizes of some dimensions.
        # NetCDF default: time=full, depth=1 (good for point-timeseries analysis).
        # Zarr default:   time=1,    depth=1 (good for reading one map at a time,
        #                 e.g. web visualisation pipelines).
        chunksizes = {}
        if output_format == 'zarr':
            chunksizes["time"] = 1
        else:
            chunksizes["time"] = ds.time.size
        if has_depth:
            chunksizes["depth"] = 1

        # For data_vars, set chunk sizes for each dimension
        # This is either the override specified in "chunksizes"
        # or the length of the dimension
        default_chunksizes = {dim: len(data_out[dim]) for dim in data_out.dims}

        # everything below is written as float32, so the type-sensitive
        # attributes have to say so too - see postprocess.retype_attrs()
        for var in data_out.data_vars:
            data_out[var].attrs = post.retype_attrs(data_out[var].attrs, 'float32')

        if output_format == 'zarr':
            from numcodecs import Blosc
            compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.SHUFFLE)
            encoding = {
                var: {
                    "dtype": "float32",
                    "chunks": [chunksizes.get(dim, default_chunksizes[dim]) for dim in data_out[var].dims],
                    "compressor": compressor,
                }
                for var in data_out.data_vars
            }
        else:
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
        if output_format == 'zarr':
            # e.g. '.nc' -> '.zarr', '.nc.2' -> '.zarr.2'
            extension = extension.replace('.nc', '.zarr')

        fname_out = os.path.abspath(os.path.join(dir_out, basename_no_extension + extension))

        # If the output already exists, we remove it to avoid permission errors.
        # Zarr is a directory; NetCDF is a single file.
        if os.path.exists(fname_out):
            if output_format == 'zarr':
                shutil.rmtree(fname_out)
            else:
                os.remove(fname_out)

        if output_format == 'zarr':
            print("Generating Zarr store")
            write_op = data_out.to_zarr(
                    fname_out,
                    encoding=encoding,
                    mode="w",
                    compute=False,
                    )
            print("Writing Zarr store")
            dask.compute(write_op)
        else:
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
