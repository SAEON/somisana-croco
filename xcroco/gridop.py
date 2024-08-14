import numpy as np
import xarray as xr
from functools import partial
import pandas as pd
from xgcm import Grid
import intake

import time
import os.path as path
import itertools
from collections import OrderedDict

from math import radians, cos, sin, asin, sqrt

import inout as io

def get_cs(model, ds, gd, vgrid):
    """get croco vertical grid  stretching
    https://www.myroms.org/wiki/Vertical_S-coordinate
    Args:
        model (class): classe of the model
        ds (DataSet): input dataset from the history files
        gd (DataSet): DataSet of the grid file
        vgrid (character): type of metrics ('r'': rho level, 'w': w level)

    Returns:
        DataArray: vertical grid  stretching
    """
    # search vtransform 
    if 'vtransform' in gd:
        vtransform = gd.vtransform.values
    elif 'vtransform' in ds:
        vtransform = ds.vtransform.values
    elif 'vtransform' in gd.attrs:
        vtransform = gd.attrs['vtransform']
    elif 'vtransform' in ds.attrs:
        vtransform = ds.attrs['vtransform']
    else:
        print("Can't find vtransform neither in filename nor in gridname  ")
        return None
    
    # search theta_s 
    if io.find_var(model,'theta_s',ds,gd) is not None: 
        theta_s = io.find_var(model,'theta_s',ds,gd).values
    else:
        print("Can't find theta_s neither in filename nor in gridname  ")
        return None
    # search theta_b 
    if io.find_var(model,'theta_b',ds,gd) is not None: 
        theta_b = io.find_var(model,'theta_b',ds,gd).values
    else:
        print("Can't find theta_b neither in filename nor in gridname  ")
        return None
    
    sc = ds.sc_r.values if vgrid=='r' else ds.sc_w.values
    
    # New coordinates
    if vtransform == 2 or vtransform.lower()=='new':
        if theta_s>0:
            csf = (1-np.cosh(theta_s*sc)) / (np.cosh(theta_s)-1.)
        else:
            csf = sc**2
        if theta_b>0:
            cs = (np.exp(theta_s*csf)-1.) / (1-np.exp(-theta_b))
        else:
            cs = csf
    # Old coordinates
    else:
        cs = (1-theta_b)*np.sinh(theta_s*sc)/np.sinh(theta_s) \
             + theta_b*0.5 \
               *(np.tanh((sc+0.5)*theta_s)/np.tanh(0.5*theta_s)-1.)
    return cs

def add_grid(model, gridname, grid_metrics=1, remove_ghost_pts=True,
             xperiodic=False, yperiodic=False):
    """from the gridname file, add the grid to the dataset and compute the XGCM grid

    Args:
        model (class): classe of the model
        ds (DataSet): input dataset from the history files
        gridname (string): name of the grid file
        grid_metrics (int, optional): type of metrics (0: no metrics, 1: horizontal, 2: 1+vertical)
        remove_ghost_pts (bool, optional): remove ghost points. Defaults to True.
        xperiodic (bool, optional): periodic in x. Defaults to False.
        yperiodic (bool, optional): periodic in y. Defaults to False.

    Returns:
        DataSet: the input dataset with the grid inside
        XGCM grid: the XGCM grid associated to the dataset
    """
        
    # open grid file
    try : 
        gd = xr.open_zarr(gridname).squeeze()
    except:
        try : 
            gd = xr.open_dataset(gridname).squeeze() 
        except :
            print('add_grid: unknown format for grid : only Netcdf or Zarr')
            
    # Rename variable according model
    gd = adjust_grid(model, gd)
    ds = model.ds

    if io.find_var(model,'hc',ds,gd) is not None: ds['hc'] = io.find_var(model,'hc',ds,gd)
    if io.find_var(model,'h',ds,gd) is not None: ds['h'] = io.find_var(model,'h',ds,gd)
    if io.find_var(model,'pm',ds,gd) is not None: ds['pm']   = io.find_var(model,'pm',ds,gd)
    if io.find_var(model,'pn',ds,gd) is not None: ds['pn']   = io.find_var(model,'pn',ds,gd)
    try:
        N = ds.dims['s']
        if 'sc_r' not in ds:
            if io.find_var(model,'sc_r',ds,gd) is not None: 
                ds['sc_r'] = io.find_var(model,'sc_r',ds,gd)
            else:
                ds['sc_r'] = xr.DataArray(np.arange(-1.+1./(N+1),0., 1/(N+1)), dims='s')  
        if 'sc_w' not in ds:      
            if io.find_var(model,'sc_w',ds,gd) is not None: 
                ds['sc_w'] = io.find_var(model,'sc_w',ds,gd)
            else:
                ds['sc_w'] = xr.DataArray(np.arange(-1.,0., 1/(N+1)), dims='s_w')
        if 'Cs_r' not in ds:
            if  io.find_var(model,'Cs_r',ds,gd) is not None: 
                ds['Cs_r'] = io.find_var(model,'Cs_r',ds,gd)
            else:
                ds['Cs_r'] = get_cs(model,ds, gd, 'r')
        if 'Cs_w' not in ds:
            if  io.find_var(model,'Cs_w',ds,gd) is not None: 
                ds['Cs_w'] = io.find_var(model,'Cs_w',ds,gd)
            else:
                ds['Cs_w'] = get_cs(model,ds, gd, 'w')
    except:
        pass        
        
    if io.find_var(model,'angle',ds,gd) is not None: ds['angle'] = io.find_var(model,'angle',ds,gd)
    if io.find_var(model,'mask',ds,gd) is not None: ds['mask'] = io.find_var(model,'mask',ds,gd)
    if io.find_var(model,'lon',ds,gd) is not None: ds['lon'] = io.find_var(model,'lon',ds,gd)
    if io.find_var(model,'lat',ds,gd) is not None: ds['lat'] = io.find_var(model,'lat',ds,gd)
    if io.find_var(model,'f',ds,gd) is not None: ds['f'] = io.find_var(model,'f',ds,gd)
    if io.find_var(model,'rho0',ds,gd) is not None: ds['rho0'] = io.find_var(model,'rho0',ds,gd)
    if io.find_var(model,'g',ds,gd) is not None: ds['g'] = io.find_var(model,'g',ds,gd)
    
    
    # coords = [c for c in ds.coords if c not in ['t','s','s_w']]
    coords = [c for c in ds.coords if c in ['t','s','s_w','lat','lon']]
    ds = ds.reset_coords()
    ds = ds.set_coords(coords)
        
    # remove ghost points
    if remove_ghost_pts:
        ds = remove_ghost_points(ds, xperiodic=xperiodic, yperiodic=yperiodic)
    model.ds = ds
    
    # On crée la grille xgcm
    ds, grid = xgcm_grid(model, grid_metrics=grid_metrics, 
                         xperiodic=xperiodic, yperiodic=yperiodic)
    
    return ds, grid

def remove_ghost_points(ds, xperiodic=False, yperiodic=False):
    """
    Remove ghost points from the DataSet
    Args:
        ds (DataSet): input dataset from the history files
        xperiodic (bool, optional): periodic in x. Defaults to False.
        yperiodic (bool, optional): periodic in y. Defaults to False.

    Returns:
        DataSet: the input dataset without any ghost points
    """
    ds = ds.isel(x=slice(1,-1),y=slice(1,-1))
    if xperiodic:
        ds = ds.isel(x_u=slice(0,-1))
    if yperiodic:
        ds = ds.isel(y_v=slice(0,-1))
    return ds

def xgcm_grid(model, grid_metrics=1, xperiodic=False, yperiodic=False):
    """
    Create the xgcm grid of the dataset.
    Args:
        model: (Model class) the model class
        grid_metrics: (integer) 0:no metrics, 1:horizontal metrics, 2:hor + vert metrics
        xperiodic: True if ds periodic in x
        yperiodic: True if ds periodic in y

    Return:
        DataSet : the dataset with the news metrics
        XGCM grid: the xgcm grid of the dataset

    """
        
    # Create xgcm grid without metrics
    coords={}
    if all(d in model.ds.dims for d in ['x','x_u']):
        if xperiodic:
            coords.update({'x': {'center':'x', 'left':'x_u'}})
        else:
            coords.update({'x': {'center':'x', 'outer':'x_u'}})
    if all(d in model.ds.dims for d in ['y','y_v']):
        if yperiodic:
            coords.update({'y': {'center':'y', 'left':'y_v'}} )
        else:
            coords.update({'y': {'center':'y', 'outer':'y_v'}} )
    if 's' in model.ds.dims:
        coords.update({'z': {'center':'s', 'outer':'s_w'}})
        
    grid = Grid(model.ds, 
              coords=coords,
              periodic=False,
              boundary='extend')
    
    if grid_metrics==0:           
        model.xgrid = grid
        return model.ds, grid
    
    # compute horizontal coordinates

    ds = model.ds
    if 'x_u' in ds.dims:
        ds['lon_u'] = grid.interp(ds.lon,'x')
        ds['lat_u'] = grid.interp(ds.lat,'x')
    if 'y_v' in ds.dims:
        ds['lon_v'] = grid.interp(ds.lon,'y')
        ds['lat_v'] = grid.interp(ds.lat,'y')
    if 'x_u' in ds.dims and 'y_v' in ds.dims: 
        ds['lon_p'] = grid.interp(ds.lon_v,'x')
        ds['lat_p'] = grid.interp(ds.lat_u,'y')
    _coords = [d for d in ds.data_vars.keys() if d.startswith(tuple(['lon','lat']))]
    ds = ds.set_coords(_coords)
    
    
    # add horizontal metrics for u, v and psi point
    if 'pm' in ds and 'pn' in ds:
        ds['dx'] = 1/ds['pm']
        ds['dy'] = 1/ds['pn']
    else: # backward compatibility, hack
        dlon = grid.interp(grid.diff(ds.lon,'x'),'x')
        dlat =  grid.interp(grid.diff(ds.lat,'y'),'y')
        ds['dx'], ds['dy'] = dll_dist(dlon, dlat, ds.lon, ds.lat)
    ds["dx_u"] = grid.interp(ds["dx"], "x")
    # ds["dy_u"] = grid.interp(ds["dy"], "x")
    # ds["dx_v"] = grid.interp(ds["dx"], "y")
    ds["dy_v"] = grid.interp(ds["dy"], "y")
    # ds["dx_psi"] = grid.interp(ds["dx_u"], "y")
    # ds["dy_psi"] = grid.interp(ds["dy_v"], "x")

    # add areas metrics for rho,u,v and psi points
    # ds['rAr'] = ds.dx_psi * ds.dy_psi
    # ds['rAu'] = ds.dx_v * ds.dy_v
    # ds['rAv'] = ds.dx_u * ds.dy_u
    # ds['rAf'] = ds.dx * ds.dy
    
    metrics={}    
    # if all(d in model.ds.dims for d in ['x','x_u']):
    # metrics.update({('x',): ['dx', 'dx_u', 'dx_v', 'dx_psi']})
    metrics.update({('x',): ['dx', 'dx_u']})
    # if all(d in model.ds.dims for d in ['y','y_v']):
    # metrics.update({('y',): ['dy', 'dy_u', 'dy_v', 'dy_psi']})
    metrics.update({('y',): ['dy', 'dy_v']})
    # if all(d in model.ds.dims for d in ['x','x_u','y','y_v']):
        # metrics.update({('x', 'y'): ['rAr', 'rAu', 'rAv', 'rAf']})
     
    if grid_metrics==1:
        # generate xgcm grid
        grid = Grid(ds,
                    coords=coords,
                    periodic=False,
                    metrics=metrics,
                    boundary='extend')
        model.xgrid = grid
        model.ds = ds
        return ds, grid
    
    # compute z coordinate at rho/w points
    if 'z_sfc' in [v for v in ds.data_vars] and \
       's' in [d for d in ds.dims.keys()] and \
        ds['s'].size>1:
        ds['is3D'] = True
        z = get_z(model, z_sfc=ds.z_sfc, xgrid=grid).fillna(0.)
        z_w = get_z(model, z_sfc=ds.z_sfc, xgrid=grid, vgrid='w').fillna(0.)
        ds['z'] = z
        ds['z_w'] = z_w
        ds['z_u'] = grid.interp(z,'x')
        ds['z_v'] = grid.interp(z,'y')
        ds['z_p'] = grid.interp(ds.z_u,'y')
        # set as coordinates in the dataset
        _coords = ['z','z_w','z_u','z_v','z_p']
        ds = ds.set_coords(_coords)
    else:
        ds['is3D'] = False

    # add vertical metrics for u, v, rho and psi points
    # if 'z' in [v for v in ds.coords]:
    if ds['is3D']:
        ds['dz'] = grid.interp(grid.diff(z,'z'),'z')
        ds['dz_w'] = grid.interp(grid.diff(z_w,'z'),'z')
        # ds['dz'] = grid.interp(grid.diff(ds.z,'z'),'z')
        # ds['dz_w'] = grid.interp(grid.diff(ds.z_w,'z'),'z')
        # ds['dz_u'] = grid.interp(grid.diff(ds.z_u,'z'),'z')
        # ds['dz_v'] = grid.interp(grid.diff(ds.z_v,'z'),'z')
        # ds['dz_p'] = grid.interp(grid.diff(ds.z_p,'z'),'z')
        
    # add coords and metrics for xgcm for the vertical direction
    # if 'z' in ds:
    if ds['is3D']:
# #             coords.update({'z': {'center':'s', 'outer':'s_w'}})
#             metrics.update({('z',): ['dz', 'dz_u', 'dz_v', 'dz_p', 'dz_w']}), # Z distances
        metrics.update({('z',): ['dz', 'dz_w']}), # Z distances
    # generate xgcm grid
    grid = Grid(ds,
                coords=coords,
                periodic=False,
                metrics=metrics,
                boundary='extend')

    model.xgrid = grid
    model.ds = ds

    return ds, grid

def fast_xgcm_grid(ds, grid_metrics=1, xperiodic=False, yperiodic=False):
    """
    Create the xgcm grid without computing any metrics. Just use those which are already
    in the dataset.
    Args:
        ds: (Xarray Dataset) the dataset to create the xgcm grid
        grid_metrics: (integer) 0:no metrics, 1:horizontal metrics, 2:hor + vert metrics
        xperiodic: True if ds periodic in x
        yperiodic: True if ds periodic in y

    Return:
        grid: the xgcm grid

    """
    
    # Create xgcm grid without metrics
    coords={}
    if all(d in ds.dims for d in ['x','x_u']):
        if xperiodic:
            coords.update({'x': {'center':'x', 'left':'x_u'}})
        else:
            coords.update({'x': {'center':'x', 'outer':'x_u'}})
    if all(d in ds.dims for d in ['y','y_v']):
        if yperiodic:
            coords.update({'y': {'center':'y', 'left':'y_v'}} )
        else:
            coords.update({'y': {'center':'y', 'outer':'y_v'}} )
    if all(d in ds.dims for d in ['s','s_w']):
        coords.update({'z': {'center':'s', 'outer':'s_w'}})
        
    grid = Grid(ds, 
              coords=coords,
              periodic=False,
              boundary='extend')

    if grid_metrics==0: return grid

    # set all lon/lat variables as coordinates
    _coords = [d for d in ds.data_vars.keys() if d.startswith(tuple(['lon','lat']))]
    ds = ds.set_coords(_coords)
             
    # set horizontal metrics
    # move horizontal metrics from global attributes to variables
    attrs = [k for k in ds.attrs.keys() if k.startswith(('dx','dy','rA'))]
    if attrs is not None:
        for k in attrs: ds[k] = ds.attrs[k]
    metrics={}   
    # add dx metrics
    if all(d in ds.dims for d in ['x','x_u']):
        dx = [v for v in ds.data_vars if v in ['dx','dx_u','dx_v','dx_psi']]
        metrics.update({('x',): dx})
    # add dy metrics
    if all(d in ds.dims for d in ['y','y_v']):
        dy = [v for v in ds.data_vars if v in ['dy','dy_u','dy_v','dy_psi']]
        metrics.update({('y',): dy})
    # add area metrics
    if all(d in ds.dims for d in ['x','x_u','y','y_v']):        
        rA = [v for v in ds.data_vars if v in ['rAr','rAu','rAv','rAf']]
        metrics.update({('x', 'y'): rA})
    
    if grid_metrics==1:
        # generate xgcm grid
        grid = Grid(ds,
                    coords=coords,
                    periodic=False,
                    metrics=metrics,
                    boundary='extend')
        return grid

    # Set z variables as coordinates
    _coords = [d for d in ds.data_vars.keys() if d in ['z','z_w','z_u','z_v','z_p']]
    ds = ds.set_coords(_coords)

    # add vertical metrics
    # move vertical metrics from global attributes to variables
    attrs = [k for k in ds.attrs.keys() if k.startswith(('dz'))]
    if attrs is not None:
        for k in attrs: ds[k] = ds.attrs[k]
    # add dz metrics
    if all(d in ds.dims for d in ['s','s_w']):        
        dz = [v for v in ds.data_vars if v in ['dz', 'dz_u', 'dz_v', 'dz_p', 'dz_w']]
        metrics.update({('z',): dz})
    
    # generate xgcm grid
    grid = Grid(ds,
                coords=coords,
                periodic=False,
                metrics=metrics,
                boundary='extend')

    return grid
    
def dll_dist(dlon, dlat, lon, lat):
    """
    Converts lat/lon differentials into distances in meters
    Args:
        dlon : xarray.DataArray longitude differentials
        dlat : xarray.DataArray latitude differentials
        lon : xarray.DataArray longitude values
        lat : xarray.DataArray latitude values
    Return:
        dx : xarray.DataArray distance inferred from dlon
        dy : xarray.DataArray distance inferred from dlat
    """
    distance_1deg_equator = 111000.0
    dx = dlon * np.cos(np.deg2rad(lat)) * distance_1deg_equator 
    dy = ((lon * 0) + 1) * dlat * distance_1deg_equator
    return dx, dy



def adjust_grid(model, ds):
    """
    Change the names in the dataset according to the model
    Args :
        model (Model class): Instance of the model class
        ds (Dataset): dataset to change
    Return :
        DataSet : changed dataset
    """
   
    for k,v in model.rename_vars.items():
        if k in ds or k in ds.dims.keys():
            if v in ds and k != v:
                ds = ds.drop(k)
            else:
                ds = ds.rename({k: v})
    # change names in attributes
    for k,v in model.rename_vars.items():
        if (k in ds.attrs and v not in ds.attrs):
            ds.attrs[v] = ds.attrs.pop(k)
    return ds
    

def get_spatial_dims(v):
    """Return an ordered dict of spatial dimensions in the s, y, x order
    Args :
        v (DataArray) : variable for which you have to guess the dimensions
    Return:
        Dictionary : ordered dimensions
    """
    if isinstance(v, xr.DataArray):
        dims = OrderedDict( (d, next((x for x in v.dims if x[0]==d), None))
                        for d in ['s','y','x'] )
    elif isinstance(v, xr.Dataset):
        # dims = OrderedDict( (d, next((x for x in v.dims if x==d), None))
        dims = OrderedDict( (d, [x for x in v.dims if x[0]==d])
                        for d in ['s','y','x'] ) 
        # convert empty list to None
        dims = {k: None if not d else d for k, d in dims.items() }
    else:
        print('get_spatial_dims: ERROR!!! the argument must be a DataArray or a Dataset')
    return dims


def get_spatial_coords(v):
    """Return an ordered dict of spatial coordinates in the z, lat, lon order
    Args:
        v (DataArray) : variable for which you have to guess the coordinates
    Return:
        Dictionary: ordered coordinates
    """
    coords = OrderedDict( (d, next((x for x in v.coords if x.startswith(d)), None))
                       for d in ['z','lat','lon'] )
    for k,c in coords.items():
        if c is not None and v.coords[c].size==1: coords[k]= None
    return coords


def order_dims(var):
    """Reorder var to typical dimensional ordering.

    Args:
        var (DataArray) : Variable to operate on.

    Returns
        DataArray : with dimensional order ['T', 'Z', 'Y', 'X'], or whatever subset of
    dimensions are present in var.

    Notes
    -----
    Do not consider previously-selected dimensions that are kept on as coordinates but
    cannot be transposed anymore. This is accomplished with `.reset_coords(drop=True)`.

    Examples
    --------
    >>> xroms.order(var)
    """

    return var.cf.transpose(
        *[
            # dim
            # for dim in ["T", "Z", "Y", "X"]
            # if dim in var.reset_coords(drop=True).cf.axes
            dim
            for dim in ["t", "s", "s_w", "y", "y_v" "x", "x_u"]
            if dim in var.dims
        ]
    )


def reorder_dims(da):    
    # reorder spatial dimensions and place them last
    sdims = list(get_spatial_dims(da).values())
    sdims = tuple(filter(None,sdims)) # delete None values
    reordered_dims = tuple(d for d in da.dims if d not in sdims) + sdims
    return da.transpose(*reordered_dims, transpose_coords=True)


def to_rho(v, grid, hboundary="extend", hfill_value=None):
    """Interpolate to rho horizontal grid

    Args:
        - v (DataArray): variable to interpolate
        - grid (xgcm.grid): grid object associated with v
        - hboundary (str, optional):
            From xgcm documentation:
            A flag indicating how to handle boundaries:
            - None: Do not apply any boundary conditions.
                    Raise an error if boundary conditions are required for the operation.
            - ‘fill’: Set values outside the array boundary to fill_value
                    (i.e. a Dirichlet boundary condition.)
            - ‘extend’: Set values outside the array to the nearest array value.
                    (i.e. a limited form of Neumann boundary condition.)
        -hfill_value: float, optional
            From xgcm documentation:
            The value to use in the boundary condition with `boundary='fill'`.

    Returns:
        DataArray: input variable interpolated on a rho horizontal point
    """
    vout = v.copy()
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    coords = {k: v for k, v in coords.items() if v is not None}

    if "x" not in v.dims and dims["x"] is not None and v[dims["x"]].size > 1:
        v = grid.interp(v, "x")
        for k, c in coords.items():
            v.coords[c] = grid.interp(vout.coords[c], "x")

    vout = v.copy()
    if "y" not in v.dims and dims["y"] is not None and v[dims["y"]].size > 1:
        v = grid.interp(v, "y")
        for k, c in coords.items():
            v.coords[c] = grid.interp(vout.coords[c], "y")

    return v


def to_u(v, grid, hboundary="extend", hfill_value=None):
    """Interpolate to u horizontal grid

    Args:
        - v (DataArray): variable to interpolate
        - grid (xgcm.grid): grid object associated with v
        - hboundary (str, optional):
            From xgcm documentation:
            A flag indicating how to handle boundaries:
            - None: Do not apply any boundary conditions.
                    Raise an error if boundary conditions are required for the operation.
            - ‘fill’: Set values outside the array boundary to fill_value
                    (i.e. a Dirichlet boundary condition.)
            - ‘extend’: Set values outside the array to the nearest array value.
                    (i.e. a limited form of Neumann boundary condition.)
        -hfill_value: float, optional
            From xgcm documentation:
            The value to use in the boundary condition with `boundary='fill'`.

    Returns:
        DataArray: input variable interpolated on a u horizontal point
    """
    vout = v.copy()
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    coords = {k: v for k, v in coords.items() if v is not None}

    if "x_u" not in v.dims and dims["x"] is not None and v[dims["x"]].size > 1:
        v = grid.interp(v, "x")
        for k, c in coords.items():
            v.coords[c] = grid.interp(vout.coords[c], "x")

    vout = v.copy()
    if "y" not in v.dims and dims["y"] is not None and v[dims["y"]].size > 1:
        v = grid.interp(v, "y")
        for k, c in coords.items():
            v.coords[c] = grid.interp(vout.coords[c], "y")
            
    coords = {k: v for k, v in coords.items() if v[0:3] in ['lon', 'lat']}
    for k,c in coords.items():
        v = v.rename({c:c[0:4]+"_u"})

    return v


def to_v(v, grid, hboundary="extend", hfill_value=None):
    """Interpolate to v horizontal grid

    Args:
        - v (DataArray): variable to interpolate
        - grid (xgcm.grid): grid object associated with v
        - hboundary (str, optional):
            From xgcm documentation:
            A flag indicating how to handle boundaries:
            - None: Do not apply any boundary conditions.
                    Raise an error if boundary conditions are required for the operation.
            - ‘fill’: Set values outside the array boundary to fill_value
                    (i.e. a Dirichlet boundary condition.)
            - ‘extend’: Set values outside the array to the nearest array value.
                    (i.e. a limited form of Neumann boundary condition.)
        -hfill_value: float, optional
            From xgcm documentation:
            The value to use in the boundary condition with `boundary='fill'`.

    Returns:
        DataArray: input variable interpolated on a v horizontal point
    """
    vout = v.copy()
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    coords = {k: v for k, v in coords.items() if v is not None}

    if "x" not in v.dims and dims["x"] is not None and v[dims["x"]].size > 1:
        v = grid.interp(v, "x")
        for k, c in coords.items():
            v.coords[c] = grid.interp(vout.coords[c], "x")

    vout = v.copy()
    if "y_v" not in v.dims and dims["y"] is not None and v[dims["y"]].size > 1:
        v = grid.interp(v, "y")
        for k, c in coords.items():
            v.coords[c] = grid.interp(vout.coords[c], "y")
            
    coords = {k: v for k, v in coords.items() if v[0:3] in ['lon', 'lat']}
    for k,c in coords.items():
        v = v.rename({c:c[0:4]+"_v"})

    return v


def to_psi(v, grid, hboundary="extend", hfill_value=None):
    """Interpolate to psi horizontal grid

    Args:
        - v (DataArray): variable to interpolate
        - grid (xgcm.grid): grid object associated with v
        - hboundary (str, optional):
            From xgcm documentation:
            A flag indicating how to handle boundaries:
            - None: Do not apply any boundary conditions.
                    Raise an error if boundary conditions are required for the operation.
            - ‘fill’: Set values outside the array boundary to fill_value
                    (i.e. a Dirichlet boundary condition.)
            - ‘extend’: Set values outside the array to the nearest array value.
                    (i.e. a limited form of Neumann boundary condition.)
        -hfill_value: float, optional
            From xgcm documentation:
            The value to use in the boundary condition with `boundary='fill'`.

    Returns:
        DataArray: input variable interpolated on a psi horizontal point
    """
    vout = v.copy()
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    coords = {k: v for k, v in coords.items() if v is not None}

    if "x_u" not in v.dims and dims["x"] is not None and v[dims["x"]].size > 1:
        v = grid.interp(v, "x")
        for k, c in coords.items():
            v.coords[c] = grid.interp(vout.coords[c], "x")

    vout = v.copy()
    if "y_v" not in v.dims and dims["y"] is not None and v[dims["y"]].size > 1:
        v = grid.interp(v, "y")
        for k, c in coords.items():
            v.coords[c] = grid.interp(vout.coords[c], "y")
            
    coords = {k: v for k, v in coords.items() if v[0:3] in ['lon', 'lat']}
    for k,c in coords.items():
        v = v.rename({c:c[0:4]+"_p"})
    return v


def to_s_rho(v, grid, vboundary="extend", vfill_value=None):
    """Interpolate to rho vertical grid

    Args:
        - v (DataArray): variable to interpolate
        - grid (xgcm.grid): grid object associated with v
        - hboundary (str, optional):
            From xgcm documentation:
            A flag indicating how to handle boundaries:
            - None: Do not apply any boundary conditions.
                    Raise an error if boundary conditions are required for the operation.
            - ‘fill’: Set values outside the array boundary to fill_value
                    (i.e. a Dirichlet boundary condition.)
            - ‘extend’: Set values outside the array to the nearest array value.
                    (i.e. a limited form of Neumann boundary condition.)
        -hfill_value: float, optional
            From xgcm documentation:
            The value to use in the boundary condition with `boundary='fill'`.

    Returns:
        DataArray: input variable interpolated on a rho vertical level
    """
    vout = v.copy()
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    coords = {k: v for k, v in coords.items() if v is not None}

    interp = False
    if "s" not in v.dims and dims["s"] is not None and v[dims["s"]].size > 1:
        interp = True
        v = grid.interp(v, "z")
        if "lon" in coords.keys():
            v.coords[coords["lon"]] = vout.coords[coords["lon"]]
        if "lat" in coords.keys():
            v.coords[coords["lat"]] = vout.coords[coords["lat"]]
        if "z" in coords.keys():
            v.coords["z"] = grid.interp(vout.coords[coords["z"]], "z")

    return v


def to_s_w(v, grid, vboundary="extend", vfill_value=None):
    """Interpolate to w vertical grid

    Args:
        - v (DataArray): variable to interpolate
        - grid (xgcm.grid): grid object associated with v
        - hboundary (str, optional):
            From xgcm documentation:
            A flag indicating how to handle boundaries:
            - None: Do not apply any boundary conditions.
                    Raise an error if boundary conditions are required for the operation.
            - ‘fill’: Set values outside the array boundary to fill_value
                    (i.e. a Dirichlet boundary condition.)
            - ‘extend’: Set values outside the array to the nearest array value.
                    (i.e. a limited form of Neumann boundary condition.)
        -hfill_value: float, optional
            From xgcm documentation:
            The value to use in the boundary condition with `boundary='fill'`.

    Returns:
        DataArray: input variable interpolated on a w vertical level
    """
    vout = v.copy()
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    coords = {k: v for k, v in coords.items() if v is not None}

    interp = False
    if "s_w" not in v.dims and dims["s"] is not None and v[dims["s"]].size > 1:
        interp = True
        v = grid.interp(v, "z")
        if "lon" in coords.keys():
            v.coords[coords["lon"]] = vout.coords[coords["lon"]]
        if "lat" in coords.keys():
            v.coords[coords["lat"]] = vout.coords[coords["lat"]]
        if "z" in coords.keys():
            v.coords["z_w"] = grid.interp(vout.coords[coords["z"]], "z")

    return v


def to_grid_point(
    var,
    xgrid,
    hcoord=None,
    vcoord=None,
    hboundary="extend",
    hfill_value=None,
    vboundary="extend",
    vfill_value=None,
    attrs=None,
):
    """Interpolate to a new grid point.

    Args
        var: DataArray or ndarray
            Variable to operate on.
        xgrid: xgcm.grid
            Grid object associated with var
        hcoord: string, optional.
            Name of horizontal grid to interpolate output to.
            Options are 'r', 'rho','p', 'psi', 'u', 'v'.
        vcoord: string, optional.
            Name of vertical grid to interpolate output to.
            Options are 's_rho', 's_w', 'rho', 'r', 'w'.
        hboundary: string, optional
            From xgcm documentation:
            A flag indicating how to handle boundaries:
            * None:  Do not apply any boundary conditions. Raise an error if
            boundary conditions are required for the operation.
            * 'fill':  Set values outside the array boundary to fill_value
            (i.e. a Neumann boundary condition.)
            * 'extend': Set values outside the array to the nearest array
            value. (i.e. a limited form of Dirichlet boundary condition.
        hfill_value: float, optional
            From xgcm documentation:
            The value to use in the boundary condition with `boundary='fill'`.
        sboundary: string, optional
            From xgcm documentation:
            A flag indicating how to handle boundaries:
            * None:  Do not apply any boundary conditions. Raise an error if
            boundary conditions are required for the operation.
            * 'fill':  Set values outside the array boundary to fill_value
            (i.e. a Neumann boundary condition.)
            * 'extend': Set values outside the array to the nearest array
            value. (i.e. a limited form of Dirichlet boundary condition.
        sfill_value: float, optional
            From xgcm documentation:
            The value to use in the boundary condition with `boundary='fill'`.

    Returns
        DataArray or ndarray interpolated onto hcoord horizontal and vcoord
        vertical point.

    Notes
    -----
    If var is already on selected grid, nothing happens.

    Examples
    --------
    >>> to_grid_point(ds.salt, xgrid, hcoord='rho', scoord='w')
    """

    if attrs is None and isinstance(var, xr.DataArray):
        attrs = var.attrs.copy()
        attrs["name"] = var.name
        attrs["units"] = attrs.setdefault("units", "units")
        attrs["long_name"] = attrs.setdefault("long_name", "var")

    if hcoord is not None:
        assert hcoord in ["rho", "r", "psi", "p", "u", "v"], (
            'hcoord should be "rho" or "r" or "psi" or "p" or "u" or "v" but is "%s"'
            % hcoord
        )
        if hcoord in ["rho", "r"]:
            var = to_rho(var, xgrid, hboundary=hboundary, hfill_value=hfill_value)
        elif hcoord in ["psi", "p"]:
            var = to_psi(var, xgrid, hboundary=hboundary, hfill_value=hfill_value)
        elif hcoord == "u":
            var = to_u(var, xgrid, hboundary=hboundary, hfill_value=hfill_value)
        elif hcoord == "v":
            var = to_v(var, xgrid, hboundary=hboundary, hfill_value=hfill_value)

    if vcoord is not None:
        assert vcoord in ["s_rho", "rho", "r", "s_w", "w"], (
            'scoord should be "s_rho", "rho", "r", "s_w", or "w" but is "%s"' % scoord
        )
        if vcoord in ["s_rho", "rho", "r"]:
            var = to_s_rho(var, xgrid, vboundary=vboundary, vfill_value=vfill_value)
        elif vcoord in ["s_w", "w"]:
            var = to_s_w(var, xgrid, vboundary=vboundary, vfill_value=vfill_value)

    if isinstance(var, xr.DataArray):
        var.attrs = attrs
        var.name = var.attrs["name"]

    return var


def get_z(model, ds=None, z_sfc=None, h=None, xgrid=None, vgrid='r',
          hgrid='r', vtransform=2):
    """Compute vertical coordinates
    Spatial dimensions are placed last, in the order: s/s_w, y, x

    Args:
        ds: xarray dataset
        z_sfc: xarray.DataArray, optional
            Sea level data, default to 0 if not provided
            If you use slices, make sure singleton dimensions are kept, i.e do:
                z_sfc.isel(x=[i])
            and not :
                z_sfc.isel(x=i)
        h: xarray.DataArray, optional
            Water depth, searche depth in grid if not provided
        vgrid: str, optional
            Vertical grid, 'r'/'rho' or 'w'. Default is 'r'
        hgrid: str, optional
            Any horizontal grid: 'r'/'rho', 'u', 'v', 'f'. Default is 'r'
        vtransform: int, str, optional
            croco vertical transform employed in the simulation.
            1="old": z = z0 + (1+z0/_h) * _z_sfc  with  z0 = hc*sc + (_h-hc)*cs
            2="new": z = z0 * (_z_sfc + _h) + _z_sfc  with  z0 = (hc * sc + _h * cs) / (hc + _h)
    Return:
        DataArray : the z coordinate
    """
    xgrid = model.xgrid if xgrid is None else xgrid
    ds = model.ds if ds is None else ds

    h = ds.h if h is None else h
    z_sfc = 0*ds.h if z_sfc is None else z_sfc

    # switch horizontal grid if needed
    if hgrid in ['u','v','p']:
        h = to_grid_point(h, xgrid, hcoord=hgrid,vcoord=vgrid)
        z_sfc = to_grid_point(z_sfc, xgrid, hcoord=hgrid,vcoord=vgrid)

    # align datasets (z_sfc may contain a slice along one dimension for example)
    h, z_sfc  = xr.align(h, z_sfc, join='inner')

    if vgrid in ['r', 'rho']:
        vgrid = 'r'
        sc = ds['sc_r']
        cs = ds['Cs_r']
    else:
        sc = ds['sc_'+vgrid]
        cs = ds['Cs_'+vgrid]

    hc = ds['hc']

    if vtransform == 1:
        z0 = hc*sc + (h-hc)*cs
        z = z0 + (1+z0/h) * z_sfc
    elif vtransform == 2:
        z0 = (hc * sc + h * cs) / (hc + h)
        z = z0 * (z_sfc + h) + z_sfc

    # reorder spatial dimensions and place them last
    sdims = list(get_spatial_dims(z).values())
    sdims = tuple(filter(None,sdims)) # delete None values
    reordered_dims = tuple(d for d in z.dims if d not in sdims) + sdims
    z = z.transpose(*reordered_dims, transpose_coords=True).rename('z_'+hgrid)
    z.name = z.name.replace('z_r','z_'+vgrid)
    
    return z.fillna(0.) #.rename('z_'+hgrid).replace('z_r','z_'+vgrid)


def rot_uv(u, v, angle, xgrid):
    """Rotate u,v to lat,lon coordinates

    Args:
        u: DataArray
            3D velocity components in XI direction

        v: DataArray
            3D velocity components in ETA direction

        angle: DataArray
            Angle [radians] between XI-axis and the direction to the EAST at RHO-points

        xgrid: xgcm.grid
            grid object associated with u and v

    Returns:
        DatArray: rotated velocities, urot/vrot at the horizontal u/v grid point
    """

    assert xgrid is not None, "Xgcm grid should be input."
    assert isinstance(xgrid, Grid), "xgrid must be `xgcm` grid object."

    # save attributes to reinstitute at end
    attrsu = u.attrs
    attrsv = v.attrs

    # change all the variables to the rho grid point if needed
    u = to_grid_point(u, xgrid, hcoord="r", vcoord="r")
    v = to_grid_point(v, xgrid, hcoord="r", vcoord="r")
    angle = to_grid_point(angle, xgrid, hcoord="r", vcoord="r")

    cosang = np.cos(angle)
    sinang = np.sin(angle)

    # Rotate velocities
    urot = u * cosang - v * sinang
    vrot = u * sinang + v * cosang

    # change velocities to their horizontal grid point
    urot = to_grid_point(urot, xgrid, hcoord="u", vcoord="r")
    vrot = to_grid_point(vrot, xgrid, hcoord="v", vcoord="r")

    # add original attributes back in
    urot.attrs = {**attrsu, **urot.attrs}
    vrot.attrs = {**attrsv, **vrot.attrs}

    # try to guess cf coordinates
    urot = urot.squeeze().cf.guess_coord_axis()
    vrot = vrot.squeeze().cf.guess_coord_axis()

    return [urot.rename("urot"), vrot.rename("vrot")]



def hgrad(
    q,
    xgrid,
    which="both",
    z=None,
    hcoord=None,
    scoord=None,
    hboundary="extend",
    hfill_value=None,
    sboundary="extend",
    sfill_value=None,
    attrs=None,
):
    """Return gradients of property q accounting for s coordinates.

    Note that you need the 3D metrics for horizontal derivatives for ROMS, so ``include_3D_metrics=True`` in ``xroms.roms_dataset()``.

    Parameters
    ----------

    q: DataArray
        Property to take gradients of.
    xgrid: xgcm.grid
        Grid object associated with q.
    which: string, optional
        Which components of gradient to return.
        * 'both': return both components of hgrad.
        * 'xi': return only xi-direction.
        * 'eta': return only eta-direction.
    z: DataArray, ndarray, optional
        Depth [m]. If None, use z coordinate attached to q.
    hcoord: string, optional
        Name of horizontal grid to interpolate output to.
        Options are 'rho', 'psi', 'u', 'v'.
    scoord: string, optional
        Name of vertical grid to interpolate output to.
        Options are 's_rho', 's_w', 'rho', 'w'.
    hboundary: string, optional
        Passed to `grid` method calls; horizontal boundary selection
        for calculating horizontal derivatives of q. This same value
        will be used for all horizontal grid changes too.
        From xgcm documentation:
        A flag indicating how to handle boundaries:
        * None:  Do not apply any boundary conditions. Raise an error if
          boundary conditions are required for the operation.
        * 'fill':  Set values outside the array boundary to fill_value
          (i.e. a Neumann boundary condition.)
        * 'extend': Set values outside the array to the nearest array
          value. (i.e. a limited form of Dirichlet boundary condition.
    hfill_value: float, optional
        Passed to `grid` method calls; horizontal boundary selection
        fill value.
        From xgcm documentation:
        The value to use in the boundary condition with `boundary='fill'`.
    sboundary: string, optional
        Passed to `grid` method calls; vertical boundary selection
        for calculating horizontal derivatives of q. This same value will
        be used for all vertical grid changes too.
        From xgcm documentation:
        A flag indicating how to handle boundaries:
        * None:  Do not apply any boundary conditions. Raise an error if
          boundary conditions are required for the operation.
        * 'fill':  Set values outside the array boundary to fill_value
          (i.e. a Neumann boundary condition.)
        * 'extend': Set values outside the array to the nearest array
          value. (i.e. a limited form of Dirichlet boundary condition.
    sfill_value: float, optional
        Passed to `grid` method calls; vertical boundary selection
        fill value.
        From xgcm documentation:
        The value to use in the boundary condition with `boundary='fill'`.
    attrs: dict, optional
        Dictionary of attributes to add to resultant arrays. Requires that
        q is DataArray.

    Returns
    -------
    DataArray(s) of dqdxi and/or dqdeta, the gradients of q
    in the xi- and eta-directions with attributes altered to reflect calculation.

    Notes
    -----
    dqdxi = dqdx*dzdz - dqdz*dzdx

    dqdeta = dqdy*dzdz - dqdz*dzdy

    Derivatives are taken in the ROMS curvilinear grid native xi- and eta- directions.

    These derivatives properly account for the fact that ROMS vertical coordinates are
    s coordinates and therefore can vary in time and space.

    The xi derivative will alter the number of points in the xi and s dimensions.
    The eta derivative will alter the number of points in the eta and s dimensions.

    Examples
    --------
    >>> dtempdxi, dtempdeta = xroms.hgrad(ds.temp, xgrid)
    """

    assert isinstance(q, xr.DataArray), "var must be DataArray"

    if not [dim for dim in q.dims if dim.startswith('s')]:
        is3D = False
    else:
        is3D = True
        
    if is3D and z is None:
        try:
            coords = list(q.coords)
            z_coord_name = coords[[coord[:2] == "z_" for coord in coords].index(True)]
            z = q[z_coord_name]
            is3D = True
        except:
            # if we get here that means that q doesn't have z coords (like zeta)
            print("!!! Missing z coordinate, only horizontal gradient")
            is3D = False

    if which in ["both", "x"]:

        if is3D:
            dqdx = xgrid.interp(
                xgrid.derivative(q, "x", boundary=hboundary, fill_value=hfill_value),
                "z",
                boundary=sboundary,
                fill_value=sfill_value,
            )
            dqdz = xgrid.interp(
                xgrid.derivative(q, "z", boundary=sboundary, fill_value=sfill_value),
                "x",
                boundary=hboundary,
                fill_value=hfill_value,
            )
            dzdx = xgrid.interp(
                xgrid.derivative(z, "x", boundary=hboundary, fill_value=hfill_value),
                "z",
                boundary=sboundary,
                fill_value=sfill_value,
            )
            # dzdz = xgrid.interp(
            #     xgrid.derivative(z, "z", boundary=sboundary, fill_value=sfill_value),
            #     "x",
            #     boundary=hboundary,
            #     fill_value=hfill_value,
            # )

            # dqdx = dqdx * dzdz - dqdz * dzdx
            dqdx = dqdx - dqdz * dzdx

        else:  # 2D variables
            dqdx = xgrid.derivative(q, "x", boundary=hboundary, fill_value=hfill_value)

        if attrs is None and isinstance(q, xr.DataArray):
            attrs = q.attrs.copy()
            attrs["name"] = "d" + q.name + "dx"
            attrs["units"] = "1/m * " + attrs.setdefault("units", "units")
            attrs["long_name"] = "horizontal xi derivative of " + attrs.setdefault(
                "long_name", "var"
            )

    if which in ["both", "y"]:

        if is3D:
            dqdy = xgrid.interp(
                xgrid.derivative(q, "y", boundary=hboundary, fill_value=hfill_value),
                "z",
                boundary=sboundary,
                fill_value=sfill_value,
            )
            dqdz = xgrid.interp(
                xgrid.derivative(q, "z", boundary=sboundary, fill_value=sfill_value),
                "y",
                boundary=hboundary,
                fill_value=hfill_value,
            )
            dzdy = xgrid.interp(
                xgrid.derivative(z, "y", boundary=hboundary, fill_value=hfill_value),
                "z",
                boundary=sboundary,
                fill_value=sfill_value,
            )
            dzdz = xgrid.interp(
                xgrid.derivative(z, "z", boundary=sboundary, fill_value=sfill_value),
                "y",
                boundary=hboundary,
                fill_value=hfill_value,
            )

            dqdy = dqdy * dzdz - dqdz * dzdy

        else:  # 2D variables
            dqdy = xgrid.derivative(
                q, "y", boundary=hboundary, fill_value=hfill_value
            )

        if attrs is None and isinstance(q, xr.DataArray):
            attrs = q.attrs.copy()
            attrs["name"] = "d" + q.name + "dy"
            attrs["units"] = "1/m * " + attrs.setdefault("units", "units")
            attrs["long_name"] = "horizontal eta derivative of " + attrs.setdefault(
                "long_name", "var"
            )

    if which == "both":
        return dqdx, dqdy
    elif which == "x":
        return dqdx
    elif which == "y":
        return dqdy
    else:
        print("nothing being returned from hgrad")
        

def get_grid_point(var):
    """Get the horizontal and vertical grid point of a variable

    Args:
        var (DataArray): variable to operate on

    Returns:
        character, character: horizontal, vertical grid point
    """
    dims = var.dims
    # horizontal point
    hpoint='r'
    if "x_u" in dims:
        if "y" in dims:
            hpoint='u'
        else:
            hpoint='f'
    elif "y_v" in dims:
        hpoint='v'
    if 's' in dims:
        vpoint='r'
    else:
        vpoint='w'
    return hpoint,vpoint
        
    
def slices(model, var, z, ds=None, xgrid=None, longitude=None, latitude=None, depth=None):
    """
    This function interpolate a 3D variable on slices at constant depths/longitude/latitude
    This function use xcgm transform method and needs xgcm.Grid to be defined over the 3 axes.
    !!! For now, it works only with curvilinear coordinates !!!

    Args:
        ds      dataset to find the grid
        var     (dataArray) Variable to process (3D matrix).
        z       (dataArray) Depths at the same point than var (3D matrix).
        longitude   (scalar,list or ndarray) longitude of the slice (scalar meters, negative).
        latitude    (scalar,list or ndarray) latitude of the slice (scalar meters, negative).
        depth       (scalar,list or ndarray) depth of the slice (scalar meters, negative).
    Return:
        vnew    (dataArray) Horizontal slice
    """
    from matplotlib.cbook import flatten
  
    xgrid = model.xgrid if xgrid is None else xgrid
    ds = model.ds if ds is None else ds
    
    if longitude is None and latitude is None and depth is None:
        "Longitude or latitude or depth must be defined"
        return None

    # check typ of longitude/latitude/depth
    # longitude = longitude.tolist() if isinstance(longitude,np.ndarray) else longitude
    # longitude = [longitude] if (isinstance(longitude,int) or isinstance(longitude,float)) else longitude
    longitude = np.asarray(longitude) if isinstance(longitude,list) else longitude
    longitude = np.asarray([longitude]) if (isinstance(longitude,int) or isinstance(longitude,float)) else longitude

    # latitude = latitude.tolist() if isinstance(latitude,np.ndarray) else latitude
    # latitude = [latitude] if (isinstance(latitude,int) or isinstance(latitude,float)) else latitude
    latitude = np.asarray(latitude) if isinstance(latitude,list) else latitude
    latitude = np.asarray([latitude]) if (isinstance(latitude,int) or isinstance(latitude,float)) else latitude

    # depth = depth.tolist() if isinstance(depth,np.ndarray) else depth
    # depth = [depth] if (isinstance(depth,int) or isinstance(depth,float)) else depth
    depth = np.asarray(depth) if isinstance(depth,list) else depth
    depth = np.asarray([depth]) if (isinstance(depth,int) or isinstance(depth,float)) else depth

     # Find dimensions and coordinates of the variable
    dims = get_spatial_dims(var)
    coords = get_spatial_coords(var)
    if dims['s'] is not None and coords['z'] is None: 
        var = var.assign_coords(coords={'z':z})
        coords = get_spatial_coords(var)
    # hgrid,vgrid = get_grid_point(var)

    if longitude is not None:
        axe = 'x'
        coord_ref = coords['lon']
        coord_x = coords['lat']
        coord_y = coords['z']
        slices_values = longitude
    elif latitude is not None:
        axe = 'y'
        coord_ref = coords['lat']
        coord_x = coords['lon']
        coord_y = coords['z']
        slices_values = latitude
    else:
        axe = 'z'
        coord_ref = coords['z']
        coord_x = coords['lon']
        coord_y = coords['lat']
        slices_values = depth

    # Recursively loop over time if needed
    if len(var.squeeze().dims) == 4:
        vnew = [slices(model, var.isel(t=t), z.isel(t=t), ds=ds, xgrid=xgrid,
                      longitude=longitude, latitude=latitude, depth=depth)
                      for t in range(len(var.t))]
        vnew = xr.concat(vnew, dim='t')
    else:
        vnew = xgrid.transform(var, axe, slices_values,
                               target_data=var[coord_ref]).squeeze()
    # Do the linear interpolation
    if not depth:
        x = xgrid.transform(var[coord_x], axe, slices_values,
                                   target_data=var[coord_ref]).squeeze() #\
                     #.expand_dims({dims['s']: len(var[dims['s']])})            
        vnew = vnew.assign_coords(coords={coord_x:x})

        #y = xgrid.transform(var[coord_y], axe, slices_values,
        if dims['s'] is not None:
            y = xgrid.transform(z, axe, slices_values,
                                       target_data=var[coord_ref]).squeeze()
            # Add the coordinates to dataArray
            vnew = vnew.assign_coords(coords={coord_y:y})
    else:
        # Add the coordinates to dataArray
        vnew = vnew.assign_coords(coords={coord_x:var[coord_x]})
        vnew = vnew.assign_coords(coords={coord_y:var[coord_y]})

#     return vnew.squeeze().unify_chunks().fillna(0.)  #unify_chunks() 
    return vnew.squeeze().fillna(0.)  #unify_chunks()   


def isoslice(var, target, xgrid, target_data=None, axis="z"):
    """Interpolate var to target.

    This wraps `xgcm` `transform` function for slice interpolation,
    though `transform` has additional functionality.

    Args:
        var: DataArray
            Variable to operate on.
        target: ndarray
            Values to interpolate to. If calculating var at fixed depths,
            target are the fixed depths, which should be negative if
            below mean sea level. If input as array, should be 1D.
        xgrid: xgcm.grid, optional
            Grid object associated with var.
        target_data: DataArray, optional
            Array that var is interpolated onto (e.g., z coordinates or
            density). If calculating var on fixed depth slices, target_data
            contains the depths [m] associated with var. In that case and
            if None, will use z coordinate attached to var. Also use this
            option if you want to interpolate with z depths constant in
            time and input the appropriate z coordinate.
        axis: str, optional
            Dimension over which to calculate isoslice. If calculating var
            onto fixed depths, `dim='z'`. Options are 'z', 'y', and 'x'.
    Return:
        DataArray of var interpolated to target. Dimensionality will be the
        same as var except with dim dimension of size of target.

    Notes
    -----
    var cannot have chunks in the dimension dim.

    cf-xarray should still be usable after calling this function.

    Examples
    --------
    To calculate temperature onto fixed depths:
    >>> isoslice(ds.temp, np.linspace(0, -30, 50))

    To calculate temperature onto salinity:
    >>> isoslice(ds.temp, np.arange(0, 36), target_data=ds.salt, axis='z')

    Calculate lat-z slice of salinity along a constant longitude value (-91.5):
    >>> isoslice(ds.salt, -91.5, target_data=ds.lon_rho, axis='x')

    Calculate slice of salt at 28 deg latitude
    >>> isoslice(ds.salt, 28, target_data=ds.lat_rho, axis='y')

    Interpolate temp to salinity values between 0 and 36 in the X direction
    >>> isoslice(ds.temp, np.linspace(0, 36, 50), target_data=ds.salt, axis='x')

    Interpolate temp to salinity values between 0 and 36 in the Z direction
    >>> isoslice(ds.temp, np.linspace(0, 36, 50), target_data=ds.salt, axis='z')

    Calculate the depth of a specific isohaline (33):
    >>> isoslice(ds.salt, 33, target_data=ds.z_rho, axis='z')


    """

    assert xgrid is not None, "Xgcm grid should be input."

    assert isinstance(xgrid, Grid), "xgrid must be `xgcm` grid object."

    attrs = var.attrs  # save to reinstitute at end

    # make sure target are array-like
    if isinstance(target, (int, float)):
        target = np.asarray([target])

    # interpolate to the z coordinates associated with var
    if target_data is None:
        key = [coord for coord in var.coords if "z_" in coord][0]
        assert (
            len(key) > 0
        ), "z coordinates associated with var could not be identified."
        target_data = var[key]
    else:
        if isinstance(target_data, xr.DataArray) and target_data.name is not None:
            key = target_data.name
        else:
            key = "z"

    # perform interpolation
    transformed = xgrid.transform(var, axis, target, target_data=target_data)

    if key not in transformed.coords:
        transformed = transformed.assign_coords({key: target_data})

    # bring along attributes for cf-xarray
    transformed[key].attrs["axis"] = axis
    # add original attributes back in
    transformed.attrs = {**attrs, **transformed.attrs}

    # save key names for later
    # perform interpolation for other coordinates if needed
    if "longitude" in var.cf.standard_names:
        # lonkey = var.cf["longitude"].name
        lonkey = var.cf.standard_names["longitude"][0]

        if lonkey not in transformed.coords:
            # this interpolation won't work for certain combinations of var[latkey] and target_data
            # without the following step
            if "T" in target_data.reset_coords(drop=True).cf.axes:
                target_data = target_data.cf.isel(T=0).drop_vars(
                    target_data.cf["T"].name, errors="ignore"
                )
            if "Z" in target_data.reset_coords(drop=True).cf.axes:
                target_data = target_data.cf.isel(Z=0).drop_vars(
                    target_data.cf["Z"].name, errors="ignore"
                )
            transformedlon = xgrid.transform(
                var[lonkey], axis, target, target_data=target_data
            )
            transformed = transformed.assign_coords({lonkey: transformedlon})

        transformed[lonkey].attrs["standard_name"] = "longitude"

    if "latitude" in var.cf.standard_names:
        # latkey = var.cf["latitude"].name
        latkey = var.cf.standard_names["latitude"][0]

        if latkey not in transformed.coords:
            # this interpolation won't work for certain combinations of var[latkey] and target_data
            # without the following step
            if "T" in target_data.reset_coords(drop=True).cf.axes:
                target_data = target_data.cf.isel(T=0).drop_vars(
                    target_data.cf["T"].name, errors="ignore"
                )
            if "Z" in target_data.reset_coords(drop=True).cf.axes:
                target_data = target_data.cf.isel(Z=0).drop_vars(
                    target_data.cf["Z"].name, errors="ignore"
                )
            transformedlat = xgrid.transform(
                var[latkey], axis, target, target_data=target_data
            )
            transformed = transformed.assign_coords({latkey: transformedlat})

        transformed[latkey].attrs["standard_name"] = "latitude"

    if "vertical" in var.cf.standard_names:
        # zkey = var.cf["vertical"].name
        zkey = var.cf.standard_names["vertical"][0]

        if zkey not in transformed.coords:
            transformedZ = xgrid.transform(
                var[zkey], axis, target, target_data=target_data
            )
            transformed = transformed.assign_coords({zkey: transformedZ})

        transformed[zkey].attrs["positive"] = "up"

    transformed = transformed.squeeze().cf.guess_coord_axis()

    # reorder back to normal ordering in case changed
    transformed = order_dims(transformed)

    return transformed



def cross_section(grid, da, lon1, lat1, lon2, lat2, dlon=None):
    """ Extract a section between 2 geographic points

    Args:
        grid (XGCM grid): the XGCM grid associated
        da (DataArray): variable to operate on
        lon1 (float): minimum longitude
        lat1 (float): minimum latitude
        lon2 (float): maximum longitude
        lat2 (float): maximum latitude
        dlon (float, optional): longitude interval. Defaults to None.

    Returns:
        DataArray: new section
    """
    
    # check input parameters
    if not isinstance(grid,Grid): print('grid must be a xgcm grid'); return None
    if not isinstance(da,xr.DataArray): print('da must be a xarray DataArray'); return None
    dims = get_spatial_dims(da)
    coords = get_spatial_coords(da)
    if coords['lon'] is None or coords['lat'] is None:
        print('da must have longitude AND latitude coordinates')
        return None
    if not isinstance(lon1,(int,float)): print('lon1 must be a float'); return None
    if not isinstance(lat1,(int,float)): print('lat1 must be a float'); return None
    if not isinstance(lon2,(int,float)): print('lon2 must be a float'); return None
    if not isinstance(lat2,(int,float)): print('lat2 must be a float'); return None
    if dlon is not None and not isinstance(dlon,(int,float)): print('dlon must be a number'); return None
           
    # compute the linear function from the two points
    a = (lat2 - lat1)/(lon2 - lon1)
    b = lat1 - a*lon1
    
    # get the longitude interval of the new grid , compute the new longitude grid
    if dlon is None:
        dlon = ((da[coords['lon']].max().values - da[coords['lon']].min().values) /
                da[coords['lon']].sizes[dims['x']])
    longrd = np.arange(lon1,lon2,dlon)
    
    # compute the latitude coordinates of the new grid with the linear function
    latgrd = a * longrd + b

    # interpolate on the regular longitude grid
    newda = auto_chunk(da, keep_complete_dim='x', wanted_chunk=200)
    newda = grid.transform(newda,'x', longrd, target_data=da[coords['lon']])
    newda = newda.rename({coords['lon']:dims['x']})
    newlat = auto_chunk(da[coords['lat']], keep_complete_dim='x', wanted_chunk=200)
    newlat = grid.transform(newlat, 'x', longrd, target_data=da[coords['lon']])              
    newlat = newlat.rename({coords['lon']:dims['x']})
    if coords['z'] is not None:
        newz = auto_chunk(da[coords['z']], keep_complete_dim='x', wanted_chunk=200)
        newz = grid.transform(newz,'x', longrd, target_data=da[coords['lon']]).fillna(0.)
        newz = newz.rename({coords['lon']:dims['x']})
        
    # interpolate on a new latitude grid
    newda = auto_chunk(newda, keep_complete_dim='y', wanted_chunk=200)
    newda = grid.transform(newda,'y',latgrd,target_data=newlat)
    newda = newda.rename({coords['lat']:dims['y']})
    if coords['z'] is not None:
        newz = auto_chunk(newz, keep_complete_dim='y', wanted_chunk=200)
        newz = grid.transform(newz,'y',latgrd,target_data=newlat).fillna(0.)
        newz = newz.rename({coords['lat']:dims['y']}).fillna(0.)
    
    # extract the cross section
    crossda = []; crossz=[]
    for lon,lat in zip(longrd,latgrd):
        crossda.append(newda.loc[{dims['x']:lon}].loc[{dims['y']:lat}])
        if coords['z'] is not None:
            crossz.append(newz.loc[{dims['x']:lon}].loc[{dims['y']:lat}])

    cross = xr.concat(crossda,dim=dims['x'])
    if coords['z'] is not None:
        crossz = xr.concat(crossz,dim=dims['x'])
    
    # assign the coordinates lon/lat/z to the section
    cross = cross.assign_coords({coords['lon']:cross[dims['x']],
                                 coords['lat']:cross[dims['y']]})
    if coords['z'] is not None:
        cross = cross.assign_coords({coords['z']:crossz})    
        
    return reorder_dims(cross.fillna(0.))

def interp_regular(da,grid,axis,tgrid,rgrid=None):
    """
    interpolate on a regular grid
    Args:
        - da (DataArray) : variable to interpolate
        - grid (xgcm grid): xgcm grid
        - axis (str): axis of the xgcm grid for the interpolation ('x', 'y' or 'z')
        - tgrid (numpy vector): target relular grid space
        - rgrid (numpy array or DataArray): reference grid of da
    Return:
        - (DataArray): regurlarly interpolated variable
    """
    
    # check axis
    if axis not in ['x','y','z']: 
        print('axis must be x, y or z')
        return None
    
    # corresponding keys between spatial coords/dims and axes of the xgcm grid
    refc = {'x':'lon', 'y':'lat', 'z':'z'}
    refd = {'x':'x', 'y':'y', 'z':'s'}
    
    # find spatial coordinates/dims of da
    coords = get_spatial_coords(da)
    coord = coords[refc[axis]]
    dims = get_spatial_dims(da)
    dim = dims[refd[axis]]
    
    # initialize the reference coordinate
    if rgrid is None and coord is not None:
        rgrid = da[coord]
    else:
        print('the reference grid is missing along the axis of interpolation')
        return None

    # interpolate da on the regular grid
    newvar = grid.transform(da,axis,tgrid,target_data=rgrid).rename({coord:dim})
    
    return reorder_dims(newvar )


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    Args:
        lon1 (float): minimum longitude
        lat1 (float): minimum latitude
        lon2 (float): maximum longitude
        lat2 (float): maximum latitude
    Returns:
        float: distance in km
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    # Radius of earth in kilometers is 6371
    km = 6371* c
    return km

    
# ----------------------------- grid rechunk -----------------------------

def auto_chunk(ds, keep_complete_dim=None, wanted_chunk=150):
    """
    Rechunk Dataset or DataArray such as each partition size is about 150Mb
    Args:
        - ds : (Dataset or DataArray) object to rechunk
        - keep_complete_dim : (character) Horizontal axe to keep with no chunk ('x','y','s')
        - wanted_chunk : (integer) size of each partition in Mb
    Return:
        - object rechunked
    """

    #check input parameters
    if not isinstance(ds, (xr.Dataset,xr.DataArray)):
        print('argument must be a xarray.DataArray or xarray.Dataset')
        return
    if keep_complete_dim is not None:
        keep_complete_dim = list(keep_complete_dim) if not isinstance(keep_complete_dim,list) else keep_complete_dim
        if not all(item in ['s','y','x'] for item in keep_complete_dim):
            print('keep_complete_dim must equal x or y or s')
            return

    # get horizontal dimensions names of the Dataset/DataArray
    dname = get_spatial_dims(ds)
    # remove None values
    dname = {k: v for k, v in dname.items() if v is not None}
    chunks_name = dname.copy()
    
    # get max dimensions sizes of the Dataset/DataArray
    chunks_size={}
    for k,v in chunks_name.items():
        if isinstance(v,list):       # for a dataset
            chunks_size[k] = max([ds.sizes[d] for d in v])
        else:
            chunks_size[k] = ds.sizes[v]

    # always chunk in time
    if 't' in ds.dims: chunks_size['t'] = 1
        
    if keep_complete_dim:
        # remove keep_complete_dim from the dimensions of the Dataset/DatAarray
        for d in keep_complete_dim: del chunks_name[d]
        
    # reduce chunks size  beginning by 's' then 'y' then 'x' if necessary
    for k in chunks_name.keys():
        for d in range(chunks_size[k],0,-1):
            # chunk_size = (chunks_size['x']*chunks_size['y']*chunks_size['s']*4 / 1.e6)
            chunk_size = 4 / 1.e6
            for chunk in chunks_size.values():
                chunk_size = chunk_size*chunk
            if chunk_size > wanted_chunk:
                chunks_size[k] = d
            else:
                break
        if chunk_size > wanted_chunk:
            break

    if isinstance(ds,xr.Dataset):
        # set chunk for all the dimensions of the dataset (ie : x and x_u)
        for c in list(itertools.product(dname.keys(), ds.dims.keys())):
            if c[1].startswith(c[0]):
                chunks_size[c[1]] = chunks_size[c[0]]
    else:
        # rename the dimension name by the right values (ie: x_u instead of x)
        for key in dname.keys():
            chunks_size[dname[key]] = chunks_size.pop(key)

        
    return ds.chunk(chunks_size)

