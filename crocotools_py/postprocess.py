import numpy as np
import xarray as xr
import dask
from datetime import timedelta, datetime
from glob import glob
from crocotools_py.define_attrs import CROCO_Attrs_RotatedVectors, CROCO_Attrs
import pandas as pd
import os
import shutil
import time
from dask.diagnostics import ProgressBar


def change_attrs(attrs,da,var_str):
    meta = getattr(attrs, var_str)
    da.attrs['long_name'] = meta.long_name
    da.attrs['units'] = meta.units
    da.attrs['standard_name'] = meta.standard_name
    
    return da

def u2rho(u):
    """
    regrid the croco u-velocity from it's native u grid to the rho grid
    u can be 2D, 3D or 4D numpy array or an xarray dataarray
    returns a numpy array (or xarray dataarray, depending on the input) of the data on the rho grid
    """
    Num_dims=len(u.shape)
    if Num_dims==4:
        # T: Num Time 
        # D: Num Depth
        # Mp: Num eta_rho
        # L:  Num xi_u
        [T, D, Mp, L] = u.shape
        
        # extend the xi axis by one on either side before doing the interpolating along the xi axis
        u_extended = np.zeros((T, D, Mp, L + 2))
        u_extended[:,:,:,1:-1] = u
        u_extended[:,:,:,0] = u[:,:,:,0] # repeat adjacent values along the extended part
        u_extended[:,:,:,-1] = u[:,:,:,-1] # repeat adjacent values along the extended part
        
        # Interpolate from u grid to rho grid by summing adjacent u points, and divide by 2
        # output is on the rho grid due to extending the xi axis
        u_rho = 0.5 * (u_extended[:, :, :, 0 : -1] + u_extended[:, :, :, 1 :])
        
    elif Num_dims==3:
        # TorD: Time or Depth
        [TorD, Mp, L] = u.shape # works if first dimension is time or depth
        
        # extend the xi axis by one on either side before doing the interpolating along the xi axis
        u_extended = np.zeros((TorD, Mp, L + 2))
        u_extended[:,:,1:-1] = u
        u_extended[:,:,0] = u[:,:,0] # repeat adjacent values along the extended part
        u_extended[:,:,-1] = u[:,:,-1] # repeat adjacent values along the extended part
        
        # Interpolate from u grid to rho grid by summing adjacent u points, and divide by 2
        # output is on the rho grid due to extending the xi axis
        u_rho = 0.5 * (u_extended[:, :, 0 : -1] + u_extended[:, :, 1 :])
        
    else: # Num_dims==2:
        # 2D grid - no time/depth information (possibly a surface level at a single time step)
        [Mp, L] = u.shape
        
        # extend the xi axis by one on either side before doing the interpolating along the xi axis
        u_extended = np.zeros((Mp, L + 2))
        u_extended[:,1:-1] = u
        u_extended[:,0] = u[:,0] # repeat adjacent values along the extended part
        u_extended[:,-1] = u[:,-1] # repeat adjacent values along the extended part
        
        # Interpolate from u grid to rho grid by summing adjacent u points, and divide by 2
        # output is on the rho grid due to extending the xi axis
        u_rho = 0.5 * (u_extended[:, 0 : -1] + u_extended[:, 1 :])
        
    return u_rho

def v2rho(v):
    """
    regrid the croco v-velocity from it's native v grid to the rho grid
    v can be 2D, 3D or 4D numpy array or an xarray dataarray
    returns a numpy array (or xarray dataarray, depending on the input) of the data on the rho grid
    """
    Num_dims=len(v.shape)
    if Num_dims==4:
        # T: Num Time 
        # D: Num Depth
        # M: Num eta_v
        # Lp:  Num xi_rho
        [T, D, M, Lp] = v.shape
        
        # extend the eta axis by one on either side before doing the interpolating along the eta axis
        v_extended = np.zeros((T, D, M + 2, Lp))
        v_extended[:,:,1:-1,:] = v
        v_extended[:,:,0,:] = v[:,:,0,:] # repeat adjacent values along the extended part
        v_extended[:,:,-1,:] = v[:,:,-1,:] # repeat adjacent values along the extended part
        
        # Interpolate from v grid to rho grid by summing adjacent v points, and divide by 2
        # output is on the rho grid due to extending the eta axis
        v_rho = 0.5 * (v_extended[:, :, 0 : -1, :] + v_extended[:, :, 1 :, :])
        
        
    elif Num_dims==3:
        [TorD, M, Lp] = v.shape # works if first dimension is time or depth
        
        # extend the eta axis by one on either side before doing the interpolating along the eta axis
        v_extended = np.zeros((TorD, M + 2, Lp))
        v_extended[:,1:-1,:] = v
        v_extended[:,0,:] = v[:,0,:] # repeat adjacent values along the extended part
        v_extended[:,-1,:] = v[:,-1,:] # repeat adjacent values along the extended part
        
        # Interpolate from v grid to rho grid by summing adjacent v points, and divide by 2
        # output is on the rho grid due to extending the eta axis
        v_rho = 0.5 * (v_extended[:, 0 : -1, :] + v_extended[:, 1 :, :])
        
    else: # Num_dims==2:
        [M, Lp] = v.shape
        
        # extend the eta axis by one on either side before doing the interpolating along the eta axis
        v_extended = np.zeros((M + 2, Lp))
        v_extended[1:-1,:] = v
        v_extended[0,:] = v[0,:] # repeat adjacent values along the extended part
        v_extended[-1,:] = v[-1,:] # repeat adjacent values along the extended part
        
        # Interpolate from v grid to rho grid by summing adjacent v points, and divide by 2
        # output is on the rho grid due to extending the eta axis
        v_rho = 0.5 * (v_extended[0 : -1, :] + v_extended[1 :, :])
        
    return v_rho

def psi2rho(var_psi):
    """
    regrid a variable on the psi grid to the rho grid
    """
     
    Num_dims=len(var_psi.shape)
    if Num_dims==2:
    
        [M,L]=var_psi.shape

        var_rho=np.zeros((M+1,L+1))
        
        var_rho[1:M, 1:L] = 0.25 * (var_psi[:-1, 1:] + var_psi[1:, 1:] + var_psi[:-1, :-1] + var_psi[1:, :-1])
        
        var_rho[:, 0] = var_rho[:, 1]
        var_rho[:, L] = var_rho[:, L - 1]
        var_rho[0, :] = var_rho[1, :]
        var_rho[M, :] = var_rho[M - 1,:]
        
    elif Num_dims==3:
    
        [T_D,M,L]=var_psi.shape

        var_rho=np.zeros((T_D,M+1,L+1))
        
        var_rho[:, 1:M, 1:L] = 0.25 * (var_psi[:, :-1, 1:] + var_psi[:, 1:, 1:] + var_psi[:, :-1, :-1] + var_psi[:, 1:, :-1])
        
        var_rho[:, :, 0] = var_rho[:, :, 1]
        var_rho[:, :, L] = var_rho[:, :, L - 1]
        var_rho[:, 0, :] = var_rho[:, 1, :]
        var_rho[:, M, :] = var_rho[:, M - 1,:]
    
    else: # Num_dims==4:
    
        [T,D,M,L]=var_psi.shape

        var_rho=np.zeros((T,D,M+1,L+1))
        
        var_rho[:, :, 1:M, 1:L] = 0.25 * (var_psi[:, :, :-1, 1:] + var_psi[:, :, 1:, 1:] + var_psi[:, :, :-1, :-1] + var_psi[:, :, 1:, :-1])
        
        var_rho[:, :, :, 0] = var_rho[:, :, :, 1]
        var_rho[:, :, :, L] = var_rho[:, :, :, L - 1]
        var_rho[:, :, 0, :] = var_rho[:, :, 1, :]
        var_rho[:, :, M, :] = var_rho[:, :, M - 1,:]
    
    return var_rho

def rho2u(var_rho):
    """
    regrid a variable on the rho grid to the u-grid
    """
    Num_dims=len(var_rho.shape)
    if Num_dims == 2:
        [Mp,Lp]=var_rho.shape
        L=Lp-1
        var_u=0.5*(var_rho[:,0:L]+var_rho[:,1:Lp])
    elif Num_dims == 3:
        [T_D,Mp,Lp]=var_rho.shape
        L=Lp-1
        var_u=0.5*(var_rho[:,:,0:L]+var_rho[:,:,1:Lp])
    else: # Num_dims == 4:
        [T,D,Mp,Lp]=var_rho.shape
        L=Lp-1
        var_u=0.5*(var_rho[:,:,:,0:L]+var_rho[:,:,:,1:Lp])
    return var_u

def rho2v(var_rho):
    """
    regrid a variable on the rho grid to the v-grid
    """
    Num_dims=len(var_rho.shape)
    if Num_dims == 2:
        [Mp,Lp]=var_rho.shape
        M=Mp-1
        var_v=0.5*(var_rho[0:M,:]+var_rho[1:Mp,:])
    elif Num_dims == 3:
        [T_D,Mp,Lp]=var_rho.shape
        M=Mp-1
        var_v=0.5*(var_rho[:,0:M,:]+var_rho[:,1:Mp,:])
    else: # Num_dims == 4:
        [T,D,Mp,Lp]=var_rho.shape
        M=Mp-1
        var_v=0.5*(var_rho[:,:,0:M,:]+var_rho[:,:,1:Mp,:])
    return var_v

def csf(sc, theta_s, theta_b):
    """
    Allows use of theta_b > 0 (July 2009)
    is required in zlevs.py
    """
    one64 = np.float64(1)

    if theta_s > 0.0:
        csrf = (one64 - np.cosh(theta_s * sc)) / (np.cosh(theta_s) - one64)
    else:
        csrf = -(sc**2)
    sc1 = csrf + one64
    if theta_b > 0.0:
        Cs = (np.exp(theta_b * sc1) - one64) / (np.exp(theta_b) - one64) - one64
    else:
        Cs = csrf

    return Cs

def z_levels(h, zeta, theta_s, theta_b, hc, N, type, vtransform):
    """
    Vectorised version of z_levels, accepts zeta with shape (T, M, L) and returns
    z with shape (T, N, M, L).
    """

    """
    this provides a 3D grid of the depths of the sigma levels
    h = 2D bathymetry of your grid
    zeta = zeta at particular timestep or also accepts zeta with shape (T, M, L) 
    theta_s = surface stretching parameter
    theta_b = bottom stretching parameter
    hc = critical depth
    N = number of sigma levels
    type = 'w' or 'rho'
    vtransform = 1 (OLD) or 2 (NEW)
    
    returns
    z with shape (T, N, M, L).

    this is adapted (by J.Veitch and G.Fearon - Feb 2022, Mar 2025) from zlevs.m in roms_tools (by P. Penven)
    """

    if zeta.ndim == 2:
        # Add time dimension if missing
        zeta = zeta[None, :, :]

    T, M, L = zeta.shape

    # Stretching coordinates
    if vtransform == 2:
        ds = 1 / N
        if type == "w":
            sc = np.linspace(-1.0, 0.0, N + 1)
            Cs = csf(sc, theta_s, theta_b)
            N += 1
        else:
            sc = ds * (np.arange(1, N + 1) - N - 0.5)
            Cs = csf(sc, theta_s, theta_b)
    else:
        if type == "w":
            sc = (np.arange(0, N + 1) - N) / N
            N += 1
        else:
            sc = (np.arange(1, N + 1) - N - 0.5) / N
        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5 * theta_s)
        Cs = (1 - theta_b) * cff1 * np.sinh(theta_s * sc) + theta_b * (
            cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5
        )

    # Safety check
    h = np.where(h == 0, 1e-2, h)  # avoid divide-by-zero
    zeta = np.maximum(zeta, 0.01 - h[None, :, :])  # ensure zeta >= Dcrit - h

    # Broadcast h to match zeta
    h_bcast = h[None, :, :]  # shape (1, M, L)

    if vtransform == 2:
        h2 = h_bcast + hc
        h2inv = 1.0 / h2
        cff = hc * sc[:, None, None] + Cs[:, None, None] * h_bcast

        z = cff * h_bcast / h2 + zeta[:, None, :, :] * (1.0 + cff * h2inv)
    else:
        hinv = 1.0 / h_bcast
        cff = hc * (sc[:, None, None] - Cs[:, None, None])
        z = cff + Cs[:, None, None] * h_bcast + zeta[:, None, :, :] * (
            1.0 + (cff + Cs[:, None, None] * h_bcast) * hinv
        )

    return z  # shape (T, N, M, L)

def hlev_xarray(var, z, depth):
    """
    This function interpolates a 4D xarray DataArray on horizontal levels of constant depth(s).

    INPUT:
        var     Variable to process (xarray DataArray: time, s_rho, eta_rho, xi_rho).
        z       Depths (m) of RHO- or W-points (xarray DataArray: time, s_rho, eta_rho, xi_rho).
        depth   Slice depth(s) (scalar or list of scalars; meters, negative).

    OUTPUT:
        vnew    Horizontal slice(s) (xarray DataArray: time, depth, eta_rho, xi_rho. if depth is a list, otherwise time, eta_rho, xi_rho).
    """
    # this function can be slow if the variable was read using dask (i.e. using open_mfdataset)
    # so we load the data into memory first
    var=var.compute()
    z=z.compute()
    
    # Convert depth to xarray DataArray if it's a scalar or a list
    if np.isscalar(depth):
        depth = xr.DataArray([depth], dims="depth")
    else:
        depth = xr.DataArray(depth, dims="depth")
        
    # Add attributes to the depth coordinate
    depth.attrs["long_name"] = "Depth"
    depth.attrs["units"] = "m"
    depth.attrs["standard_name"] = "depth"
    depth.attrs["positive"] = "up"
    
    # Determine the nearest vertical levels where z brackets each depth
    below_depth = z < depth
    levs = below_depth.sum(dim="s_rho")  # Find levels below depth for each case
    levs = levs.clip(1, z.sizes["s_rho"] - 1)  # Ensure valid indices

    # Mask where the depth is outside the range of z
    # Outside range: depth < z_bottom or depth > z_top
    z_bottom = z.isel(s_rho=0)
    mask = xr.where((depth < z_bottom), np.nan, 1)

    # Extract the bracketing levels
    z_up = z.isel(s_rho=levs)
    z_down = z.isel(s_rho=levs - 1)
    v_up = var.isel(s_rho=levs)
    v_down = var.isel(s_rho=levs - 1)

    # Linear interpolation along the depth dimension
    vnew = mask * (
        ((v_up - v_down) * depth + v_down * z_up - v_up * z_down) / (z_up - z_down)
    )

    # Assign the depth coordinate and reorder dimensions
    vnew = vnew.assign_coords(depth=depth)
    vnew = vnew.transpose("time", "depth", "eta_rho", "xi_rho")
    
    # Drop the s_rho coordinate if it exists
    vnew = vnew.drop_vars("s_rho", errors="ignore")

    return vnew.squeeze("depth") if len(depth) == 1 else vnew

def hlev(var,z,depth):
    """
    This function interpolate a 3D variable on a horizontal level of constant depth
    
    INPUT:
        var     Variable to process (3D matrix).
        z       Depths (m) of RHO- or W-points (3D matrix).
        depth   Slice depth (scalar; meters, negative).

    OUTPUT: 
        vnew    Horizontal slice (2D matrix). 
    
    """

    # Identify the grid dimensions
    # N (number of vertical levels), M (number of rows), L (number of columns) 
    [N,M,L]=np.shape(z)
    i1=np.arange(0,L)
    j1=np.arange(0,M)
    
    # Determine nearest vertical levels
    # Logic: For each horizontal position (row j and column i), identify which vertical levels (z) bracket the target depth
    # adjust the levels at the boundaries to avoid indexing errors.
    a=np.int_(z<depth)
    levs=np.sum(a,axis=0)
    levs[np.where(levs==N)] = N-1
    
    # Handle Masking
    # Creates a mask to ignore invalid values (i.e. locations where the target depth is outside the valid range of z).
    mask=0.*levs + 1.
    mask[np.where(levs==0)]=np.nan
    
    # Locate the bracketing levels
    # Converts the 3D indices into linear indices for easier slicing of z and var (avoids having to loop over each index point).
    [i2,j2]=np.meshgrid(i1,j1)
    pos_up  = L*M*levs + L*j2 + i2
    pos_down= L*M*(levs-1) + L*j2 + i2
    
    # Extract the Bracketing Values
    # Reshape z and var into 1D arrays for efficient access using pos_up and pos_down indices defined in the previous step
    # Extract the depth values (z_up and z_down) and variable values (v_up and v_down) at the bracketing levels.
    za=np.reshape(z,N*M*L)
    z_up=za[pos_up]
    z_down=za[pos_down]

    va=np.reshape(var,N*M*L)
    v_up=va[pos_up]
    v_down=va[pos_down]
    
    # Linear Interpolation between two depths
    vnew=mask * ( (v_up-v_down)*depth + v_down*z_up - v_up*z_down ) / (z_up-z_down)

    return vnew

def get_ds(fname,var_str='h'):
    '''
    flexible method to get the xarray dataset for either a
    single or multiple CROCO files 
    '''
    # print('Opening dataset: ' + fname)
    if ('*' in fname) or ('?' in fname) or ('[' in fname):
        # this approach borrowed from OpenDrift's reader_ROMS_native.py
        # our essential vars are the 'var_str' (obviously) plus some other 
        # static vars we need to keep:
        static_vars=['s_rho', 's_w', 'sc_r', 'sc_w', 'Cs_r', 'Cs_w', 
                        'hc', 'angle', 'h', 'f', 'pn', 'pm',
                        'Vtransform','theta_s','theta_b',
                        'lon_rho', 'lat_rho', 'mask_rho',
                        'lon_u', 'lat_u', 'lon_v', 'lat_v',
                        'eta_rho', 'xi_rho', 'eta_v', 'xi_u']
        if var_str in static_vars:
            # No need for open_mfdataset, which can be slower.
            # This is here just in case you want to use get_var() and not
            # have to change fname just to read a static variable
            fname=glob(fname)[0]
            ds = xr.open_dataset(fname, decode_times=False)
        else:
            # let's use open_mfdataset, but drop non-essential vars
            essential_vars=static_vars+[var_str,'time', 'zeta']
            def drop_non_essential_vars_pop(ds):
                dropvars = [v for v in ds.variables if v not in
                            essential_vars]
                ds = ds.drop_vars(dropvars)
                return ds
            ds = xr.open_mfdataset(fname,
                # chunks={'time': 1000}, # limited tests show using chunks can be slower 
                compat='override', 
                decode_times=False,
                preprocess=drop_non_essential_vars_pop,
                data_vars='minimal', 
                coords='minimal', 
                # parallel=True # can actually slow it down in some limited tests!
                )
    else:
        ds = xr.open_dataset(fname, decode_times=False)
    return ds

def get_depths(ds):
    '''
        extract the depth levels (in metres, negative downward) of the sigma levels in a CROCO file(s)
        ds = xarray dataset object read in from CROCO output file(s)
        see get_var()
        the time dimension must be in the ds, even if it is length 1
    '''
    
    ssh=ds.zeta.values
    h = ds.h.values
    
    # get the variables used to calculate the sigma levels
    # CROCO uses these params to determine how to deform the vertical grid
    theta_s = ds.theta_s
    theta_b = ds.theta_b
    hc = ds.hc.values
    N = np.shape(ds.s_rho)[0]
    type_coordinate = "rho"
    vtransform = ds.Vtransform.values
    if not vtransform == 1 and not vtransform == 2:
        raise Exception("Unexpected value for vtransform (" + vtransform + ")")

    depth_rho = z_levels(h, ssh, theta_s, theta_b, hc, N, type_coordinate, vtransform)
    
    # depth_rho = np.squeeze(depth_rho)
    
    # we need to create a new dataarray for the depths of the sigma levels
    depth_da = xr.DataArray(depth_rho, coords={'time': ds['time'].values,
                                         's_rho': ds['s_rho'].values,
                                         'eta_rho': ds['eta_rho'].values, 
                                         'xi_rho': ds['xi_rho'].values
                                         },
                                  dims=['time', 's_rho', 'eta_rho', 'xi_rho'])
    depth_da.attrs['long_name'] = 'Depth of sigma levels of the rho grid (centred in grid cells)'
    depth_da.attrs['units'] = 'meter'
    depth_da.attrs['positive'] = 'up'
    
    return depth_da

def find_nearest_time_indx(dt,dts):
    '''
    dt : array of datetimes
    dts : list of datetimes for which we want to return the nearest indices
    returns corresponding indices

    '''
    
    # dts needs to be list, even if it's a single datetime
    # so the enumerate() loop below will always work
    if isinstance(dts, datetime):
        dts = [dts]
    
    indx_out = np.zeros_like(dts)
    for t, dts_t in enumerate(dts):
        indx_out[t] = np.argmin(np.abs(np.array(dt)-dts_t))

    return indx_out.astype(int)

def get_time(fname,ref_date=None):
    ''' 
        fname = CROCO output file (or file pattern to use when opening with open_mfdataset())
                fname can also be an xarray dataset for enhanced functionality
        ref_date = reference date for the croco run as a datetime object
        
    '''
    if isinstance(fname, xr.Dataset):
        ds = fname.copy()
    else:
        ds = get_ds(fname)

    time_ds = ds.time.values
    
    # convert time from floats to datetimes, if not already converted in the input ds object
    if all(isinstance(item, float) for item in np.atleast_1d(time_ds)):
        
        if ref_date is None:
            print('ref_date is not defined - using default of 2000-01-01')
            ref_date=datetime(2000,1,1)
    
        # convert 'time_ds' (in seconds since ref_date) to a list of datetimes
        time_dt = []
        for t in time_ds:
            date_now = ref_date + timedelta(seconds=np.float64(t))
            time_dt.append(date_now)
    else:
        time_dt = time_ds.astype('datetime64[s]').astype(datetime)
    
    ds.close()
    return time_dt

def get_grd_var(fname,var_str,
                   eta_rho=slice(None),
                   xi_rho=slice(None)):
    '''
    Extract a grid variable on the rho grid
    '''
    
    if isinstance(fname, xr.Dataset) or isinstance(fname, xr.DataArray):
        ds = fname.copy()
    else:
        # for effeciency we shouldn't use open_mfdataset for this function 
        # only use the first file
        if ('*' in fname) or ('?' in fname) or ('[' in fname):
            fname=glob(fname)[0]
        ds = get_ds(fname)
    
    # subset in space
    ds = ds.isel(eta_rho=eta_rho,  xi_rho=xi_rho)
    
    # extract the requested variable - da is a dataarray
    da = ds[var_str]
    
    ds.close()
    
    return da
    
def get_lonlatmask(fname,type='r',
                   eta_rho=slice(None),
                   xi_rho=slice(None)):
    
    lon = get_grd_var(fname,'lon_rho',eta_rho=eta_rho,xi_rho=xi_rho)
    lat = get_grd_var(fname,'lat_rho',eta_rho=eta_rho,xi_rho=xi_rho)
    mask = get_grd_var(fname,'mask_rho',eta_rho=eta_rho,xi_rho=xi_rho)
    
    mask = mask.where(mask != 0, np.nan)
    [Mp,Lp]=mask.shape
    
    # croco output files don't have lon_u,lat_u, mask_u written to them.
    # We could get these from the grid file but I'm rather computing them
    # from lon_rho, lat_rho and mask_rho for the convenience of not having to 
    # specify two input files in functions like get_var()
    
    if type=='u':
        lon=0.5*(lon[:,0:Lp-1]+lon[:,1:Lp])
        lat=0.5*(lat[:,0:Lp-1]+lat[:,1:Lp])
        mask=mask[:,0:Lp-1]*mask[:,1:Lp]
        
    if type=='v':
        lon=0.5*(lon[0:Mp-1,:]+lon[1:Mp,:])
        lat=0.5*(lat[0:Mp-1,:]+lat[1:Mp,:])
        mask=mask[0:Mp-1,:]*mask[1:Mp,:]
        
    return lon,lat,mask

def time_to_slice(time_croco, time):
    '''
    Take the input to get_var, and return a slice object to be used to
    subset the dataset using ds.isel()
    see get_var() for how this is used
    '''
    # check if time input is(are) datetime object(s), 
    # in which case convert it/them into the correct time index/indices
    if isinstance(np.atleast_1d(time)[0],datetime):
        time = find_nearest_time_indx(time_croco,time)
        
    # get the time indices for input to ds.isel()
    if not isinstance(time,slice):
        if isinstance(time,int):
            # make sure time is a slice, even if it's a single integer
            # this is a hack to make sure we keep the time dimension 
            # after the ds.isel() step below, even though it's a single index
            # https://stackoverflow.com/questions/52190344/how-do-i-preserve-dimension-values-in-xarray-when-using-isel
            time = slice(time,time+1) 
        elif len(time)==1:
            # so time is a list with length 1
            time = slice(time[0],time[0]+1) # this will be a slice with a singe number
    
        elif len(time)==2:
            # convert the start and end limits into a slice
            time = slice(time[0],time[1]+1) # +1 to make indices inclusive 
    
    return time

def domain_to_slice(eta_rho,eta_v,xi_rho,xi_u,subdomain,grdname,var_str):
    '''
    Take the domain reltaed input to get_var, and return slice objects to be used to
    subset the dataset using ds.isel()
    '''
    if subdomain is None:
        # as per time, make sure we keep the eta_rho/xi dimensions after the ds.isel() step, even if we specify a single value
        # this greatly simplifies further functions for depth interpolation 
        # as we know the number of dimensions, even if some of them are single length
        # https://stackoverflow.com/questions/52190344/how-do-i-preserve-dimension-values-in-xarray-when-using-isel
        eta_rho = [eta_rho] if not isinstance(eta_rho, slice) else eta_rho
        xi_rho = [xi_rho] if not isinstance(xi_rho, slice) else xi_rho
        eta_v = [eta_v] if not isinstance(eta_v, slice) else eta_v
        xi_u = [xi_u] if not isinstance(xi_u, slice) else xi_u

    else:
        # using subdomain input to do the spatial subset - [lon0,lon1,lat0,lat1]
        # get the indices corresponding to the four corners of the requested subdomain
        j_bl,i_bl = find_nearest_point(grdname, subdomain[0], subdomain[2])
        j_br,i_br = find_nearest_point(grdname, subdomain[1], subdomain[2])
        j_tr,i_tr = find_nearest_point(grdname, subdomain[1], subdomain[3])
        j_tl,i_tl = find_nearest_point(grdname, subdomain[0], subdomain[3])
        # get extreme indices for slicing (using all to be safe in the case of a curvilinear grid)
        j_min = min(j_bl,j_br,j_tr,j_tl)-1
        j_max = max(j_bl,j_br,j_tr,j_tl)+1
        i_min = min(i_bl,i_br,i_tr,i_tl)-1
        i_max = max(i_bl,i_br,i_tr,i_tl)+1
        # create the slice objects
        eta_rho = slice(j_min,j_max+1) # adding one as slice is non-inclusive of the end index
        xi_rho = slice(i_min,i_max+1) # adding one as slice is non-inclusive of the end index
        eta_v = slice(j_min,j_max) # not adding one to keep the correct v-grid indices
        xi_u = slice(i_min,i_max) # not adding one to keep the correct u-grid indices
        
    return eta_rho,eta_v,xi_rho,xi_u # all as slice objects

def level_to_slice(level):
    '''
    Take the input to get_var, and return a slice object to be used to
    subset the dataset using ds.isel()
    see get_var() for how this is used
    '''
    # it gets a bit convoluted for the vertical levels 
    # as we have the option of a constant z level which needs interpolation...
    # so we start by getting a variable 'level_for_isel' which is as it sounds
    if not isinstance(level,slice):
        # so level is a single number or a list of numbers
        level=np.atleast_1d(level).astype('float32') # makes life easier for handling both profiles and time-series if they're both arrays
        if np.mean(level) >= 0: 
            # so we're extracting a single sigma layer or a list of sigma layers
            level_for_isel = slice(int(level[0]), int(level[-1]+1))
            level = level_for_isel
        else:
            # so we're extracting one or more z levels
            # so we'll need to do vertical interpolations later 
            # for this we'll need to initially extract all the sigma levels
            level_for_isel = slice(None)
    else:
        level_for_isel = level # a slice object by definition of the logic
    
    return level_for_isel, level

def get_var(fname,var_str,
            grdname=None,
            time=slice(None),
            level=slice(None),
            eta_rho=slice(None),
            eta_v=slice(None),
            xi_rho=slice(None),
            xi_u=slice(None),
            subdomain=None,
            ref_date=None,
            nc_out=None):
    '''
        extract a variable from a CROCO file
        fname = CROCO output file name (or file pattern to be used with open_mfdataset())
                fname can also be a previously extracted xarray dataset for enhanced functionality
        var_str = variable name (string) in the CROCO output file(s)
        grdname = optional name of your croco grid file (only needed if the grid info is not in fname)
        time = time step indices to extract 
               If slice(None), then all time-steps are extracted
               time can be a single integer (starting at zero) or a datetime object, in which case the nearest time index is extracted
               time can betwo integers/datetimes in a list e.g. [dt1,dt2], in which case the range between the two is extracted
                
        level = vertical level to extract
                if slice(None), then all sigma levels are extracted, and the depths of the levels are provided as an additional variable
                if a positve integer or a list of positive integers, then those sigma levels are extracted (zero denotes the bottom layer, going upward to the surface)
                if a negative number, or a list of negative numers then data are interpolated to those z levels
                (if zero is contained in a list of negative numbers, then it will be treated as the surface)
        eta_rho = index/indices of the eta_rho axis
              If slice(None), then all indices are extracted
        eta_v = index/indices of the eta_v axis
              If slice(None), then all indices are extracted
              this is only needed for extracting a subset of 'v'
        xi_rho = index/indices of the eta_rho axis
              If slice(None), then all indices are extracted
        xi_u = index/indices of the eta_rho axis
              If slice(None), then all indices are extracted
              this is only needed for extracting a subset of 'u'
        subdomain = extents used to do a spatial subset, as a list in format [lon0,lon1,lat0,lat1]
              If None, then no subsetting will get done
        ref_date = reference datetime used in croco runs
        nc_out = option to write a netcdf file from the output dataset
        
        Retruns an xarray dataset object of the requested data
    '''
    
    print('extracting the data from croco file(s) - ' + var_str)

    # Get an xarray dataset of the croco file(s)
    if isinstance(fname, xr.Dataset): # handles the case of using an already extracted dataset as input
        ds = fname.copy()
    else:
        ds = get_ds(fname,var_str)
    
    # Write grid dimensions and variables from grid file, if specfied
    if grdname is not None:
        ds_grd = get_ds(grdname)
        ds = ds.assign_coords(lon_rho = ds_grd['lon_rho'])
        ds = ds.assign_coords(lat_rho = ds_grd['lat_rho'])
        ds = ds.assign_coords(eta_rho = ds_grd['eta_rho'])
        ds = ds.assign_coords(xi_rho = ds_grd['xi_rho'])
        ds['h'] = ds_grd['h']
        ds['mask_rho'] = ds_grd['mask_rho']
        ds_grd.close()
    
    # get the time as a list of datetimes
    time_dt = get_time(ds, ref_date)
    ds = ds.assign_coords(time=time_dt)
    
    # for each of the input dimensions we check the format of the input 
    # and construct the appropriate slice to extract
    time = time_to_slice(time_dt, time)
    eta_rho,eta_v,xi_rho,xi_u = domain_to_slice(eta_rho,eta_v,xi_rho,xi_u,subdomain,ds,var_str)
    level_for_isel,level = level_to_slice(level)
    
    # subset the dataset
    ds = ds.isel(time=time,
                       s_rho=level_for_isel,
                       s_w=level_for_isel,
                       eta_rho=eta_rho,
                       xi_rho=xi_rho,
                       xi_u=xi_u,
                       eta_v=eta_v,
                       missing_dims='ignore' # handle case where input is a previously extracted dataset
                       )
    
    # get dataarrays of the data we want
    da = ds[var_str]
    if len(da.shape)==4:
        var_is_2d=False
    else:
        var_is_2d=True
    h = ds['h']
    zeta = ds['zeta']
    
    # regrid u/v vector components onto the rho grid if needed
    if var_str in ['u','sustr','bustr','ubar'] or var_str in ['v','svstr','bvstr','vbar']:
        print('regridding '+var_str+' onto rho grid')
        if var_str in ['u','sustr','bustr','ubar']:
            data_rho=u2rho(da)   
        if var_str in ['v','svstr','bvstr','vbar']:
            data_rho=v2rho(da) 
    
        # Create a new xarray DataArray with correct dimensions
        # now that u/v data is on the rho grid
        if var_is_2d:
            da_rho = xr.DataArray(data_rho, 
                                  coords={
                                      'time': ds['time'].values,
                                      'eta_rho': ds['eta_rho'].values,
                                      'xi_rho': ds['xi_rho'].values
                                      },
                                  dims=['time', 'eta_rho', 'xi_rho']
                                  )
        else:
            da_rho = xr.DataArray(data_rho, 
                                  coords={
                                      'time': ds['time'].values,
                                      's_rho': ds['s_rho'].values,
                                      'eta_rho': ds['eta_rho'].values,
                                      'xi_rho': ds['xi_rho'].values,
                                      },
                                  dims=['time', 's_rho', 'eta_rho', 'xi_rho']
                                  )

        da = da_rho.copy()
   
    # Do vertical interpolations if needed
    if not var_is_2d and not isinstance(level,slice):
        if np.mean(np.atleast_1d(level)) < 0: # we can't put this in the line above as you can't use '<' on a slice, so at least here we know 'level' is not a slice
            
            print('doing vertical interpolations - ' + var_str)
            # given the above checks in the code, here we should be dealing with a 3D variable 
            # and we want a hz slice at a constant depth level
            z=get_depths(ds) # have to use ds, not da, as we need zeta and h for this
            da_out=hlev_xarray(da, z, level)
            # use the same attributes as the original da
            da_out.attrs = da.attrs
            # update da to be the data for the specified level
            da=da_out.copy()
        
    # Masking
    print('applying the mask - ' + var_str)
    mask = ds.mask_rho
    mask_nan = mask.where(mask != 0, np.nan)
    da = da.squeeze() * mask_nan
    zeta = zeta.squeeze() * mask_nan
    h = h.squeeze() * mask_nan

    # include the depths of the sigma levels in the output
    if 's_rho' in da.coords: # this will include 1 sigma layer - is this an issue?       
        print('computing depths of sigma levels...')
        depths = get_depths(ds).squeeze() * mask_nan
        print('making the output dataset for get_var()...')
        var_data, depth_data, zeta_data, h_data = dask.compute(da, depths, zeta, h)
        ds_out = xr.Dataset({var_str: var_data, 'depth': depth_data, 'zeta': zeta_data, 'h': h_data, 'mask':mask})
        ds_out['s_rho'].attrs.pop('formula_terms', None)
    else:
        print('making the output dataset for get_var()...')
        var_data, zeta_data, h_data = dask.compute(da, zeta, h)
        ds_out = xr.Dataset({var_str: var_data, 'zeta': zeta_data, 'h': h_data, 'mask':mask})

    # remove singleton dimensions
    ds_out = ds_out.squeeze()
    
    # change the attributes to make the dataset cf compliant
    attrs = CROCO_Attrs()
    ds_out['eta_rho'] = change_attrs(attrs,ds_out.eta_rho,'eta_rho')
    ds_out['xi_rho']  = change_attrs(attrs,ds_out.xi_rho,'xi_rho')
    ds_out['lon_rho'] = change_attrs(attrs,ds_out.lon_rho,'lon_rho')
    ds_out['lat_rho']  = change_attrs(attrs,ds_out.lat_rho,'lat_rho')
    ds_out['h'] = change_attrs(attrs,ds_out.h,'h')
    ds_out['mask'] = change_attrs(attrs,ds_out.mask,'mask')
    ds_out['zeta'] = change_attrs(attrs,ds_out.zeta,'zeta')
    ds_out[var_str] = change_attrs(attrs,ds_out[var_str],var_str)
    
    if nc_out is not None:
        print('writing the netcdf file')
        ds_out.to_netcdf(nc_out)
    
    ds.close()
    
    return ds_out

def get_uv(fname,
           grdname=None,
           time=slice(None),
           level=slice(None),
           eta_rho=slice(None),
           eta_v=slice(None),
           xi_rho=slice(None),
           xi_u=slice(None),
           subdomain=None,
           ref_date=None,
           var_u='u', # could also be sustr, bustr, ubar
           var_v='v', # could also be svstr, bvstr, vbar
           nc_out=None):
    '''
    extract u and v components from a CROCO output file(s), regrid onto the 
    rho grid and rotate from grid-aligned to east-north components
    
    see get_var() for a description of the inputs
    
    in addition to the standard get_var() inputs, there are options to sepcify 
    var_u and var_v and strings. Default is to extract the baroclinic velocities
    but you can also extract other variables on the u and v grids
    
    returns xarray dataarrays for both u and v data
    
    '''

    u=get_var(fname,var_u,
              grdname=grdname,
              time=time,
              level=level,
              eta_rho=eta_rho,
              eta_v=eta_v,
              xi_rho=xi_rho,
              xi_u=xi_u,
              subdomain=subdomain,
              ref_date=ref_date)
    v=get_var(fname,var_v,
              grdname=grdname,
              time=time,
              level=level,
              eta_rho=eta_rho,
              eta_v=eta_v,
              xi_rho=xi_rho,
              xi_u=xi_u,
              subdomain=subdomain,
              ref_date=ref_date)
    # get the dataarrays from the datasets
    u_da=u[var_u]
    v_da=v[var_v]
    
    # regridding from the u and v grids to the rho grid is now handled inside 
    # get_var() which allows us to more easily do the vertical interpolation 
    # inside get_var() using the depth levels which are defined on the rho grid
    
    # -------------------
    # Rotate the vectors
    # -------------------
    print('rotating u/v vector components to be east/north components')
    # grid angle
    if grdname is None:
        grdname = fname
    eta_rho,eta_v,xi_rho,xi_u = domain_to_slice(eta_rho,eta_v,xi_rho,xi_u,subdomain,grdname,'angle')
    angle=get_grd_var(grdname, 'angle',
                  eta_rho=eta_rho,
                  xi_rho=xi_rho,
                  )
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    #
    # Refer to https://en.wikipedia.org/wiki/Rotation_matrix
    # although 'angle' is 2D, numpy and xarray are clever enough for this to work even if u_rho and v_rho are 3D or 4D
    u_out = u_da*cos_a - v_da*sin_a
    v_out = v_da*cos_a + u_da*sin_a

    attrs = CROCO_Attrs_RotatedVectors()
    u_out = change_attrs(attrs,u_out,var_u)
    v_out = change_attrs(attrs,v_out,var_v)

    # create a dataset containing both u and v
    ds_out=u # just using u as the basis for the output dataset
    # then overwrite var_u and add var_v
    ds_out = ds_out.assign({var_u: u_out, var_v: v_out})
    
    if nc_out is not None:
        print('writing the netcdf file')
        ds_out.to_netcdf(nc_out)
    
    return ds_out

def get_vort(fname,
             grdname=None,
             time=slice(None),
             level=slice(None),
             ref_date=None):
    '''
    extract the relative vorticity from a CROCO output file:
    dv/dx - du/dy
    
    see get_var() for a description of the inputs   
    
    subsetting in space not perimitted for this. get_var does in fact allow for you to compute
    vorticity on a subset easily... just need to implement here
    
    get_var() functionality has been massively updated since this function was last used
    so it is guarenteed to need some edits to get it to work again
    
    '''
    
    # start by getting u and v
    # and we'll leave them on their native grids for this calc
    # (i.e. intentionally not regridding to the rho grid)
    u=get_var(fname,'u',grdname=grdname,time=time,level=level,ref_date=ref_date)
    v=get_var(fname,'v',grdname=grdname,time=time,level=level,ref_date=ref_date)
    if grdname is None:
        grdname = fname
    pm=get_grd_var(grdname, 'pm') # 1/dx on the rho grid
    pn=get_grd_var(grdname, 'pn') # 1/dy on the rho grid
    
    # this code was taken largely from croco_tools-v1.1/croco_pyvisu/derived_variables.py
    #
    # interpolate pm from rho grid onto psi grid
    dxm1 = 0.25 * (pm[:-1, 1:] + pm[1:, 1:] + pm[:-1, :-1] + pm[1:, :-1]) 
    # Compute d(v)/d(xi) on the psi grid
    # (this should work regardless of whether v is 4D, 3D or 2D)
    dvdxi = np.diff(v, n=1, axis=-1) * dxm1 # axis = -1 will always be the xi dimension
    
    # interpolate pn from rho grid onto psi grid
    dym1 = 0.25 * (pn[:-1, 1:] + pn[1:, 1:] + pn[:-1, :-1] + pn[1:, :-1]) 
    # Compute d(u)/d(eta) on the psi grid
    # (this should work regardless of whether v is 4D, 3D or 2D)
    dudeta = np.diff(u, n=1, axis=-2) * dym1 # axis = -2 will always be the eta dimension
    
    # Vortivity on the psi grid
    vort=dvdxi-dudeta
    
    # regrid to be on the rho grid
    vort=psi2rho(vort)
    
    # TODO: rather return an xarray dataarray 
    
    return vort

def get_boundary(grdname):
    '''
    Return lon,lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    lon_rho,lat_rho,_=get_lonlatmask(grdname,type='r')
    lon = np.hstack((lon_rho[0:, 0], lon_rho[-1, 1:-1],
                     lon_rho[-1::-1, -1], lon_rho[0, -2::-1]))
    lat = np.hstack((lat_rho[0:, 0], lat_rho[-1, 1:-1],
                     lat_rho[-1::-1, -1], lat_rho[0, -2::-1]))
    return lon, lat

def find_nearest_point(grdname, Longi, Latit, Bottom=None):
    """
    Find the nearest indices of the model rho grid to a specified lon, lat coordinate:
            
    Parameters:
    - fname :filename of the model, or the grid file
    - Longi :longitude
    - Latit :latitude
    - Bottom (positive value): if the model bathy is slightly different. This Option to find nearest
      lat and lon in water that is as deep as reference. If == None then
      this looks for only the closest horizontal point.

    Returns:
    - j :the nearest eta index
    - i :the nearest xi index
    
    j,i can be used in xarrays built-in xr.isel() function to extract data at this grid point
    
    (Note j,i aren't the eta_rho,xi_rho values! 
     The j,i indices are one less than the eta_rho,xi_rho values
     because eta_rho,xi_rho are 1 based)

    
    """
    
    lon_rho = get_grd_var(grdname, 'lon_rho').values
    lat_rho = get_grd_var(grdname, 'lat_rho').values
    h = get_grd_var(grdname, 'h').values

    # Calculate the distance between (Longi, Latit) and all grid points
    distance = ((lon_rho - Longi) ** 2 +
                (lat_rho - Latit) ** 2) ** 0.5
    
    if Bottom is None:
        mask=h/h

    else:
        mask=h
        mask[mask<Bottom]=10000
        mask[mask<10000]=1
    
    distance_mask=distance*mask

    # Find the indices of the minimum distance
    # unravel_index method Converts a flat index or array of flat indices into a tuple of coordinate 
    # arrays: https://numpy.org/doc/stable/reference/generated/numpy.unravel_index.html
    min_index = np.unravel_index(distance_mask.argmin(), distance_mask.shape)

    j, i = min_index

    return j, i

def dist_spheric(lat1,lon1,lat2,lon2,R=6367442.76):
    l=np.abs(lon2-lon1)
    if ((np.size(lat1)>1) | (np.size(lat2)>1)):
        l[l>=180]=360-l[l>=180]
    else:
        if l>180:
            l=360-l
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    l=np.radians(l)
    dist=R*np.arctan2(np.sqrt((np.sin(l)*np.cos(lat2))**2 +(np.sin(lat2)*np.cos(lat1)
                     -np.cos(lat2)*np.sin(lat1)*np.cos(l))**2 ),
                          np.sin(lat2)*np.sin(lat1)+ np.cos(lat2)*np.cos(lat1)*np.cos(l)
                          )
    return dist

def find_fractional_eta_xi(grdname,lons,lats):
    '''
    extract the fractional eta_rho,xi_rho indices corresponding to 
    the provided lon,lat coordinates
    
    inputs:
        grdname - croco grid file
        lons  - input longitudes
        lats  - input latitudes
    
    Unlike find_nearest_point(), here we return the eta_rho, xi_rho values, not the indices
    This is useful as we can use these as input to xarrays built-in xr.interp() function
    
    '''
    
    lons, lats = np.atleast_1d(lons), np.atleast_1d(lats)  # Ensure inputs are arrays
    
    lon_rho = get_grd_var(grdname, 'lon_rho').values
    lat_rho = get_grd_var(grdname, 'lat_rho').values
    eta_rho = get_grd_var(grdname, 'eta_rho').values
    xi_rho = get_grd_var(grdname, 'xi_rho').values
    angle = get_grd_var(grdname, 'angle').values
    dx = 1/get_grd_var(grdname, 'pm').values
    dy = 1/get_grd_var(grdname, 'pn').values
    eta_fracs = []
    xi_fracs = []
    for lon, lat in zip(lons, lats):
        # Compute distance to find the nearest grid cell
        lonlat_diff = np.sqrt((lon_rho - lon) ** 2 + (lat_rho - lat) ** 2)
        j, i = np.unravel_index(np.argmin(lonlat_diff), lonlat_diff.shape)  # Closest grid point
        
        eta_nearest=eta_rho[j]
        xi_nearest=xi_rho[i]
        
        # find the residual lon,lat around this point
        dlon = lon - lon_rho[j,i]
        dlat = lat - lat_rho[j,i]
        
        # convert dlon,dlat to a distance in m and an angle relative to TN
        distance = dist_spheric(lat_rho[j,i],lon_rho[j,i],lat,lon)
        angle_TN = (270 - np.rad2deg(np.arctan2(-dlat, -dlon))) % 360  # Angle in radians measured clockwise from true north, toward which the vector is pointing
            
        # change the angle so it is now relative to the rotated model axis (the input angle needs to be added to the TN angle)
        angle_grid = np.rad2deg(angle[j,i]) + angle_TN # this is now the angle clockwise from the eta_rho axis
        
        # use the distance and the angle relative to the grid to compute the grid aligned distance vector components
        distance_eta = distance * np.cos(np.deg2rad(angle_grid))
        distance_xi = distance * np.sin(np.deg2rad(angle_grid))
        
        # compute the fractional eta, xi offsets, based on the grid size dx, dy
        deta = distance_eta / dy[j,i]
        dxi = distance_xi / dx[j,i]
        
        # compute the fractional eta and xi to return
        eta_fracs.append(eta_nearest + deta)
        xi_fracs.append(xi_nearest + dxi)

    return eta_fracs, xi_fracs


def get_ts_multivar(fname, lon, lat, ref_date, 
                grdname=None,
                vars = ['temp','salt'],
                i_shift=0, j_shift=0, 
                time=slice(None),
                level=slice(None),
                nc_out=None
                ):
    """
    Convenience function to get multiple variables of interest into a single xarray dataset/ nc file.
    By default zeta, temp, salt, u and v are extracted, but you can add 'rho' variables to the 'vars' input if you like
                   
    Parameters:
    - see get_ts() for a description of the inputs - they're the same
            
    Returns:
    - ds, an xarray dataset containing the time-series or profile data
    """
    
    # Initialize an empty list to store datasets
    all_datasets = []
    for var in vars:
        ds_var = get_ts(fname, var, lon, lat, ref_date, 
                        grdname=grdname,
                        i_shift=i_shift, j_shift=j_shift, 
                        time=time,
                        level=level)
        all_datasets.append(ds_var)
    # add u,v
    ds_uv = get_ts_uv(fname, lon, lat, ref_date, 
                    grdname=grdname,
                    i_shift=i_shift, j_shift=j_shift, 
                    time=time,
                    level=level)
    all_datasets.append(ds_uv)
    
    # merge into a single dataset
    ds_all = xr.merge(all_datasets)
    
    if nc_out is not None:
        print('writing the netcdf file')
        ds_all.to_netcdf(nc_out)
    
    return ds_all

def get_ts(fname, var_str, lon, lat, ref_date,
                grdname=None,
                i_shift=0, j_shift=0,
                time=slice(None),
                level=slice(None),
                nc_out=None,
                Bottom=None
                ):
    """
           Extract a ts from the model:
                   
            Parameters:
            see get_var() for a description of the common inputs. This function has a few additional ones:
            - lat               :latitude of time-series
            - lon               :longitude of time-series
            - i_shift           :number of grid cells to shift along the xi axis, useful if input lon,lat is on land mask or if input depth is deeper than model depth 
            - j_shift           :number of grid cells to shift along the eta axis, (similar utility to i_shift)
            - Bottom (positive value): if the model bathy is slightly different. This Option to find nearest
              lat and lon in water that is as deep as reference. If == None then
              this looks for only the closest horizontal point.
            
            Returns:
            - ds, an xarray dataset containing the ts data
    """
    if var_str in ['u','sustr','bustr','ubar'] or var_str in ['v','svstr','bvstr','vbar']:
        print('WARNING: '+var_str+' will be the grid aligned vector component')
        print('rather use get_ts_uv() for extracting a time-series of u/v data which represents east/north components')
    
    #find_nearest_point finds the nearest point in the model to the model grid lon, lat extracted from the model grid input.
    if grdname is None:
        grdname = fname
    j, i = find_nearest_point(grdname, lon, lat, Bottom) 
    
    # apply the shifts along the xi and eta axis
    i = i+i_shift
    j = j+j_shift
    
    ds = get_var(fname, var_str,
                          grdname=grdname,
                          time=time,
                          level=level,
                          eta_rho=j,
                          xi_rho=i,
                          ref_date=ref_date,
                          nc_out=nc_out)
        
    return ds

def get_ts_uv(fname, lon, lat, ref_date, 
                grdname=None,
                i_shift=0, j_shift=0, 
                time=slice(None),
                level=slice(None),
                default_to_bottom=False,
                var_u='u', # could also be sustr, bustr, ubar
                var_v='v', # could also be svstr, bvstr, vbar
                nc_out=None,
                Bottom=None
                ):
    """
           Extract a time-series or profile of u,v from the model
           u,v are rotated from grid-aligned to east-north components
           Data are extracted from the nearest 'rho' grid point and u,v data 
           either side of this rho grid point are interpolated to the rho grid point
                              
            Parameters:
            - see get_ts() for description of inputs
            - the only difference is 'var' is not defined here - u and v are extracted by definition
                    
            Returns:
            - ds, an xarray dataset containing the time-series or profile data
              
    """
    
    # finds the rho grid indices nearest to the input lon, lat
    if grdname is None:
        grdname = fname
    j, i = find_nearest_point(grdname, lon, lat,Bottom) 
    
    # apply the shifts along the xi and eta axis
    i = i+i_shift
    j = j+j_shift
    
    # j,i are the eta,xi indices of the rho grid for our time-series extraction
    # extracting u and v at this location is kind of complicated by the u, and v grids
    # to get around this, let's consider a 3x3 block of rho grid points 
    # with j,i in the centre. Then let's define the indices we'll need to extract
    # for the rho grid, u grid and v grid within this 3x3 block of rho grid cells    
    i_rho=slice(i-1,i+2) # 3 indices for the xi_rho axis, with i in the middle
    j_rho=slice(j-1,j+2) # 3 indices for the eta_rho axis, with j in the middle
    i_u=slice(i-1,i+1) # 2 indices for the xi_u axis, either side of i
    j_v=slice(j-1,j+1) # 2 incidces for the eta_v axis, either side of j 
            
    ds = get_uv(fname,
                grdname=grdname,
                time=time,
                level=level,
                eta_rho=j_rho,
                xi_rho=i_rho,
                eta_v=j_v,
                xi_u=i_u,
                var_u=var_u,
                var_v=var_v,
                ref_date=ref_date
                )
    
    # pull out the middle data point from our 3x3 block of rho grid points
    # this is by definition the grid cell we are interested in
    ds = ds.isel(eta_rho=1,xi_rho=1).squeeze()
    
    if nc_out is not None:
        print('writing the netcdf file')
        ds.to_netcdf(nc_out)
    
    return ds

def get_section_coords(lon0, lat0, lon1, lat1, dgc, R=6367442.76):
    """Computes points along a great-circle path with approximately equal distances.

    Args:
        lon0, lat0: Start point in degrees.
        lon1, lat1: End point in degrees.
        dgc: Distance between points in meters.
        R: Earth radius in meters.

    Returns:
        lons, lats, distances: Arrays of longitudes, latitudes, and cumulative distances.
    """
    # Convert inputs to radians
    lat0, lon0, lat1, lon1 = map(np.radians, [lat0, lon0, lat1, lon1])

    # Compute total great-circle distance using the Haversine formula
    delta_lat = lat1 - lat0
    delta_lon = lon1 - lon0
    a = np.sin(delta_lat / 2)**2 + np.cos(lat0) * np.cos(lat1) * np.sin(delta_lon / 2)**2
    sigma_total = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))  # Central angle
    dist_total = R * sigma_total  # Convert to meters

    # Determine the number of segments
    npsec = int(dist_total // dgc) + 1  # Number of points

    # Compute the initial bearing
    X = np.cos(lat1) * np.sin(delta_lon)
    Y = np.cos(lat0) * np.sin(lat1) - np.sin(lat0) * np.cos(lat1) * np.cos(delta_lon)
    initial_bearing = np.arctan2(X, Y)

    # Generate points along the great-circle path
    lons, lats, distances = [lon0], [lat0], [0]
    
    for i in range(1, npsec):
        sigma = i * dgc / R  # Angular distance along the sphere

        lat_new = np.arcsin(np.sin(lat0) * np.cos(sigma) +
                            np.cos(lat0) * np.sin(sigma) * np.cos(initial_bearing))
        lon_new = lon0 + np.arctan2(np.sin(initial_bearing) * np.sin(sigma) * np.cos(lat0),
                                    np.cos(sigma) - np.sin(lat0) * np.sin(lat_new))

        lats.append(lat_new)
        lons.append(lon_new)
        distances.append(i * dgc)  # Keep cumulative distance

    # Convert back to degrees
    lons, lats = np.degrees(lons), np.degrees(lats)

    return np.array(lons), np.array(lats), np.array(distances)

def get_section(fname,
                var_str,
                section_start,
                section_end,
                grdname=None,
                time=slice(None),
                level=slice(None),
                ref_date=None,
                res=None,
                nc_out=None,
                ):
    """
    Extract a vertical section from a CROCO output file(s) 
    The transect can be in any direction. Multiple files can be loaded in.

    Inputs:
    see get_var() for a description of some of the inputs. In addition to the get_var inputs there is:
    section_start = Start point of transect (list; eg. section_start = [lon0, lat0])
    section_end = End points of transect (list; eg. section_end = [lon1, lat1])
    res = Horizontal resolution of the section in meters (eg. res = 300). Default is None in which case 
         it takes the smallest grid size as the resolution. 
    ...
           
    
    Returns:
    - ds, an xarray dataset containing the section data for the variable
    """
    
    # if the gridname is not provided, we use the fname for the gridname 
    if grdname is None:
        grdname = fname
    
    # Define the start and end points of the transect
    lon0,lat0 = section_start[0],section_start[1]
    lon1,lat1 = section_end[0],section_end[1]
    
    # extract the data for the defined subdomain which covers the section extents
    subdomain = [lon0,lon1,lat0,lat1]
    ds = get_var(fname,var_str,
                 grdname=grdname,
                 time=time,
                 level=level,
                 subdomain=subdomain,
                 ref_date=ref_date)
    
    # if the grid resolution is not provided, we use the smallest grid size as the resolution. 
    if res is None:
        dy_min=np.min(1/get_grd_var(grdname, 'pn').values[:])
        dx_min=np.min(1/get_grd_var(grdname, 'pm').values[:])
        res = min(dy_min,dx_min)
    
    # Make the transect on which the interpolation will take place    
    section_lons,section_lats,section_dist = get_section_coords(lon0,lat0,lon1,lat1,res)
    
    # Compute fractional eta_rho/xi_rho indices for all lon/lat pairs in the section
    print('mapping section lon,lat pairs to fractional eta_rho,xi_rho indices...')
    eta_fracs, xi_fracs = find_fractional_eta_xi(grdname,section_lons, section_lats)
    
    # Interpolate the variable along the line
    print('interpolating along the section...')
    ds = ds.interp(eta_rho=("points", eta_fracs), xi_rho=("points", xi_fracs))
    # I'm aware that there is a slight mismatch between the ds.lon_rho, ds.lat_rho and
    # section_lons, section_lats, while theoretically they should be identical
    # this is due to how the fractional eta, xi are interpolated in find_fractional_eta_xi
    # (I tried many different approaches to minimise the error... it is a bit of a head scratcher and maybe could be improved?)
    # The error is however much less than the model grid size, so I am not too bothered by this
    
    # add the section distance to the ds
    ds["distance"] = xr.DataArray(
        section_dist, dims=("points",),
        coords={"points": ds.coords["points"]},
        name="distance_along_section",
        attrs={
            "units": "meter",
            "standard_name": "distance",
        })
    
    if nc_out is not None:
        print('')
        print(f'File saved to: {nc_out}')
        print('')
        ds.to_netcdf(nc_out)
    
    return ds


def compute_mhw(fname_clim, fname_in, fname_out,
                ref_date="2000-01-01",
                use_constant_clim=False):
    """
    Detect marine heatwaves (MHW) using the 90th percentile threshold from monthly climatology.
    Outputs:
    - temp_mhw_mask (bool)
    - temp_anomaly (CROCO-style)
    - temp (SST) original
    
    Parameters
    ----------
    fname_clim : str
        Path to the monthly climatology NetCDF file (12 time steps).
    fname_in : str
        Path to high-frequency CROCO model output NetCDF file (e.g., hourly).
    fname_out : str
        Output file path for the MHW binary mask and anomaly.
    ref_date : str
        Reference date used to convert HF time units.
    use_constant_clim : bool
        If True, uses constant climatology value interpolated to midpoint of HF time.
        If False, interpolates across all HF time steps.
    """

    start_time = time.time()
    print("🔹 Loading climatology and high-frequency files...")
    ds_clim = xr.open_dataset(fname_clim, decode_times=False)["temp"]
    ds_hf = xr.open_dataset(fname_in, decode_times=False)
    ref_hf = np.datetime64(ref_date)

    HF_t = ds_hf.time.values.astype("float64")
    hf_dates = ref_hf + HF_t.astype("timedelta64[s]")
    ds_hf["time"] = hf_dates

    if len(ds_clim.time) != 12:
        raise ValueError("Climatology must have exactly 12 time steps (monthly)")

    print("🔹 Extending climatology for time interpolation...")
    ds_clim = ds_clim.transpose("time", ...)
    start = ds_clim.isel(time=0)
    end = ds_clim.isel(time=-1)
    ds_clim_ext = xr.concat([end, ds_clim, start], dim="time")

    hf_year = pd.to_datetime(hf_dates[0]).year
    clim_time = pd.date_range(start=f"{hf_year - 1}-12-15", periods=14, freq="MS")
    clim_seconds = ((clim_time - ref_hf) / np.timedelta64(1, "s")).astype("float64")
    ds_clim_ext = ds_clim_ext.assign_coords(time=clim_seconds)

    print("🔹 Computing 90th percentile threshold from climatology...")
    clim_p90 = ds_clim_ext.quantile(0.9, dim="time")

    print("🔹 Interpolating threshold to HF time axis...")
    if use_constant_clim:
        middle_time = ds_hf.time.isel(time=int(len(ds_hf.time) / 2))
        clim_interp = clim_p90.expand_dims(time=[middle_time])
        clim_interp = clim_interp.broadcast_like(ds_hf["temp"])
    else:
        clim_interp = clim_p90.expand_dims(time=ds_hf.time)

    print("🔹 Detecting marine heatwaves...")
    mhw_mask = ds_hf["temp"] > clim_interp
    mhw_mask.name = "temp_mhw_mask"
    mhw_mask.attrs.update({
        "long_name": "Marine Heatwave Binary Mask",
        "standard_name": "marine_heatwave_event",
        "units": "1 (true/false)"
    })

    print("🔹 Computing SST anomaly...")
    sst_anom = ds_hf["temp"] - clim_interp
    sst_anom.name = "temp_anom"

    # Apply CROCO attributes to anomaly
    print("🔹 Adding CROCO-style attributes to anomaly...")
    croco_attrs = CROCO_Attrs()
    sst_anom = change_attrs(croco_attrs, sst_anom, "temp_anom")

    print("🔹 Preparing output dataset...")
    ds_out = xr.Dataset(coords=ds_hf.coords)
    ds_out["temp"] = ds_hf["temp"]  # Save original SST
    # ds_out["temp"] = ds_hf["temp"].isel(s_rho=-1)  # assuming s_rho is the vertical dim
    # ds_out["temp"] = ds_out["temp"].expand_dims(dim={"s_rho": [ds_hf.s_rho.values[-1]]})  # add singleton s_rho dimension
    ds_out["temp_anom"] = sst_anom
    ds_out["temp_mhw_mask"] = mhw_mask

    # Add vertical grid and metadata
    add_vars = ['theta_s', 'theta_b', 'hc', 'Vtransform', 'h', 'zeta', 'mask_rho']
    for var in add_vars:
        if var in ds_hf:
            ds_out[var] = ds_hf[var]
        elif var in ds_hf.attrs:
            ds_out.attrs[var] = ds_hf.attrs[var]

    print("💾 Writing output file...")
    encoding = {var: {"dtype": "float32"} for var in ds_out.data_vars}
    ds_out.to_netcdf(fname_out, encoding=encoding, mode="w")

    print(f"✅ Done! Output saved to: {fname_out}")
    print(f"⏱️ Total time elapsed: {time.time() - start_time:.2f} seconds")







def compute_anomaly(fname_clim, fname_in, fname_out,
                    ref_date="2000-01-01",
                    varlist=["temp", "u", "v", "salt", "zeta"],
                    use_constant_clim=False):
    """
    Compute anomalies by subtracting monthly climatology from high-frequency CROCO output.

    Parameters:
    -----------
    fname_clim : str
        Path to the NetCDF file containing 12 monthly climatology time steps.
    fname_in : str
        Path to the NetCDF file containing high-frequency model output (e.g., hourly).
    fname_out : str
        Path to output NetCDF file containing the anomalies
    ref_date : str, optional
        Reference date (e.g., "2000-01-01") used to define the high-frequency time axis.
    varlist : list of str, optional
        List of variables for which anomalies should be computed.
     use_constant_clim : bool, optional
        If True, use a constant climatology (interpolated to midpoint of HF time axis insted of to the entire HF time axis).
        This is useful when the high frequency file is very large, but a single climatology value is sufficient
        An alternative approach to this problem is to chunk the data. This works to speed up the interpolation step
        but leads to memory issue when trying to write the output (at least I couldn't solve them without bypassing the problem)

    Returns:
    --------
    None. Saves output NetCDF with anomaly variables added.
    """

    start_time = time.time()
    print("Loading climatology and high-frequency files...")
    ds_clim = xr.open_dataset(fname_clim, decode_times=False)
    ds_clim = ds_clim[varlist] # subset to only the variables we need
    ds_hf = xr.open_dataset(fname_in, decode_times=False)
    ref_hf = np.datetime64(ref_date)

    # Ensure climatology file has 12 time steps
    if len(ds_clim.time) != 12:
        raise ValueError("ERROR: Provided climatology file must have exactly 12 monthly time steps.")

    print("Extending climatology time axis...")
    start = ds_clim.isel(time=0)
    end = ds_clim.isel(time=-1)
    ds_clim = xr.concat([end, ds_clim, start], dim="time")
    ds_clim = ds_clim.transpose('time', ...) # ensure time is the first dimension
    
    print("Set climatology time axis to align with high frequency file...")
    HF_t = ds_hf.time.values.astype("float64")
    hf_dates = ref_hf + HF_t.astype("timedelta64[s]")
    hf_year = pd.to_datetime(hf_dates[0]).year
    clim_time = pd.date_range(start=f'{hf_year-1}-12-15', periods=14, freq='MS')
    clim_seconds = ((clim_time - ref_hf) / np.timedelta64(1, "s")).to_numpy(dtype=np.float64) # seconds since initialisation i.e. as per the high frequency file
    ds_clim = ds_clim.assign_coords(time=clim_seconds)
    if use_constant_clim:
        print("(Using constant climatology based on midpoint of HF time axis)")
        middle_time = ds_hf.time.isel(time=int(len(ds_hf.time) / 2))
        ds_clim = ds_clim.interp(time=middle_time, method="linear")
    else:
        ds_clim = ds_clim.interp(time=ds_hf.time, method="linear")
    
    print("Computing anomalies...")
    ds_anom = xr.Dataset(coords=ds_hf.coords)
    croco_attrs = CROCO_Attrs()
    for var in varlist:
        print(f"{var}")
        anom = ds_hf[var] - ds_clim[var]
        anom_name = f"{var}_anom"
        anom = change_attrs(croco_attrs, anom, anom_name)
        ds_anom[anom_name] = anom
    ds_hf.close()
    ds_clim.close()
    
    # add variables from ds_hf related to the vertical grid
    # (this is needed if you want to plot or extract data at specific vertical levels later)
    add_vars = ['theta_s','theta_b','hc','Vtransform','h','zeta','mask_rho']
    for add_var in add_vars:
        if add_var in ds_hf:
            ds_anom[add_var] = ds_hf[add_var]
        elif add_var in ds_hf.attrs: # handle 'theta_s' and 'theta_b' which are stored as global attributes not variables
            ds_anom.attrs[add_var] = ds_hf.attrs[add_var]
    
    print("Writing output file...")   
    encoding = {
        var: {"dtype": "float32"}
        for var in ds_anom.data_vars
    }
    ds_anom.to_netcdf(fname_out, encoding=encoding, mode='w')

    end_time = time.time()
    print(f"Total time elapsed: {end_time - start_time:.2f} seconds")