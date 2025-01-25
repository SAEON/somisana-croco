import numpy as np
from datetime import timedelta
import xarray as xr
import dask
from datetime import timedelta, datetime
from glob import glob
import sys
from scipy.interpolate import RegularGridInterpolator

def u2rho(u):
    """
    regrid the croco u-velocity from it's native u grid to the rho grid
    u can be 2D, 3D or 4D numpy array or an xarray dataarray
    returns a numpy array of the data on the rho grid
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
    returns a numpy array of the data on the rho grid
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
    this provides a 3D grid of the depths of the sigma levels
    h = 2D bathymetry of your grid
    zeta = zeta at particular timestep that you are interested in
    theta_s = surface stretching parameter
    theta_b = bottom stretching parameter
    hc = critical depth
    N = number of sigma levels
    type = 'w' or 'rho'
    vtransform = 1 (OLD) or 2 (NEW)

    this is adapted (by J.Veitch - Feb 2022) from zlevs.m in roms_tools (by P. Penven)
    """

    [M, L] = np.shape(h)

    sc_r = np.zeros((N, 1))
    Cs_r = np.zeros((N, 1))
    sc_w = np.zeros((N + 1, 1))
    Cs_w = np.zeros((N + 1, 1))

    if vtransform == 2:
        ds = 1 / N

        if type == "w":
            sc_r[0, 0] = -1.0
            sc_w[N, 0] = 0
            Cs_w[0, 0] = -1.0
            Cs_w[N, 0] = 0

            sc_w[1:-1, 0] = ds * (np.arange(1, N, 1) - N)

            Cs_w = csf(sc_w, theta_s, theta_b)
            N = N + 1
        else:
            sc = ds * (np.arange(1, N + 1, 1) - N - 0.5)
            Cs_r = csf(sc, theta_s, theta_b)
            sc_r = sc

    else:
        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5 * theta_s)

        if type == "w":
            sc = (np.arange(0, N + 1, 1) - N) / N
            N = N + 1
        else:
            sc = (np.arange(1, N + 1, 1) - N - 0.5) / N

        Cs = (1.0 - theta_b) * cff1 * np.sinh(theta_s * sc) + theta_b * (
            cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5
        )

    h[h == 0] = 1e-2
    Dcrit = 0.01
    zeta[zeta < (Dcrit - h)] = Dcrit - h[zeta < (Dcrit - h)]
    hinv = 1 / h

    z = np.zeros((N, M, L))

    if vtransform == 2:
        if type == "w":
            cff1 = Cs_w
            cff2 = sc_w + 1
            sc = sc_w
        else:
            cff1 = Cs_r
            cff2 = sc_r + 1
            sc = sc_r
        h2 = h + hc
        cff = hc * sc
        h2inv = 1 / h2

        for k in np.arange(N, dtype=int):
            z0 = cff[k] + cff1[k] * h
            z[k, :, :] = z0 * h / (h2) + zeta * (1.0 + z0 * h2inv)
    else:
        cff1 = Cs
        cff2 = sc + 1
        cff = hc * (sc - Cs)
        cff2 = sc + 1
        for k in np.arange(N, dtype=int):
            z0 = cff[k] + cff1[k] * h
            z[k, :, :] = z0 + zeta * (1.0 + z0 * hinv)

    return z

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
    # Convert depth to xarray DataArray if it's a scalar or a list
    if np.isscalar(depth):
        depth = xr.DataArray([depth], dims="depth")
    else:
        depth = xr.DataArray(depth, dims="depth")
        
    # Add attributes to the depth coordinate
    depth.attrs["long_name"] = "water depth from free surface"
    depth.attrs["units"] = "meters"
    depth.attrs["positive"] = "up"
    depth.attrs["bottom"] = "-99999 denotes the bottom layer of the model"
    
    # Ensure z and var have the same coordinate details
    # (we need to ensure this so that xarray can map the two dataarrays properly)
    # TODO: this should in theory already be handled in get_depths - the z
    # input to this function should already have identical coordinates to var
    z = z.assign_coords(time=var.coords["time"],
                        eta_rho=var.coords["eta_rho"],
                        xi_rho=var.coords["xi_rho"])

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

def get_ds(fname,var_str=''):
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

    T,M,L = np.shape(ssh)
    depth_rho = np.zeros((T,N,M,L))
    for x in np.arange(T):
        depth_rho[x, :, :, :] = z_levels(
            h, ssh[x, :, :], theta_s, theta_b, hc, N, type_coordinate, vtransform
        )
    
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

def get_time(fname,ref_date=None,time_lims=slice(None)):
    ''' 
        fname = CROCO output file (or file pattern to use when opening with open_mfdataset())
                fname can also be an xarray dataset for enhanced functionality
        ref_date = reference date for the croco run as a datetime object
        time_lims = optional list of two datetimes i.e. [dt1,dt2], which define the range of times to extract
                    If slice(None), then all time-steps are extracted
    '''
    if isinstance(fname, xr.Dataset) or isinstance(fname, xr.DataArray):
        ds = fname.copy()
    else:
        ds = get_ds(fname)

    time = ds.time.values
    
    # convert time from floats to datetimes, if not already converted in the input ds object
    if all(isinstance(item, float) for item in np.atleast_1d(time)):
        
        if ref_date is None:
            print('ref_date is not defined - using default of 2000-01-01')
            ref_date=datetime(2000,1,1)
    
        # convert 'time' (in seconds since ref_date) to a list of datetimes
        time_dt = []
        for t in time:
            date_now = ref_date + timedelta(seconds=np.float64(t))
            time_dt.append(date_now)
    else:
        time_dt = time.astype('datetime64[s]').astype(datetime)
    
    # subset based in time_lims input
    if not isinstance(time_lims,slice):
        time_lims = tstep_to_slice(fname, time_lims, ref_date)
        
    time_dt = time_dt[time_lims]
    
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
    
    lon = get_grd_var(fname,'lon_rho',eta_rho=eta_rho,xi_rho=xi_rho).values
    lat = get_grd_var(fname,'lat_rho',eta_rho=eta_rho,xi_rho=xi_rho).values
    mask = get_grd_var(fname,'mask_rho',eta_rho=eta_rho,xi_rho=xi_rho).values
    
    mask[np.where(mask == 0)] = np.nan
    [Mp,Lp]=mask.shape;
    
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

def tstep_to_slice(fname, tstep, ref_date):
    '''
    Take the input to get_var, and return a slice ofbject to be used to
    subset the dataset using ds.isel()
    see get_var() for how this is used
    '''
    # check if tstep input is instance of datetime, 
    # in which case convert it/them into the correct time index/indices
    if isinstance(np.atleast_1d(tstep)[0],datetime):
        if ref_date is None:
            print('ref_date is not defined - using default of 2000-01-01')
            ref_date=datetime(2000,1,1)
        time_croco = get_time(fname,ref_date) # get_time actually calls tstep_to_slice (CIRCULAR!), but only inside an if statement which won't be entered with this input. MESSY. Should do better
        tstep = find_nearest_time_indx(time_croco,tstep)
        
    # get the time indices for input to ds.isel()
    if not isinstance(tstep,slice):
        if isinstance(tstep,int):
            # make sure tstep is a slice, even if it's a single integer
            # this is a hack to make sure we keep the time dimension 
            # after the ds.isel() step below, even though it's a single index
            # https://stackoverflow.com/questions/52190344/how-do-i-preserve-dimension-values-in-xarray-when-using-isel
            tstep = slice(tstep,tstep+1) 
        elif len(tstep)==1:
            # so tstep is a list with length 1
            tstep = slice(tstep[0],tstep[0]+1) # this will be a slice with a singe number
    
        elif len(tstep)==2:
            # convert the start and end limits into a slice
            tstep = slice(tstep[0],tstep[1]+1) # +1 to make indices inclusive 
    
    return tstep

def domain_to_slice(eta_rho,eta_v,xi_rho,xi_u,subdomain,grdname,var_str):
    '''
    Take the domain reltaed input to get_var, and return slice objects to be used to
    subset the dataset using ds.isel()
    '''
    if subdomain is None:
        eta_rho = eta_or_xi_to_slice(eta_rho,var_str)
        xi_rho = eta_or_xi_to_slice(xi_rho,var_str)
        eta_v = eta_or_xi_to_slice(eta_v,var_str)
        xi_u = eta_or_xi_to_slice(xi_u,var_str)
    else:
        # using subdomain input to do the spatial subset - [lon0,lon1,lat0,lat1]
        # get the indices corresponding to the four corners of the requested subdomain
        j_bl,i_bl = find_nearest_point(grdname, subdomain[0], subdomain[2])
        j_br,i_br = find_nearest_point(grdname, subdomain[1], subdomain[2])
        j_tr,i_tr = find_nearest_point(grdname, subdomain[1], subdomain[3])
        j_tl,i_tl = find_nearest_point(grdname, subdomain[0], subdomain[3])
        # get extreme indices for slicing (using all to be safe in the case of a curvilinear grid)
        j_min = min(j_bl,j_br,j_tr,j_tl)
        j_max = max(j_bl,j_br,j_tr,j_tl)
        i_min = min(i_bl,i_br,i_tr,i_tl)
        i_max = max(i_bl,i_br,i_tr,i_tl)
        # create the slice objects
        eta_rho = slice(j_min,j_max+1) # adding one as slice is non-inclusive of the end index
        xi_rho = slice(i_min,i_max+1) # adding one as slice is non-inclusive of the end index
        eta_v = slice(j_min,j_max) # not adding one to keep the correct v-grid indices
        xi_u = slice(i_min,i_max) # not adding one to keep the correct u-grid indices
        
    return eta_rho,eta_v,xi_rho,xi_u # all as slice objects

def eta_or_xi_to_slice(eta_or_xi,var_str):
    '''
    simple function used in domain_to_slice()
    '''
    # as per time, make sure we keep the eta_rho/xi dimensions after the ds.isel() step, even if we specify a single value
    # this greatly simplifies further functions for depth interpolation 
    # as we know the number of dimensions, even if some of them are single length
    # https://stackoverflow.com/questions/52190344/how-do-i-preserve-dimension-values-in-xarray-when-using-isel
    if not isinstance(eta_or_xi,slice):
        eta_or_xi = [eta_or_xi]
        if var_str=='u' or var_str=='v':
            print('rather use get_ts_uv() for extracting a time-series of u/v data')
            sys.exit()
    return eta_or_xi

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
        level=np.atleast_1d(level)
        if np.mean(level) >= 0: 
            # so we're extracting a single sigma layer
            level_for_isel = slice(level[0], level[-1]+1) 
        else:
            # sp we'll need to do vertical interpolations later 
            # for this we'll need to initially extract all the sigma levels
            level_for_isel = slice(None)
    else:
        level_for_isel = level # a slice object by definition of the logic
    
    return level_for_isel

def get_var(fname,var_str,
            grdname=None,
            tstep=slice(None),
            level=slice(None),
            eta_rho=slice(None),
            eta_v=slice(None),
            xi_rho=slice(None),
            xi_u=slice(None),
            subdomain=None,
            ref_date=None):
    '''
        extract a variable from a CROCO file
        fname = CROCO output file name (or file pattern to be used with open_mfdataset())
                fname can also be a previously extracted xarray dataset for enhanced functionality
        var_str = variable name (string) in the CROCO output file(s)
        grdname = option name of your croco grid file (only needed if the grid info is not in fname)
        tstep = time step indices to extract 
                it can be a single integer (starting at zero) or datetime
                or two values in a list e.g. [dt1,dt2], in which case the range between the two is extracted
                If slice(None), then all time-steps are extracted
        level = vertical level to extract
                If >= 0 then a sigma level is extracted 
                If <0 then a z level in meters is extracted
                If slice(None), then all sigma levels are extracted
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
        
        Retruns an xarray dataarray object of the requested data
    '''
    
    # ---------------------------
    # Get a dataset for the grid
    # ---------------------------
    if grdname is None:
        grdname = fname
    # using 'lon_rho' as var_str input to get_ds, as a hack to ensure that
    # get_ds uses open_dataset, not open_mfdataset 
    # since get_ds only uses var_str for this purpose
    ds_grd = get_ds(grdname,var_str='lon_rho') 
        
    print('extracting the data from croco file(s) - ' + var_str)
    # ----------------------------------------------
    # Prepare indices for slicing in ds.isel() below
    # ----------------------------------------------
    #
    # for each of the input dimensions we check the format of the input 
    # and construct the appropriate slice to extract using ds.isel() below
    tstep = tstep_to_slice(fname, tstep, ref_date)
    eta_rho,eta_v,xi_rho,xi_u = domain_to_slice(eta_rho,eta_v,xi_rho,xi_u,subdomain,grdname,var_str)
    level_for_isel = level_to_slice(level)

    # -------------------------
    # Get a subset of the data
    # -------------------------
    #
    if isinstance(fname, xr.Dataset) or isinstance(fname, xr.DataArray):
        ds = fname.copy()
    else:
        ds = get_ds(fname,var_str)
    ds = ds.isel(time=tstep,
                       s_rho=level_for_isel,
                       s_w=level_for_isel,
                       eta_rho=eta_rho,
                       xi_rho=xi_rho,
                       xi_u=xi_u,
                       eta_v=eta_v,
                       missing_dims='ignore' # handle case where input is a previously extracted dataset/dataarray
                       )
    ds_grd = ds_grd.isel(eta_rho=eta_rho,
                       xi_rho=xi_rho,
                       xi_u=xi_u,
                       eta_v=eta_v,
                       missing_dims='ignore' # handle case where input is a previously extracted dataset/dataarray
                       )
    if not isinstance(fname, xr.DataArray): # handle the case where input is a previously extracted dataarray (this might be the case if you want to extract a subset after doing a full extraction)
        # extract data for the requested variable
        # da is a dataarray object
        # all dimensions not related to var_str are dropped in da
        da = ds[var_str]
    else:
        da = ds.copy()
    # da = da.values # avoiding this at all costs as it's slow!!!
    
    # replace the time dimension with a list of datetimes
    if 'time' in da.dims: # handles static variables like 'h', 'angle' etc
        time_dt = get_time(fname, ref_date, time_lims=tstep)
        da = da.assign_coords(time=time_dt)
    
    # regrid u/v data onto the rho grid
    if var_str in ['u','sustr','bustr','ubar'] or var_str in ['v','svstr','bvstr','vbar']:
        if not isinstance(fname, xr.DataArray): # handle the case where input is a previously extracted dataarray (this might be the case if you want to extract a subset after doing a full extraction)
            if var_str in ['u','sustr','bustr','ubar']:
                data_rho=u2rho(da)   
            if var_str in ['v','svstr','bvstr','vbar']:
                data_rho=v2rho(da) 
        else:
            data_rho=da.copy() 
        # Create a new xarray DataArray with correct dimensions
        # now that u/v data is on the rho grid
        if len(data_rho.shape)==4:
            da_rho = xr.DataArray(data_rho, coords={'time': da['time'].values, # NB to use da not ds here!
                                                 's_rho': ds['s_rho'].values, 
                                                 'eta_rho': ds_grd['eta_rho'].values, 
                                                 'xi_rho': ds_grd['xi_rho'].values,
                                                 'lon_rho': (('eta_rho', 'xi_rho'), ds_grd['lon_rho'].values),
                                                 'lat_rho': (('eta_rho', 'xi_rho'), ds_grd['lat_rho'].values)
                                                 },
                                          dims=['time', 's_rho', 'eta_rho', 'xi_rho'])
        else: # the case where a single sigma level is extracted
            da_rho = xr.DataArray(data_rho, coords={'time': da['time'].values, # NB to use da not ds here!
                                                 'eta_rho': ds_grd['eta_rho'].values, 
                                                 'xi_rho': ds_grd['xi_rho'].values,
                                                 'lon_rho': (('eta_rho', 'xi_rho'), ds_grd['lon_rho'].values),
                                                 'lat_rho': (('eta_rho', 'xi_rho'), ds_grd['lat_rho'].values)
                                                 },
                                          dims=['time', 'eta_rho', 'xi_rho'])
        # use the same attributes
        da_rho.attrs = da.attrs
        # update da to be the data on the rho grid
        da = da_rho
    
    # ------------------------------------
    # Do vertical interpolations if needed
    # ------------------------------------
    #
    if len(da.shape)==4 and not isinstance(level,slice): # the len(da.shape)==4 check is to exclude 2D variables
        if np.mean(np.atleast_1d(level)) < 0: # we can't put this in the line above as you can't use '<' on a slice, so at least here we know 'level' is not a slice
            
            print('doing vertical interpolations - ' + var_str)
            # given the above checks in the code, here we should be dealing with a 3D variable 
            # and we want a hz slice at a constant depth level
            z=get_depths(ds) # have to use ds as we need zeta and h for this
            da_out=hlev_xarray(da, z, level)
            # use the same attributes as the original da
            da_out.attrs = da.attrs
            # update da to be the data for the specified level
            da=da_out.copy()
        
    # --------
    # Masking
    # --------
    print('applying the mask - ' + var_str)
    if isinstance(eta_rho,slice) and isinstance(xi_rho,slice) and not isinstance(fname, xr.DataArray):
            _,_,mask=get_lonlatmask(grdname,type='r', # u and v are already regridded to the rho grid so can spcify type='r' here
                                    eta_rho=eta_rho,
                                    xi_rho=xi_rho)
    else:
        mask=1
    # it looks like xarray and numpy are clever enough to use the 2D mask on a 3D or 4D variable
    # that's useful!
    da_masked=da.squeeze()*mask
    # masking throws away the attributes, so let's keep those
    da_masked.attrs = da.attrs
    da = da_masked.copy()
    
    ds.close()
    ds_grd.close()
    return da

def get_uv(fname,
           grdname=None,
           tstep=slice(None),
           level=slice(None),
           eta_rho=slice(None),
           eta_v=slice(None),
           xi_rho=slice(None),
           xi_u=slice(None),
           subdomain=None,
           ref_date=None,
           var_u='u', # could also be sustr, bustr, ubar
           var_v='v' # could also be svstr, bvstr, vbar
           ):
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
              tstep=tstep,
              level=level,
              eta_rho=eta_rho,
              eta_v=eta_v,
              xi_rho=xi_rho,
              xi_u=xi_u,
              subdomain=subdomain,
              ref_date=ref_date)
    v=get_var(fname,var_v,
              grdname=grdname,
              tstep=tstep,
              level=level,
              eta_rho=eta_rho,
              eta_v=eta_v,
              xi_rho=xi_rho,
              xi_u=xi_u,
              subdomain=subdomain,
              ref_date=ref_date)
    
    # regridding from the u and v grids to the rho grid is now handled inside 
    # get_var() which allows us to more easily do the vertical interpolation 
    # inside get_var() using the depth levels which are defined on the rho grid
    
    # -------------------
    # Rotate the vectors
    # -------------------
    
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
    u_out = u*cos_a - v*sin_a
    v_out = v*cos_a + u*sin_a
    
    # add attributes for u_out, v_out - now east,north components
    # Define a dictionary of the attributes for potential variables
    attributes = {
        'u': {
            'long_name': 'Eastward component of baroclinic velocity',
            'units': 'meters per second',
            'standard_name': 'baroclinic_eastward_sea_water_velocity'
        },
        'sustr': {
            'long_name': 'Eastward component of surface stress',
            'units': 'Newton per meter squared',
            'standard_name': 'surface_eastward_stress'
        },
        'bustr': {
            'long_name': 'Eastward component of bottom stress',
            'units': 'Newton per meter squared',
            'standard_name': 'bottom_eastward_stress'
        },
        'ubar': {
            'long_name': 'Eastward component of barotropic velocity',
            'units': 'meters per second',
            'standard_name': 'barotropic_eastward_sea_water_velocity'
        },
        'v': {
            'long_name': 'Northward component of baroclinic velocity',
            'units': 'meters per second',
            'standard_name': 'baroclinic_northward_sea_water_velocity'
        },
        'svstr': {
            'long_name': 'Northward component of surface stress',
            'units': 'Newton per meter squared',
            'standard_name': 'surface_northward_stress'
        },
        'bvstr': {
            'long_name': 'Northward component of bottom stress',
            'units': 'Newton per meter squared',
            'standard_name': 'bottom_northward_stress'
        },
        'vbar': {
            'long_name': 'Northward component of barotropic velocity',
            'units': 'meters per second',
            'standard_name': 'barotropic_northward_sea_water_velocity'
        }
    }
    
    u_out.attrs = attributes[var_u]
    v_out.attrs = attributes[var_v]
    
    # create a dataset containing both u and v
    # preferring not to do this as it makes downstream code a little easier and I'm too lazy to change it
    # ds = xr.Dataset({'u': u_out, 'v': v_out})
    
    return u_out, v_out

def get_vort(fname,
             grdname=None,
             tstep=slice(None),
             level=slice(None),
             ref_date=None):
    '''
    extract the relative vorticity from a CROCO output file:
    dv/dx - du/dy
    
    see get_var() for a description of the inputs   
    
    subsetting in space not perimitted for this. makes no sense for a single
    point, and doing it on a subset of the domain is a proper edge case.
    Actually, the subdomain input to get_var does in fact allow for you to compute
    vorticity on a subset easily... just need to implement here
    '''
    
    # start by getting u and v
    # and we'll leave them on their native grids for this calc
    # (i.e. intentionally not regridding to the rho grid)
    u=get_var(fname,'u',grdname=grdname,tstep=tstep,level=level,ref_date=ref_date)
    v=get_var(fname,'v',grdname=grdname,tstep=tstep,level=level,ref_date=ref_date)
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

def find_nearest_point(fname, Longi, Latit, Bottom=None):
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
    """
    
    if isinstance(fname, xr.Dataset) or isinstance(fname, xr.DataArray):
        ds = fname.copy()
    else:
        # for effeciency we shouldn't use open_mfdataset for this function 
        # only use the first file
        if ('*' in fname) or ('?' in fname) or ('[' in fname):
            fname=glob(fname)[0]
        ds = get_ds(fname)

    # Calculate the distance between (Longi, Latit) and all grid points
    distance = ((ds['lon_rho'].values - Longi) ** 2 +
                (ds['lat_rho'].values - Latit) ** 2) ** 0.5
    
    if Bottom is None:
        mask=ds['h'].values/ds['h'].values

    else:
        mask=ds['h'].values
        mask[mask<Bottom]=10000
        mask[mask<10000]=1
    
    distance_mask=distance*mask

    # Find the indices of the minimum distance
    # unravel_index method Converts a flat index or array of flat indices into a tuple of coordinate 
    # arrays: https://numpy.org/doc/stable/reference/generated/numpy.unravel_index.html
    min_index = np.unravel_index(distance_mask.argmin(), distance_mask.shape)

    j, i = min_index

    ds.close()

    return j, i

def get_ts_multivar(fname, lon, lat, ref_date, 
                grdname=None,
                vars = ['temp','salt'],
                i_shift=0, j_shift=0, 
                time_lims=slice(None),
                depths=slice(None),
                write_nc=False,
                fname_nc='ts.nc'
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
                        time_lims=time_lims,
                        depths=depths)
        all_datasets.append(ds_var)
    # add u,v
    ds_uv = get_ts_uv(fname, lon, lat, ref_date, 
                    grdname=grdname,
                    i_shift=i_shift, j_shift=j_shift, 
                    time_lims=time_lims,
                    depths=depths)
    all_datasets.append(ds_uv)
    # add zeta (not added to 'vars' above as we don't want to add 'depths' 
    # as input to get_ts() for obvious reasons)
    ds_zeta = get_ts(fname, 'zeta', lon, lat, ref_date, 
                    grdname=grdname,
                    i_shift=i_shift, j_shift=j_shift, 
                    time_lims=time_lims)
    all_datasets.append(ds_zeta)
    
    # merge into a single dataset
    ds_all = xr.merge(all_datasets)
    
    # write a netcdf file if specified
    if write_nc:
        ds_all.to_netcdf(fname_nc, mode='w')
    
    return ds_all

def preprocess_profile_depths(depths,default_to_bottom,h):
    # handle the different kinds of 'depths' input,
    # specifically when negative z level(s) is (are) defined
    # I'm sticking this in it's own function as we need it for both get_profile() and get_profile_uv()
    # 
    if np.mean(depths)<0:
        # we're extracting z levels
        # in which case we'll want a depth of 0 to represent the surface layer
        depths[depths==0]=-0.001 # 1mm below the surface will automatically default to the surface layer as it'll be above the top layer, so no  will get done
        # by definition a value of -99999 represents the bottom layer
        depths[depths==-99999]=0
        if default_to_bottom:
            # option to set depths deeper than the model depth to the bottom sigma layer
            # this won't actually work for an array as we now extract all z levels at once!!
            # So would need to change this
            depths[depths<-h.values] = 0
    return depths

def get_ts(fname, var, lon, lat, ref_date,
                grdname=None,
                i_shift=0, j_shift=0,
                time_lims=slice(None),
                depths=slice(None),
                default_to_bottom=False,
                write_nc=False,
                fname_nc='ts.nc',
                Bottom=None
                ):
    """
           Extract a ts from the model:
                   
            Parameters:
            - fname             :CROCO output file name (or file pattern to be used with open_mfdataset())
                                 fname can also be an xarray dataset to enhance functionality
            - var               :variable name (string) in the CROCO output file(s)
                                 (not intended for use with u,v variables - rather use get_ts_uv())
            - lat               :latitude of time-series
            - lon               :longitude of time-series
            - ref_date          :reference datetime used in croco runs
            - grdname           :optional grid file input - only needed if grid info isn't in the croco output file(s)
            - i_shift           :number of grid cells to shift along the xi axis, useful if input lon,lat is on land mask or if input depth is deeper than model depth 
            - j_shift           :number of grid cells to shift along the eta axis, (similar utility to i_shift)
            - time_lims         :time step indices to extract 
                                 it can be a single integer (starting at zero) or datetime
                                 or two values in a list e.g. [dt1,dt2], in which case the range between the two is extracted
                                 If slice(None), then all time-steps are extracted
            - depths            :if slice(None), then all sigma levels are extracted, and the depths of the levels are provided as an additional variable
                                 if a positve integer or a slice of positive integers, then those sigma levels are extracted (zero denotes the bottom layer, going upward to the surface)
                                 if a negative number, or a list of negative numers then data interpolated to those z levels are extracted 
                                 if zero is contained in a list of negative numbers, then it is treated as the surface layer
                                 if -99999 is contained in a list of negative numbers, it is assumed to represent the bottom layer in the model 
            - default_to_bottom :flag to extract data for the bottom layer if input z level is below the model seafloor (True/False)
            - write_nc          :write a netcdf file? (True/False)
            - fname_nc          :netcdf file name. Only used if write_nc = True
            - Bottom (positive value): if the model bathy is slightly different. This Option to find nearest
              lat and lon in water that is as deep as reference. If == None then
              this looks for only the closest horizontal point.
            
            Returns:
            - ds, an xarray dataset containing the ts data
    """
    if var=='u' or var=='v':
        print('rather use get_ts_uv() for extracting a time-series/ profile of u/v data')
        sys.exit()
    
    #find_nearest_point finds the nearest point in the model to the model grid lon, lat extracted from the model grid input.
    if grdname is None:
        grdname = fname
    j, i = find_nearest_point(grdname, lon, lat, Bottom) 
    
    # apply the shifts along the xi and eta axis
    i = i+i_shift
    j = j+j_shift
    
    # get the model depth and grid
    h = get_grd_var(grdname,"h",eta_rho=j,xi_rho=i)
    
    if not isinstance(depths,slice):
        # we're extracting data at specified z level(s) or a single sigma level
        depths=np.atleast_1d(depths).astype('float32') # makes life easier for handling both profiles and time-series if they're both arrays
        depths=preprocess_profile_depths(depths,default_to_bottom,h)
        
    ts_da = get_var(fname, var,
                          grdname=grdname,
                          tstep=time_lims,
                          level=depths,
                          eta_rho=j,
                          xi_rho=i,
                          ref_date=ref_date)
    
    if 's_rho' in ts_da.coords:
        # in this case we want to include the depths of the sigma levels in the output
        # we use the get_depths() function, which takes the dataset as input        
        # so we need to extract the dataset again here unfortunately
        if isinstance(fname, xr.Dataset) or isinstance(fname, xr.DataArray):
            ds = fname.copy()
        else:
            ds = get_ds(fname,var)
        ds = ds.isel(time=time_lims,
                            s_rho=depths,
                            eta_rho=slice(j,j+1), # making it a slice to maintain the spatial dimensions for input to get_depths()
                            xi_rho=slice(i,i+1))
        
        depths_da = get_depths(ds)
    
        # create a new dataset with the extracted profile and depths of the sigma levels at this grid cell 
        ds = xr.Dataset({var: ts_da, 'depth': depths_da, 'h': h})
    else:
        ds = xr.Dataset({var: ts_da, 'h': h})
    
    # write a netcdf file if specified
    if write_nc:
        ds.to_netcdf(fname_nc)
    
    return ds

def get_ts_uv(fname, lon, lat, ref_date, 
                grdname=None,
                i_shift=0, j_shift=0, 
                time_lims=slice(None),
                depths=slice(None),
                default_to_bottom=False,
                write_nc=False,
                fname_nc='ts_uv.nc',
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
    
    # get the model depth at this location
    h = get_grd_var(grdname,"h",eta_rho=j,xi_rho=i)
    
    #-------
    if not isinstance(depths,slice):
        # we're extracting data at specified z level(s) or a single sigma level
        depths=np.atleast_1d(depths).astype('float32') # makes life easier for handling both profiles and time-series if they're both arrays
        depths=preprocess_profile_depths(depths,default_to_bottom,h)
        
    u_ts_da,v_ts_da = get_uv(fname,
                          grdname=grdname,
                          tstep=time_lims,
                          level=depths,
                          eta_rho=j_rho,
                          xi_rho=i_rho,
                          eta_v=j_v,
                          xi_u=i_u,
                          ref_date=ref_date)
    
    # pull out the middle data point from our 3x3 block of rho grid points
    # this is by definition the grid cell we are interested in
    u_ts_da = u_ts_da[:,:,1,1]
    v_ts_da = v_ts_da[:,:,1,1]
    
    if 's_rho' in u_ts_da.coords:
        # in this case we want to include the depths of the sigma levels in the output
        # we use the get_depths() function, which takes the dataset as input        
        # so we need to extract the dataset again here unfortunately
        if isinstance(fname, xr.Dataset) or isinstance(fname, xr.DataArray):
            ds = fname.copy()
        else:
            ds = get_ds(fname)
        ds = ds.isel(time=time_lims,
                            s_rho=depths,
                            eta_rho=slice(j,j+1), # making it a slice to maintain the spatial dimensions for input to get_depths()
                            xi_rho=slice(i,i+1))
        
        depths_da = get_depths(ds)
    
        # create a new dataset with the extracted profile and depths of the sigma levels at this grid cell 
        ds = xr.Dataset({'u': u_ts_da,'v': v_ts_da, 'depth': depths_da, 'h': h})
    else:
        ds = xr.Dataset({'u': u_ts_da,'v': v_ts_da, 'h': h})
    
    # write a netcdf file if specified
    if write_nc:
        ds.to_netcdf(fname_nc)
    
    return ds


'''
An idea for get_section, which I think will be super fast
get_section can do the hz interpolation at every point by mapping the lon, lat inputs to 
eta_rho, xi_rho indices (or rather decimal indices)
then you can simply use xarrays interp. function to extract the section
quick chatgpt search suggests from scipy.interpolate import RegularGridInterpolator
We essentially use the regular grid of eta_rho, xi_rho, with lat_rho,lon_rho defined on that grid
Then the inputs lons/lats get corresponding eta_rho,xi_rho decimal indices
which get input to xarrays .interp functionality, after we've extracted our data 
over a subdomain with get_var

'''
