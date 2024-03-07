import numpy as np
from datetime import timedelta
import xarray as xr
import dask
from datetime import timedelta, datetime
from glob import glob
import sys

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

def hlev(var, z, depth):
    """
    this extracts a horizontal slice
    
    (TODO: DEFINITELY SCOPE TO IMPROVE EFFICIENCY
     AT THE MOMENT WE LOOP THROUGH ALL eta,xi INDICES 
     AND INTERPOLATE ON EACH - YIKES!)

    var = 3D extracted variable of interest (assuming mask is already nan - use get_var() method in this file)
    z = depths (in m) of sigma levels, also 3D array (use get_depths() method in this file)
    depth = the horizontal depth you want (should be negative)

    Adapted (by J.Veitch and G.Fearon) from vinterp.m in roms_tools (by P.Penven)
    """
    [N, Mp, Lp] = np.shape(z)

    a = z.copy()
    a[a >= depth] = 0
    a[a < depth] = 1
    
    # 2D variable containing the sigma level index directly above the constant depth level
    levs = np.sum(a, axis=0)
    # values of zero indicate the depth level is below our sigma levels so set to nan
    levs[levs==0] = np.nan
    # make levs the sigma level index directly below the constant depth level
    levs = levs -1

    vnew = np.zeros((Mp, Lp))
    vnew[np.isnan(levs)]=np.nan
    
    # looping through every horizontal grid point makes this slow
    for m in np.arange(Mp):
        for l in np.arange(Lp):
            
            if not np.isnan(levs[m,l]):
            
                ind1 = int(levs[m, l])
                ind2 = int(levs[m, l]) + 1
                
                if ind1 == N-1: # there is no level above to interpolate between
                    # so I'd rather use the surface layer than extrapolate
                    vnew[m, l] = var[ind1, m, l]
                else:
    
                    v1 = var[ind1, m, l]
                    v2 = var[ind2, m, l]
        
                    z1 = z[ind1, m, l]
                    z2 = z[ind2, m, l]
        
                    vnew[m, l] = ((v1 - v2) * depth + v2 * z1 - v1 * z2) / (z1 - z2)

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
    
    return depth_rho

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
        ref_date = reference date for the croco run as a datetime object
        time_lims = optional list of two datetimes i.e. [dt1,dt2], which define the range of times to extract
                    If slice(None), then all time-steps are extracted
    '''
    ds = get_ds(fname)
    time = ds.time.values
    
    if ref_date is None:
        print('ref_date is not defined - using default of 2000-01-01')
        ref_date=datetime(2000,1,1)

    # convert 'time' (in seconds since ref_date) to a list of datetimes
    time_dt = []
    for t in time:
        date_now = ref_date + timedelta(seconds=np.float64(t))
        time_dt.append(date_now)
    
    # subset based in time_lims input
    if not isinstance(time_lims,slice):
        time_lims = tstep_to_slice(fname, time_lims, ref_date)
        
    time_dt = time_dt[time_lims]
    
    ds.close()
    return time_dt
    
def get_lonlatmask(fname,type='r',
                   eta_rho=slice(None),
                   xi_rho=slice(None)):
    
    # for effeciency we shouldn't use open_mfdataset for this function 
    # only use the first file
    if ('*' in fname) or ('?' in fname) or ('[' in fname):
        fname=glob(fname)[0]
    
    ds = get_ds(fname)
    ds = ds.isel(eta_rho=eta_rho,xi_rho=xi_rho)
    
    lon = ds.lon_rho.values 
    lat = ds.lat_rho.values 
    mask = ds.mask_rho.values
    
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
    Take the input to get_var, and return a slice object to be used to
    subset the dataset using ds.isel()
    see get_var() for how this is used
    '''
    # check if tstep input is instance of datetime, 
    # in which case convert it/them into the correct time index/indices
    if isinstance(np.atleast_1d(tstep)[0],datetime):
        if ref_date is None:
            print('ref_date is not defined - using default of 2000-01-01')
            ref_date=datetime(2000,1,1)
        time_croco = get_time(fname,ref_date) # get_time actually calls tstep_to_slice (CIRCULAR!), but only inside an if statement which won't be entered no time_lims input. MESSY. Should do better
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

def eta_or_xi_to_slice(eta_or_xi,var_str):
    '''
    Take the input to get_var, and return a slice object to be used to
    subset the dataset using ds.isel()
    see get_var() for how this is used
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
        # so level is a single number
        if level >= 0: 
            # so we're extracting a single sigma layer
            level_for_isel = slice(level, level+1) # this gets a single level 
        else:
            # sp we'll need to do vertical interpolations later 
            # for this we'll need to initially extract all the sigma levels
            level_for_isel = slice(None)
    else:
        level_for_isel = level # a slice object by definition of the logic
    
    return level_for_isel

def get_var(fname,var_str,
            tstep=slice(None),
            level=slice(None),
            eta_rho=slice(None),
            eta_v=slice(None),
            xi_rho=slice(None),
            xi_u=slice(None),
            ref_date=None):
    '''
        extract a variable from a CROCO file
        fname = CROCO output file name (or file pattern to be used with open_mfdataset())
        var_str = variable name (string) in the CROCO output file(s)
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
        ref_date = reference datetime used in croco runs
        
        Retruns an xarray dataarray object of the requested data
    '''
    
    print('extracting the data from croco file(s) - ' + var_str)
    # ----------------------------------------------
    # Prepare indices for slicing in ds.isel() below
    # ----------------------------------------------
    #
    # for each of the input dimensions we check the format of the input 
    # and construct the appropriate slice to extract using ds.isel() below
    tstep = tstep_to_slice(fname, tstep, ref_date)
    eta_rho = eta_or_xi_to_slice(eta_rho,var_str)
    xi_rho = eta_or_xi_to_slice(xi_rho,var_str)
    level_for_isel = level_to_slice(level)

    # -------------------------
    # Get a subset of the data
    # -------------------------
    #
    ds = get_ds(fname,var_str)
    ds = ds.isel(time=tstep,
                       s_rho=level_for_isel,
                       s_w=level_for_isel,
                       eta_rho=eta_rho,
                       xi_rho=xi_rho,
                       xi_u=xi_u,
                       eta_v=eta_v
                       )
    # extract data for the requested variable
    # da is a dataarray object
    # all dimensions not related to var_str are dropped in da
    da = ds[var_str]
    # da = da.values # avoiding this at all costs as it's slow!!!
    
    # replace the time dimension with a list of datetimes
    if 'time' in da.dims: # handles static variables like 'h', 'angle' etc
        time_dt = get_time(fname, ref_date, time_lims=tstep)
        da = da.assign_coords(time=time_dt)
    
    # regrid u/v data onto the rho grid
    if var_str == 'u' or var_str == 'v':
        if var_str=='u':
            data_rho=u2rho(da)   
        if var_str=='v':
            data_rho=v2rho(da) 
        # Create a new xarray DataArray with correct dimensions
        # now that u/v data is on the rho grid
        if len(data_rho.shape)==4:
            da_rho = xr.DataArray(data_rho, coords={'time': da['time'].values, # NB to use da not ds here!
                                                 's_rho': ds['s_rho'].values, 
                                                 'eta_rho': ds['eta_rho'].values, 
                                                 'xi_rho': ds['xi_rho'].values,
                                                 'lon_rho': (('eta_rho', 'xi_rho'), ds['lon_rho'].values),
                                                 'lat_rho': (('eta_rho', 'xi_rho'), ds['lat_rho'].values)
                                                 },
                                          dims=['time', 's_rho', 'eta_rho', 'xi_rho'])
        else: # the case where a single sigma level is extracted
            da_rho = xr.DataArray(data_rho, coords={'time': da['time'].values, # NB to use da not ds here!
                                                 'eta_rho': ds['eta_rho'].values, 
                                                 'xi_rho': ds['xi_rho'].values,
                                                 'lon_rho': (('eta_rho', 'xi_rho'), ds['lon_rho'].values),
                                                 'lat_rho': (('eta_rho', 'xi_rho'), ds['lat_rho'].values)
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
        if level < 0: # we can't put this in the line above as you can't use '<' on a slice, so at least here we know 'level' is not a slice
            print('getting numpy array of data for doing vertical interpolations - ' + var_str)
            # unavoidable here unfortunately...
            # this step can be a bit slow, but if we don't do this the vertical interpolation becomes snail pace
            # TODO: we need a more efficient way of doing the vertical interpolation
            data = da.values
            print('doing vertical interpolations - ' + var_str)
            # given the above checks in the code, here we should be dealing with a 3D variable 
            # and we want a hz slice at a constant depth level
            z=get_depths(ds) # have to use ds as we need zeta and h for this
            # z is on the rho grid, but u and v are already regridded to the rho grid above so no need to catch anything here
            T,D,M,L=da.shape
            data_out=np.zeros((T,M,L))
            for t in np.arange(T):
                data_out[t,:,:]=hlev(data[t,::], z[t,::], level)
            # create a new dataarray for the data for this level
            da_out = xr.DataArray(data_out, coords={'time': da['time'].values, # NB to use da not ds here!
                                                 'eta_rho': ds['eta_rho'].values, 
                                                 'xi_rho': ds['xi_rho'].values,
                                                 'lon_rho': (('eta_rho', 'xi_rho'), ds['lon_rho'].values),
                                                 'lat_rho': (('eta_rho', 'xi_rho'), ds['lat_rho'].values)
                                                 },
                                          dims=['time', 'eta_rho', 'xi_rho'])
            # use the same attributes
            da_out.attrs = da.attrs
            # update da to be the data for the specified level
            da=da_out.copy()
        
    # --------
    # Masking
    # --------
    print('applying the mask - ' + var_str)
    if isinstance(eta_rho,slice) and isinstance(xi_rho,slice):
        _,_,mask=get_lonlatmask(fname,type='r', # u and v are already regridded to the rho grid so can spcify type='r' here
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
    
    return da

def get_uv(fname,
           tstep=slice(None),
           level=slice(None),
           eta_rho=slice(None),
           eta_v=slice(None),
           xi_rho=slice(None),
           xi_u=slice(None),
           ref_date=None):
    '''
    extract u and v components from a CROCO output file(s), regrid onto the 
    rho grid and rotate from grid-aligned to east-north components
    
    see get_var() for a description of the inputs  
    
    returns xarray dataarrays for both u and v data
    
    '''
    
    u=get_var(fname,'u',
              tstep=tstep,
              level=level,
              eta_rho=eta_rho,
              eta_v=eta_v,
              xi_rho=xi_rho,
              xi_u=xi_u,
              ref_date=ref_date)
    v=get_var(fname,'v',
              tstep=tstep,
              level=level,
              eta_rho=eta_rho,
              eta_v=eta_v,
              xi_rho=xi_rho,
              xi_u=xi_u,
              ref_date=ref_date)
    
    # regridding from the u and v grids to the rho grid is now handled inside 
    # get_var() which allows us to more easily do the vertical interpolation 
    # inside get_var() using the depth levels which are defined on the rho grid
    
    # -------------------
    # Rotate the vectors
    # -------------------
    
    # grid angle
    angle=get_var(fname, 'angle',
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
    # u
    u_out.attrs['long_name'] = 'Eastward velocity'
    u_out.attrs['units'] = 'meters per second'
    u_out.attrs['standard_name'] = 'eastward_sea_water_velocity'
    # v
    v_out.attrs['long_name'] = 'Northward velocity'
    v_out.attrs['units'] = 'meters per second'
    v_out.attrs['standard_name'] = 'northward_sea_water_velocity'
    
    # create a dataset containing both u and v
    # preferring not to do this as it makes downstream code a little easier and I'm too lazy to change it
    # ds = xr.Dataset({'u': u_out, 'v': v_out})
    
    return u_out, v_out

def get_vort(fname,
             tstep=slice(None),
             level=slice(None),
             ref_date=None):
    '''
    extract the relative vorticity from a CROCO output file:
    dv/dx - du/dy
    
    see get_var() for a description of the inputs   
    
    subsetting in space not perimitted for this. makes no sense for a single
    point, and doing it on a subset of the domain is a proper edge case
    '''
    
    # start by getting u and v
    # and we'll leave them on their native grids for this calc
    # (i.e. intentionally not regridding to the rho grid)
    u=get_var(fname,'u',tstep=tstep,level=level,ref_date=ref_date)
    v=get_var(fname,'v',tstep=tstep,level=level,ref_date=ref_date)
    pm=get_var(fname, 'pm') # 1/dx on the rho grid
    pn=get_var(fname, 'pn') # 1/dy on the rho grid
    
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
    
    return vort

def get_boundary(fname):
    '''
    Return lon,lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    lon_rho,lat_rho,_=get_lonlatmask(fname,type='r')
    lon = np.hstack((lon_rho[0:, 0], lon_rho[-1, 1:-1],
                     lon_rho[-1::-1, -1], lon_rho[0, -2::-1]))
    lat = np.hstack((lat_rho[0:, 0], lat_rho[-1, 1:-1],
                     lat_rho[-1::-1, -1], lat_rho[0, -2::-1]))
    return lon, lat

def find_nearest_point(fname, Longi, Latit):
    """
            Find the nearest indices of the model rho grid to a specified lon, lat coordinate:
            
            Parameters:
            - fname :filename of the model
            - Longi :longitude
            - Latit :latitude

            Returns:
            - j :the nearest eta index
            - i :the nearest xi index
    """

    # for effeciency we shouldn't use open_mfdataset for this function  
    # only use the first file if fname is specified as a pattern of file names
    if ('*' in fname) or ('?' in fname) or ('[' in fname):
        fname = glob(fname)[0]

    with get_ds(fname) as ds:
        # Calculate the distance between (Longi, Latit) and all grid points
        distance = ((ds['lon_rho'].values - Longi) ** 2 +
                    (ds['lat_rho'].values - Latit) ** 2) ** 0.5

    # Find the indices of the minimum distance
    # unravel_index method Converts a flat index or array of flat indices into a tuple of coordinate 
    # arrays: https://numpy.org/doc/stable/reference/generated/numpy.unravel_index.html
    min_index = np.unravel_index(distance.argmin(), distance.shape)

    j, i = min_index

    return j, i

def get_ts_OLD(fname, var, lon, lat, ref_date, depth=-1, 
           i_shift=0, j_shift=0, 
           time_lims=slice(None),
           write_nc=False,
           fname_nc='ts.nc'
           ):
    """
           This function has been made redundant for the new get_ts() which also handles profiles            
           I'm just keeping it here in case we want to check anything from the old function
           Ultimately this function will be deleted once we're comfortable we don't need it           
           
           Extract a timeseries from the model:
                   
            Parameters:
            - fname             :CROCO output file name (or file pattern to be used with open_mfdataset())
            - var               :variable name (string) in the CROCO output file(s)
                                 (not intended for use with u,v variables - rather use get_ts_uv())
            - lat               :latitude of time-series
            - lon               :longitude of time-series
            - ref_date          :reference datetime used in croco runs
            - depth             :vertical level to extract
                                 If >= 0 then a sigma level is extracted 
                                 If <0 then a z level in meters is extracted
            - i_shift         :number of grid cells to shift along the xi axis, useful if input lon,lat is on land mask or if input depth is deeper than model depth 
            - j_shift         :number of grid cells to shift along the eta axis, (similar utility to i_shift)
            - time_lims         :time step indices to extract 
                                 two values in a list e.g. [dt1,dt2], in which case the range between the two is extracted
                                 If slice(None), then all time-steps are extracted
            - write_nc          :write a netcdf file? (True/False)
            - fname_nc          :netcdf file name. Only used if write_nc = True
            
            Returns:
            - ds, an xarray dataset containing the extracted time-series as well as 'h', the model depth at the model grid cell
    """
    
    if var=='u' or var=='v':
        print('rather use get_ts_uv() for extracting a time-series of u/v data')
        sys.exit()
    
    # finds the rho grid indices nearest to the input lon, lat
    j, i = find_nearest_point(fname, lon, lat) 
    
    # apply the shifts along the xi and eta axis
    i = i+i_shift
    j = j+j_shift
    
    # get the data from the model
    # But first check the model depth against the input depth
    h = get_var(fname,"h",eta_rho=j,xi_rho=i)
    if h<-depth:
        print('The model depth is shallower than input/(in situ) depth!!')
        print('Were rather extracting the bottom sigma layer of the model.')
        da = get_var(fname, var,
                             tstep=time_lims,
                             level=0, # extracts bottom layer
                             eta_rho=j,
                             xi_rho=i,
                             ref_date=ref_date)
    else:
        da = get_var(fname, var,
                             tstep=time_lims,
                             level=depth,
                             eta_rho=j,
                             xi_rho=i,
                             ref_date=ref_date)
    
    # create a new dataset with the extracted time-series, and add the model depth at this grid cell 
    ds = xr.Dataset({var: da, 'h': h})
    
    # write a netcdf file if specified
    if write_nc:
        ds.to_netcdf(fname_nc)
    
    return ds

def get_ts_uv_OLD(fname, lon, lat, ref_date, depth=-1, 
              i_shift=0, j_shift=0, 
              time_lims=slice(None),
              write_nc=False,
              fname_nc='ts_uv.nc'
              ):
    """
            This function has been made redundant for the new get_ts_uv() which also handles profiles            
            I'm just keeping it here in case we want to check anything from the old function
            Ultimately this function will be deleted once we're comfortable we don't need it 
            
            Extract a timeseries of u,v from the model:
                   
            Parameters:
            - fname             :CROCO output file name (or file pattern to be used with open_mfdataset())
            - lat               :latitude of time-series
            - lon               :longitude of time-series
            - ref_date          :reference datetime used in croco runs
            - depth             :vertical level to extract
                                 If >= 0 then a sigma level is extracted 
                                 If <0 then a z level in meters is extracted
            - i_shift         :number of grid cells to shift along the xi axis, useful if input lon,lat is on land mask or if input depth is deeper than model depth 
            - j_shift         :number of grid cells to shift along the eta axis, (similar utility to i_shift)
            - time_lims         :time step indices to extract 
                                 two values in a list e.g. [dt1,dt2], in which case the range between the two is extracted
                                 If slice(None), then all time-steps are extracted
            - write_nc          :write a netcdf file? (True/False)
            - fname_nc          :netcdf file name. Only used if write_nc = True
                     
            Returns:
            - ds, an xarray dataset containing the extracted time-series as well as 'h', the model depth at the model grid cell
              (u,v variables are rotated from grid-aligned to east-north components)
    """
    
    # finds the rho grid indices nearest to the input lon, lat
    j, i = find_nearest_point(fname, lon, lat) 
    
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
    
    # get the data from the model
    # But first check the model depth against the input depth
    h = get_var(fname,"h",eta_rho=j,xi_rho=i)
    if h<-depth:
        print('The model depth is shallower than input/(in situ) depth!!')
        print('Were rather extracting the bottom sigma layer of the model.')
        u,v = get_uv(fname,
                              tstep=time_lims,
                              level=0, # extracts bottom layer
                              eta_rho=j_rho,
                              xi_rho=i_rho,
                              eta_v=j_v,
                              xi_u=i_u,
                              ref_date=ref_date)
    else:
        u,v = get_uv(fname,
                             tstep=time_lims,
                             level=depth,
                             eta_rho=j_rho,
                             xi_rho=i_rho,
                             eta_v=j_v,
                             xi_u=i_u,
                             ref_date=ref_date)
    
    # now that u and v on our 3x3 subset of the rho grid, 
    # let's pull out the time-series from the middle grid cell
    u = u[:,1,1]
    v = v[:,1,1]
    
    # create a new dataset with the extracted time-series, and add the model depth at this grid cell 
    ds = xr.Dataset({'u': u, 'v': v, 'h': h})
    
    # write a netcdf file if specified
    if write_nc:
        ds.to_netcdf(fname_nc)
    
    return ds

def get_ts_multivar(fname, lon, lat, ref_date, 
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
                        i_shift=i_shift, j_shift=j_shift, 
                        time_lims=time_lims,
                        depths=depths)
        all_datasets.append(ds_var)
    # add u,v
    ds_uv = get_ts_uv(fname, lon, lat, ref_date, 
                    i_shift=i_shift, j_shift=j_shift, 
                    time_lims=time_lims,
                    depths=depths)
    all_datasets.append(ds_uv)
    # add zeta (not added to 'vars' above as we don't want to add 'depths' 
    # as input to get_ts() for obvious reasons)
    ds_zeta = get_ts(fname, 'zeta', lon, lat, ref_date, 
                    i_shift=i_shift, j_shift=j_shift, 
                    time_lims=time_lims)
    all_datasets.append(ds_zeta)
    
    # merge into a single dataset
    ds_all = xr.merge(all_datasets)
    
    # write a netcdf file if specified
    if write_nc:
        ds_all.to_netcdf(fname_nc)
    
    return ds_all

def preprocess_profile_depths(depths,default_to_bottom,h):
    # handle the different kinds of 'depths' input,
    # specifically when negative z level(s) is (are) defined
    # I'm sticking this in it's own function as we need it for both get_profile() and get_profile_uv()
    # 
    depths=np.atleast_1d(depths) # makes life easier for handling both profiles and time-series
    if np.mean(depths)<0:
        # we're extracting z levels
        depths=depths.astype('float64')
        # in which case we'll want a depth of 0 to represent the surface layer
        depths[depths==0]=-0.001 # 1mm below the surface will automatically default to the surface layer as it'll be above the top layer, so no  will get done
        # by definition a value of -99999 represents the bottom layer
        depths[depths==-99999]=0
        if default_to_bottom:
            # option to set depths deeper than the model depth to the bottom sigma layer
            depths[depths<-h.values] = 0
    return depths

def get_ts(fname, var, lon, lat, ref_date, 
                i_shift=0, j_shift=0, 
                time_lims=slice(None),
                depths=slice(None),
                default_to_bottom=False,
                write_nc=False,
                fname_nc='ts.nc'
                ):
    """
           Extract a ts from the model:
                   
            Parameters:
            - fname             :CROCO output file name (or file pattern to be used with open_mfdataset())
            - var               :variable name (string) in the CROCO output file(s)
                                 (not intended for use with u,v variables - rather use get_ts_uv())
            - lat               :latitude of time-series
            - lon               :longitude of time-series
            - ref_date          :reference datetime used in croco runs
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
            
            Returns:
            - ds, an xarray dataset containing the ts data
    """
    if var=='u' or var=='v':
        print('rather use get_ts_uv() for extracting a time-series/ profile of u/v data')
        sys.exit()
        
    # get the model time
    time_model = get_time(fname, ref_date, time_lims=time_lims)
    time_lims = tstep_to_slice(fname, time_lims, ref_date)
    
    #find_nearest_point finds the nearest point in the model to the model grid lon, lat extracted from the model grid input.
    j, i = find_nearest_point(fname, lon, lat) 
    
    # apply the shifts along the xi and eta axis
    i = i+i_shift
    j = j+j_shift
    
    # get the model depth at this location
    h = get_var(fname,"h",eta_rho=j,xi_rho=i)
    
    if isinstance(depths,slice):
        # so we're extracting some or all of the sigma levels
        # (note the case of extracting a single sigma level is handled in the else: below)
        ts_da = get_var(fname, var,
                              tstep=time_lims,
                              level=depths,
                              eta_rho=j,
                              xi_rho=i,
                              ref_date=ref_date)
        
        if 's_rho' in ts_da.coords:
            # we explicitly check for existance of 's_rho' so we can handle 2D variables
            
            # get the depths of the sigma levels using the get_depths() function, which takes the dataset as input
            ds = get_ds(fname, var)
            ds = ds.isel(time=time_lims,
                                s_rho=depths,
                                eta_rho=slice(j,j+1), # making it a slice to maintain the spacial dimensions for input to get_depths()
                                xi_rho=slice(i,i+1))
            depths_out = np.squeeze(get_depths(ds))
            
            # we need to create a new dataarray for the depths of the sigma levels
            depths_da = xr.DataArray(depths_out, coords={'time': ts_da['time'].values,
                                                 'eta_rho': ts_da['eta_rho'].values, 
                                                 'xi_rho': ts_da['xi_rho'].values,
                                                 's_rho': ts_da['s_rho'].values,
                                                 'lon_rho': ts_da['lon_rho'].values,
                                                 'lat_rho': ts_da['lat_rho'].values
                                                 },
                                          dims=['time', 's_rho'])
            depths_da.attrs['long_name'] = 'Depth of sigma levels of the rho grid (centred in grid cells)'
            depths_da.attrs['units'] = 'meter'
            depths_da.attrs['positive'] = 'up'
        
            # create a new dataset with the extracted profile and depths of the sigma levels at this grid cell 
            ds = xr.Dataset({var: ts_da, 'depth': depths_da, 'h': h})
        else:
            # handle the case of a 2D variable like 'zeta', where 'depths' won't be (or shouldn't be) a specified input
            ds = xr.Dataset({var: ts_da, 'h': h})
            
    else:
        # we're extracting data at specified z level(s) or a single sigma level
        depths_in=depths # keep a record of the user input for writing to the dataarray
        depths=preprocess_profile_depths(depths,default_to_bottom,h)
        # set up an empty array of the correct size, which we populate in a loop through depths
        ts = np.zeros((len(time_model),len(depths)))
        for index, depth in enumerate(depths):
            # extract a time-series for this z level (or sigma level)
            if depth>=0:
                depth=int(depth) # must be an integer for a sigma level
            ts_i = get_var(fname, var,
                              tstep=time_lims,
                              level=depth,
                              eta_rho=j,
                              xi_rho=i,
                              ref_date=ref_date)
            # and populate the ts array
            ts[:,index]=ts_i.values
        
        # create an xarray dataarray for ts
        # We can use the last 'ts_i' dataarray to extract some of the coordinate information
        if len(depths) > 1:
            # this is a profile extraction so we need a 'depth' dimension
            ts_da = xr.DataArray(ts, coords={'time': ts_i['time'].values,
                                                 'eta_rho': ts_i['eta_rho'].values, 
                                                 'xi_rho': ts_i['xi_rho'].values,
                                                 'depth': depths_in,
                                                 'lon_rho': ts_i['lon_rho'].values,
                                                 'lat_rho': ts_i['lat_rho'].values
                                                 },
                                          dims=['time', 'depth'])
            ts_da.coords['depth'].attrs['long_name'] = 'water depth from free surface'
            ts_da.coords['depth'].attrs['units'] = 'meters'
            ts_da.coords['depth'].attrs['postive'] = 'up'
        else: 
            # we're extracting a time-series (either a sigma level or a z level), so drop the 'depth' dimension
            ts_da = xr.DataArray(ts.squeeze(), coords={'time': ts_i['time'].values,
                                                 'eta_rho': ts_i['eta_rho'].values, 
                                                 'xi_rho': ts_i['xi_rho'].values,
                                                 'lon_rho': ts_i['lon_rho'].values,
                                                 'lat_rho': ts_i['lat_rho'].values
                                                 },
                                          dims=['time'])
        ts_da.attrs = ts_i.attrs # add the attributes for this variable
        
        # create a new dataset with the extracted time-series at this grid cell 
        ds = xr.Dataset({var: ts_da, 'h': h})
    
    # write a netcdf file if specified
    if write_nc:
        ds.to_netcdf(fname_nc)
    
    return ds

def get_ts_uv(fname, lon, lat, ref_date, 
                i_shift=0, j_shift=0, 
                time_lims=slice(None),
                depths=slice(None),
                default_to_bottom=False,
                write_nc=False,
                fname_nc='ts_uv.nc'
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
    
    # get the model time
    time_model = get_time(fname, ref_date, time_lims=time_lims)
    time_lims = tstep_to_slice(fname, time_lims, ref_date)
    
    # finds the rho grid indices nearest to the input lon, lat
    j, i = find_nearest_point(fname, lon, lat) 
    
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
    h = get_var(fname,"h",eta_rho=j,xi_rho=i)
    
    if isinstance(depths,slice):
        # so we're extracting some or all of the sigma levels
        # (note the case of extracting a single sigma level is actually handled in the else: below)
        u_ts_da,v_ts_da = get_uv(fname,
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
        
        # note here we don't check to see if 's_rho' exists like we did in get_ts(), 
        # as we are dealing with u,v which by definition have the 's_rho' coordinate
        
        # get the depths of the sigma levels using the get_depths() function, which takes the dataset as input
        ds = get_ds(fname)
        ds = ds.isel(time=time_lims,
                            s_rho=depths,
                            eta_rho=slice(j,j+1), # making it a slice to maintain the spacial dimensions for input to get_depths()
                            xi_rho=slice(i,i+1))
        depths_out = np.squeeze(get_depths(ds))
        
        # we need to create a new dataarray for the depths of the sigma levels
        depths_da = xr.DataArray(depths_out, coords={'time': u_ts_da['time'].values,
                                             'eta_rho': u_ts_da['eta_rho'].values, 
                                             'xi_rho': u_ts_da['xi_rho'].values,
                                             's_rho': u_ts_da['s_rho'].values,
                                            'lon_rho': u_ts_da['lon_rho'].values,
                                            'lat_rho': u_ts_da['lat_rho'].values
                                             },
                                      dims=['time', 's_rho'])
        depths_da.attrs['long_name'] = 'Depth of sigma levels of the rho grid (centred in grid cells)'
        depths_da.attrs['units'] = 'meter'
        depths_da.attrs['positive'] = 'up'
        
        # create a new dataset with the extracted profile and depths of the sigma levels at this grid cell 
        ds = xr.Dataset({'u': u_ts_da,'v': v_ts_da, 'depth': depths_da, 'h': h})
    
    else:
        # we're extracting data at specified z level(s) or a single sigma level
        depths_in=depths # keep a record of the user input for writing to the dataarray
        depths=preprocess_profile_depths(depths,default_to_bottom,h)
        # set up empty arrays of the correct size, which we populate in a loop through depths
        u_ts = np.zeros((len(time_model),len(depths)))
        v_ts = np.zeros((len(time_model),len(depths)))
        for index, depth in enumerate(depths):
            # extract a time-series for this z level (or sigma level)
            if depth>=0: 
                depth=int(depth) # must be an integer for a sigma level
            u_ts_i,v_ts_i = get_uv(fname,
                                 tstep=time_lims,
                                 level=depth,
                                 eta_rho=j_rho,
                                 xi_rho=i_rho,
                                 eta_v=j_v,
                                 xi_u=i_u,
                                 ref_date=ref_date)
            
            # pull out the middle data point from our 3x3 block of rho grid points
            # this is by definition the grid cell we are interested in
            u_ts_i = u_ts_i[:, 1, 1] 
            v_ts_i = v_ts_i[:, 1, 1]
            # populate our output arrays
            u_ts[:, index] = u_ts_i.values
            v_ts[:, index] = v_ts_i.values
            
        # create xarray dataarrays for u_ts and v_ts
        if len(depths) > 1:
            # this is a profile extraction so we need a 'depth' dimension
            u_ts_da = xr.DataArray(u_ts, coords={'time': u_ts_i['time'].values,
                                                 'eta_rho': u_ts_i['eta_rho'].values, 
                                                 'xi_rho': u_ts_i['xi_rho'].values,
                                                 'depth': depths_in,
                                                 'lon_rho': u_ts_i['lon_rho'].values,
                                                 'lat_rho': u_ts_i['lat_rho'].values
                                                 },
                                          dims=['time', 'depth'])
            u_ts_da.coords['depth'].attrs['long_name'] = 'water depth from free surface'
            u_ts_da.coords['depth'].attrs['units'] = 'meters'
            u_ts_da.coords['depth'].attrs['postive'] = 'up'
            
            v_ts_da = xr.DataArray(v_ts, coords={'time': v_ts_i['time'].values,
                                                 'eta_rho': v_ts_i['eta_rho'].values, 
                                                 'xi_rho': v_ts_i['xi_rho'].values,
                                                 'depth': depths_in,
                                                 'lon_rho': v_ts_i['lon_rho'].values,
                                                 'lat_rho': v_ts_i['lat_rho'].values
                                                 },
                                          dims=['time', 'depth'])
            v_ts_da.coords['depth'].attrs['long_name'] = 'water depth from free surface'
            v_ts_da.coords['depth'].attrs['units'] = 'meters'
            v_ts_da.coords['depth'].attrs['postive'] = 'up'
        else:
            # we're extracting a time-series (either a sigma level or a z level), so drop the depth dimension
            u_ts_da = xr.DataArray(u_ts.squeeze(), coords={'time': u_ts_i['time'].values,
                                                 'eta_rho': u_ts_i['eta_rho'].values, 
                                                 'xi_rho': u_ts_i['xi_rho'].values,
                                                 'lon_rho': u_ts_i['lon_rho'].values,
                                                 'lat_rho': u_ts_i['lat_rho'].values
                                                 },
                                          dims=['time'])
            
            v_ts_da = xr.DataArray(v_ts.squeeze(), coords={'time': v_ts_i['time'].values,
                                                 'eta_rho': v_ts_i['eta_rho'].values, 
                                                 'xi_rho': v_ts_i['xi_rho'].values,
                                                 'lon_rho': v_ts_i['lon_rho'].values,
                                                 'lat_rho': v_ts_i['lat_rho'].values
                                                 },
                                          dims=['time'])
        # add the attributes for u,v (note these are set in get_var())
        u_ts_da.attrs = u_ts_i.attrs 
        v_ts_da.attrs = v_ts_i.attrs
        
        # create a new dataset with the extracted time-series at this grid cell 
        ds = xr.Dataset({'u': u_ts_da,'v': v_ts_da, 'h': h})
    
    # write a netcdf file if specified
    if write_nc:
        ds.to_netcdf(fname_nc)
        
    return ds
