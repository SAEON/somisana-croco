import numpy as np
from datetime import timedelta
import xarray as xr
import dask
from datetime import timedelta, datetime
from glob import glob
import sys

def hour_rounder(t):
    """
    Rounds to nearest hour by adding a timedelta hour if minute >= 30
    """
    return t.replace(second=0, microsecond=0, minute=0, hour=t.hour) + timedelta(
        hours=t.minute // 30
    )

def u2rho(u):
    """
    regrid the croco u-velocity from it's native u grid to the rho grid
    u can be 2D, 3D or 4D
    """
    Num_dims=len(u.shape)
    if Num_dims==4:
        # T: Time 
        # D: Depth
        # Mp: Grid y
        # L:  Grid x
        [T, D, Mp, L] = u.shape

        # The u grid is length - 1 compared to rho grid
        # Because u values represent mid points between the rho grid
        # Create a new grid for interpolated u values that is equivalent to the rho grid
        u_rho = np.zeros((T, D, Mp, L + 1))
        
        # Interpolate from u grid to rho grid by summing adjacent u points, and divide by 2
        u_rho[:, :, :, 1 : L] = 0.5 * (u[:, :, :, 0 : L - 1] + u[:, :, :, 1 : L])

        # On the edges of the new u grid, you can't interpolate values.
        # So just copy the closest values for the cells on the edge of the grid
        u_rho[:, :, :, 0] = u_rho[:, :, :, 1]
        u_rho[:, :, :, L] = u_rho[:, :, :, L - 1]
        
    elif Num_dims==3:
        # TorD: Temperature or depth
        [TorD, Mp, L] = u.shape # works if first dimension is time or depth
        u_rho = np.zeros((TorD, Mp, L + 1))
        u_rho[:, :, 1 : L] = 0.5 * (u[:, :, 0 : L - 1] + u[ :, :, 1 : L])
        
        u_rho[:, :, 0] = u_rho[:, :, 1]
        u_rho[:, :, L] = u_rho[:, :, L - 1]
        
    else: # Num_dims==2:
        # 2D grid - no time/depth information (possibly a surface level at a single time step)
        [Mp, L] = u.shape
        u_rho = np.zeros((Mp, L + 1))
        u_rho[:, 1 : L] = 0.5 * (u[:, 0 : L - 1] + u[ :, 1 : L])
        u_rho[:, 0] = u_rho[:, 1]
        u_rho[:, L] = u_rho[:, L - 1]
        
    return u_rho

def v2rho(v):
    """
    regrid the croco v-velocity from it's native v grid to the rho grid
    v can be 2D, 3D or 4D

    Refer to u2rho for helpful descriptions and comments
    """
    Num_dims=len(v.shape)
    if Num_dims==4:
        [T, D, M, Lp] = v.shape
        v_rho = np.zeros((T, D, M + 1, Lp))

        v_rho[:, :, 1 : M, :] = 0.5 * (v[:, :, 0 : M - 1, :] + v[:, :, 1 : M, :])
        v_rho[:, :, 0, :] = v_rho[:, :, 1, :]
        v_rho[:, :, M, :] = v_rho[:, :, M - 1, :]
        
    elif Num_dims==3:
        [TorD, M, Lp] = v.shape # works if first dimension is time or depth
        v_rho = np.zeros((TorD, M + 1, Lp))

        v_rho[:, 1 : M, :] = 0.5 * (v[:, 0 : M - 1, :] + v[:, 1 : M, :])
        v_rho[:, 0, :] = v_rho[:, 1, :]
        v_rho[:, M, :] = v_rho[:, M - 1, :]

        
    else: # Num_dims==2:
        [M, Lp] = v.shape
        v_rho = np.zeros((M + 1, Lp))

        v_rho[1 : M, :] = 0.5 * (v[0 : M - 1, :] + v[1 : M, :])
        v_rho[0, :] = v_rho[1, :]
        v_rho[M, :] = v_rho[M - 1, :]
        
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
        # our essential vars are the 'var_str' (obviously) plus some others 
        # we need for things like vertical interp, rotating vectors, computing vorticity etc
        # all the lon* lat* and mar* variables are automatically included in the drop_non_essential_vars_pop sub-function
        static_vars=['s_rho', 's_w', 'sc_r', 'sc_w', 'Cs_r', 'Cs_w', 
                        'hc', 'angle', 'h', 'f', 'pn', 'pm',
                        'Vtransform','theta_s','theta_b',
                        'lon_rho', 'lat_rho', 'mask_rho',
                        'lon_u', 'lat_u', 'lon_v', 'lat_v']
        if var_str in static_vars:
            # no need for open_mfdataset, which can be slow
            # this is here just in case you want to use get_var() and not
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
    s_rho = ds.s_rho.values  # Vertical levels
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

def get_time(fname,ref_date,time_lims=slice(None)):
    ''' 
        fname = CROCO output file (or file pattern to use when opening with open_mfdataset())
        ref_date = reference date for the croco run as a datetime object
        time_lims = optional list of two datetimes i.e. [dt1,dt2], which define the range of times to extract
                    If slice(None), then all time-steps are extracted
    '''
    ds = get_ds(fname)
    
    time = ds.time.values

    time_dt = []
    for t in time:
        date_now = ref_date + timedelta(seconds=np.float64(t))
        time_dt.append(date_now)
    
    if not isinstance(time_lims,slice):
        indx_lims = find_nearest_time_indx(time_dt,time_lims)
        time_lims = slice(indx_lims[0],indx_lims[-1]+1) # +1 to make indices inclusive
        
    time_dt = time_dt[time_lims]
    
    ds.close()
    return time_dt

def tstep_to_slice(fname, tstep, ref_date):
    
    # check if tstep input is instance of datetime, 
    # in which case convert it/them into the correct time index/indices
    if isinstance(np.atleast_1d(tstep)[0],datetime):
        if ref_date is None:
            print('ref_date is not defined - using default of 2000-01-01')
            ref_date=datetime(2000,1,1)
        time_croco = get_time(fname,ref_date)
        tstep = find_nearest_time_indx(time_croco,tstep)
        
    # get the time indices for input to ds.isel()
    if not isinstance(tstep,slice):
        if isinstance(tstep,int):
            # make sure tstep is a list, even if it's a single integer
            # this is a hack to make sure we keep the time dimension 
            # after the ds.isel() step below, even though it's a single index
            # https://stackoverflow.com/questions/52190344/how-do-i-preserve-dimension-values-in-xarray-when-using-isel
            tstep = [tstep]  
        elif len(tstep)==2:
            # convert the start and end limits into a slice
            tstep = slice(tstep[0],tstep[1]+1) # +1 to make indices inclusive 
    
    return tstep

    
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
    '''
    
    # ----------------------------------------------
    # Prepare indices for slicing in ds.isel() below
    # ----------------------------------------------
    #
    # for each of the input dimensions we check the format of the input 
    # and construct the appropriate slice to extract
    #
    
    # get a slice object for time if not already a slice
    tstep = tstep_to_slice(fname, tstep, ref_date)
    
    # as per time, make sure we keep the eta_rho/xi dimensions after the ds.isel() step below
    # this greatly simplifies further functions for depth interpolation 
    # as we know the number of dimensions, even if some of them are single length
    # https://stackoverflow.com/questions/52190344/how-do-i-preserve-dimension-values-in-xarray-when-using-isel
    if not isinstance(eta_rho,slice):
        eta_rho = [eta_rho]
        if var_str=='u' or var_str=='v':
            print('rather use get_ts_uv() for extracting a time-series of u/v data')
            sys.exit()
    if not isinstance(xi_rho,slice):
        xi_rho = [xi_rho]

    # it gets a bit convoluted for the vertical levels 
    # as we have the option of a constant z level which needs interpolation...
    # start by defining variable 'level_for_isel' which is as it sounds
    if not isinstance(level,slice):
        # so level is a single number
        if level >= 0:
            # I'm intentionally not putting 'level' in [] as per the other single indices above
            # this means the ds.isel() step will drop the vertical dimension and 
            # it will be like we're extracting a 2D variable
            level_for_isel = level 
        else:
            # we'll need to do vertical interpolations later so we'll need to initially extract all the sigma levels
            level_for_isel = slice(None)
    else:
        level_for_isel = level # the full slice by definition of the logic
    
    # -----------------
    # Extract the data
    # -----------------
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
    
    # this step can take some time
    # TODO: can we not do this here, and work purely with xarray datasets?
    print('extracting data for '+var_str)
    var = ds[var_str].values
    
    # regrid u, v data onto the rho grid
    if var_str=='u':
        var=u2rho(var)
    if var_str=='v':
        var=v2rho(var)
    
    # ---------------------------
    # Do vertical interpolations
    # ---------------------------
    #
    if len(var.shape)==4 and not isinstance(level,slice):
        print('doing vertical interpolations')
        # given the above checks in the code, here we should be dealing with a 3D variable 
        # and we want a hz slice at a constant depth level
        z=get_depths(ds)
        # z is on the rho grid, but u and v are already regridded to the rho grid above so no need to catch anything here
        T,D,M,L=var.shape
        var_out=np.zeros((T,M,L))
        for t in np.arange(T):
            var_out[t,:,:]=hlev(var[t,::], z[t,::], level)
        var=var_out
    
    # --------
    # Masking
    # --------
    if isinstance(eta_rho,slice) and isinstance(xi_rho,slice):
        _,_,mask=get_lonlatmask(fname,type='r', # u and v are already regridded to the rho grid so can spcify type='r' here
                                eta_rho=eta_rho,
                                xi_rho=xi_rho) # var_str will define the mask type (u,v, or rho)
    else:
        mask=1
    # it looks like numpy is clever enough to use the 2D mask on a 3D or 4D variable
    # that's useful!
    var=var.squeeze()*mask
    
    return var

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
    # u=u2rho(u)
    # v=v2rho(v)
    
    # grid angle
    angle=get_var(fname, 'angle',
                  eta_rho=eta_rho,
                  xi_rho=xi_rho,
                  ) 
    
    # Use the grid angle to rotate the vectors
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)

    # Refer to https://en.wikipedia.org/wiki/Rotation_matrix
    # although 'angle' is 2D, numpy is clever enough for this to work even if u_rho and v_rho are 3D or 4D
    u_out = u*cos_a - v*sin_a
    v_out = v*cos_a + u*sin_a
    
    # Return east / north vector vector components instead of x / y components
    return u_out,v_out

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

def get_ts(fname, var, lon, lat, ref_date, depth=-1, i_shifted=0, j_shifted=0, time_lims=slice(None)):
    """
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
            - i_shifted         :number of grid cells to shift along the xi axis, useful if input lon,lat is on land mask or if input depth is deeper than model depth 
            - j_shifted         :number of grid cells to shift along the eta axis, (similar utility to i_shifted)
            - time_lims         :time step indices to extract 
                                 two values in a list e.g. [dt1,dt2], in which case the range between the two is extracted
                                 If slice(None), then all time-steps are extracted

            Returns:
            - time_model, data_model,lat_mod,lon_mod,h
    """
    
    if var=='u' or var=='v':
        print('rather use get_ts_uv() for extracting a time-series of u/v data')
        sys.exit()
    
    # get the model time
    time_model = get_time(fname, ref_date, time_lims=time_lims)
    
    # finds the rho grid indices nearest to the input lon, lat
    j, i = find_nearest_point(fname, lon, lat) 
    
    # apply the shifts along the xi and eta axis
    i = i+i_shifted
    j = j+j_shifted
    
    # get the data from the model
    # But first check the model depth against the input depth
    h = get_var(fname,"h",eta_rho=j,xi_rho=i)
    if h<-depth:
        print('The model depth is shallower than input/(in situ) depth!!')
        print('Were rather extracting the bottom sigma layer of the model.')
        data_model = get_var(fname, var,
                             tstep=time_lims,
                             level=0, # extracts bottom layer
                             eta_rho=j,
                             xi_rho=i,
                             ref_date=ref_date)
    else:
        data_model = get_var(fname, var,
                             tstep=time_lims,
                             level=depth,
                             eta_rho=j,
                             xi_rho=i,
                             ref_date=ref_date)
        
    # get the model lon, lat data for the time-series we just extracted
    # (useful for comparing against the input lon, lat values)
    lat_mod =  get_var(fname,"lat_rho",eta_rho=j,xi_rho=i)
    lon_mod =  get_var(fname,"lon_rho",eta_rho=j,xi_rho=i)
        
    return time_model, data_model,lat_mod,lon_mod,h

def get_ts_uv(fname, lon, lat, ref_date, depth=-1, i_shifted=0, j_shifted=0, time_lims=slice(None)):
    """
           Extract a timeseries of u,v from the model:
                   
            Parameters:
            - fname             :CROCO output file name (or file pattern to be used with open_mfdataset())
            - lat               :latitude of time-series
            - lon               :longitude of time-series
            - ref_date          :reference datetime used in croco runs
            - depth             :vertical level to extract
                                 If >= 0 then a sigma level is extracted 
                                 If <0 then a z level in meters is extracted
            - i_shifted         :number of grid cells to shift along the xi axis, useful if input lon,lat is on land mask or if input depth is deeper than model depth 
            - j_shifted         :number of grid cells to shift along the eta axis, (similar utility to i_shifted)
            - time_lims         :time step indices to extract 
                                 two values in a list e.g. [dt1,dt2], in which case the range between the two is extracted
                                 If slice(None), then all time-steps are extracted
                                 
            Returns:
            - time_model, u, v, lat_mod, lon_mod, h
              (u,v are rotated from grid-aligned to east-north components)
    """
    # get the model time
    time_model = get_time(fname, ref_date, time_lims=time_lims)
    
    # finds the rho grid indices nearest to the input lon, lat
    j, i = find_nearest_point(fname, lon, lat) 
    
    # apply the shifts along the xi and eta axis
    i = i+i_shifted
    j = j+j_shifted
    
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
    
    # get the model lon, lat data for the time-series we just extracted
    # (useful for comparing against the input lon, lat values)
    lat_mod =  get_var(fname,"lat_rho",eta_rho=j,xi_rho=i)
    lon_mod =  get_var(fname,"lon_rho",eta_rho=j,xi_rho=i)
        
    return time_model, u, v, lat_mod, lon_mod, h

# %%

def get_profile(fname, var, lon, lat, ref_date, time_lims=slice(None),depths=None):
    """
           Extract a timeseries from the model:
                   
            Parameters:
            - fname             :CROCO output file name (or file pattern to be used with open_mfdataset())
            - var               :variable name (string) in the CROCO output file(s)
                                 (not intended for use with u,v variables - rather use get_ts_uv())
            - lat               :latitude of time-series
            - lon               :longitude of time-series
            - ref_date          :reference datetime used in croco runs
            - i_shifted         :number of grid cells to shift along the xi axis, useful if input lon,lat is on land mask or if input depth is deeper than model depth 
            - j_shifted         :number of grid cells to shift along the eta axis, (similar utility to i_shifted)
            - time_lims         :time step indices to extract 
                                 it can be a single integer (starting at zero) or datetime
                                 or two values in a list e.g. [dt1,dt2], in which case the range between the two is extracted
                                 If slice(None), then all time-steps are extracted

            Returns:
            - time_model, profile ,lat_mod ,lon_mod, h
    """
    time_model = get_time(fname, ref_date, time_lims=time_lims)
    #find_nearest_point finds the nearest point in the model to the model grid lon, lat extracted from the model grid input.
    j, i = find_nearest_point(fname, lon, lat) 
    
    # get a slice object for time if not already a slice
    time_lims = tstep_to_slice(fname, time_lims, ref_date)

    ds = get_ds(fname, var)
    ds = ds.isel(time=time_lims,
                        eta_rho=slice(j,j+1), # making it a slice to maintain the spacial dimensions for input to get_depths()
                        xi_rho=slice(i,i+1)
                        )

    if depths is None:
        depths_out = np.squeeze(get_depths(ds))
        profile = get_var(fname, var,
                              tstep=time_lims,
                              eta_rho=j,
                              xi_rho=i,
                              ref_date=ref_date)
    else:
        depths_out = depths
        profiles = np.zeros((len(time_model),len(depths)))
        for index, depth in enumerate(depths):
            profile = get_var(fname, var,
                              tstep=time_lims,
                              level=depth,
                              eta_rho=j,
                              xi_rho=i,
                              ref_date=ref_date)
            profiles[:,index]=profile
        
    print(profiles)
    lat_mod =  get_var(fname,"lat_rho",eta_rho=j,xi_rho=i)
    lon_mod =  get_var(fname,"lon_rho",eta_rho=j,xi_rho=i)

    return time_model, depths_out, profile, lat_mod, lon_mod
