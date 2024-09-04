import numpy as np
import xarray as xr

from scipy.sparse import lil_matrix
from pyamg import ruge_stuben_solver, solve
from pyamg.krylov import bicgstab
from pyamg.util.linalg import norm

import xrft

import gridop as gop

###################################################

def density(temp, salt, z, rho0=1025, split=True):
    """Calculate the density [kg/m^3] as calculated in CROCO.

    Args
    - temp: (DataArray) tempemperature [Celsius]
    - salt: (DataArray) Salinity
    - z: (DataArray) Depth [m]. 
    - rho0: (float, optional, default=1025) reference density
    - split: (Boolean, optional, default=True) Activate EOS splitting of seawater compressibility effect
    
    Return
    -------
    - DataArray  of calculated density on rho/rho grids. Output is `[temp,Z,Y,X]`.

    Notes
    -----
    Equation of state based on CROCO Nonlinear/rho_eos.F

    Examples
    --------
    >>> density(temp, salt, z)
    """
    
    r00=999.842594
    r01=6.793952E-2
    r02=-9.095290E-3
    r03=1.001685E-4
    r04=-1.120083E-6
    r05=6.536332E-9
    r10=0.824493
    r11=-4.08990E-3
    r12=7.64380E-5
    r13=-8.24670E-7
    r14=5.38750E-9
    rS0=-5.72466E-3
    rS1=1.02270E-4
    rS2=-1.65460E-6
    r20=4.8314E-4
    K00=19092.56
    K01=209.8925
    K02=-3.041638
    K03=-1.852732e-3
    K04=-1.361629e-5
    K10=104.4077
    K11=-6.500517
    K12=0.1553190
    K13=2.326469e-4
    KS0=-5.587545
    KS1=+0.7390729
    KS2=-1.909078e-2
    B00=0.4721788
    B01=0.01028859
    B02=-2.512549e-4
    B03=-5.939910e-7
    B10=-0.01571896
    B11=-2.598241e-4
    B12=7.267926e-6
    BS1=2.042967e-3
    E00=+1.045941e-5
    E01=-5.782165e-10
    E02=+1.296821e-7
    E10=-2.595994e-7
    E11=-1.248266e-9
    E12=-3.508914e-9
    
    g = 9.81
    Tref=3.8
    Sref=34.5
    qp2=1.72e-05
    
    def _K0(temp,salt):
        return temp*( K01+temp*( K02+temp*( K03+temp*K04 ))) + salt*( K10+temp*( K11+temp*( K12+temp*K13 )) + np.sqrt(salt)*( KS0+temp*( KS1+temp*KS2 )))

    def _K1(temp,salt):
        return  B00+temp*(B01+temp*(B02+temp*B03)) +salt*( B10+temp*( B11 + temp*B12 )+np.sqrt(salt)*BS1 )

    def _K2( temp,salt):
        return  E00+temp*(E01+temp*E02) +salt*(E10+temp*(E11+temp*E12))

    def _K( temp,salt,z):
        return K00 + _K0(temp,salt) + _K1(temp,salt) * z + _K2(temp,salt) * z**2

    def _qp1(T, S, Tref, Sref):
        return 0.1*(rho0+_rho1(T,S))* (1/(K00 + _K0(T,S)) - 1/(K00 + _K0(Tref,Sref)))

    def _rho1(temp,salt):
        return  r00-rho0 + temp*( r01+temp*( r02+ temp*( r03+ temp*( r04+ temp*r05 ))))+salt*( r10+ temp*( r11+ temp*( r12+ temp*(r13+temp*r14 )))+np.sqrt(salt)*(rS0+ temp*(rS1+ temp*rS2 ))+ salt*r20 )
    
    if split:
        rho = _rho1(temp,salt) + _qp1(temp, salt, Tref, Sref) * z *(1. - qp2 * z)   
    else:
        rho = (_rho1(temp,salt) ) / (1 - 0.1 * z / _K(temp,salt,z))
        
    return rho.rename("rho")
    
###################################################
    
# relative vorticity
def relative_vorticity_z(model, ds=None, xgrid=None, u=None, v=None, f=None):
    
    """
    Compute relative vorticity normalized by f
    input: 
        model : instance of the Model class defined in the model.py module
        ds    : xarray DataSet
        xgrid : xgcm grid
        u     : xarray DataArray x-current component
        v     : xarray DataArray y-current component
    """
    
    xgrid = model.xgrid if xgrid is None else xgrid
    ds = model.ds if ds is None else ds    
    
    u = ds.u if u is None else u
    hgrid, vgrid = gop.get_grid_point(u)
    if hgrid != 'u': u = gop.to_u(u, xgrid)
        
    v = ds.v if v is None else v
    hgrid, vgrid = gop.get_grid_point(v)
    if hgrid != 'v': v = gop.to_v(v, xgrid)
        
    f = ds.f if f is None else f
    hgrid, vgrid = gop.get_grid_point(f)
    if hgrid != 'f': f = gop.to_psi(f, xgrid)
    
    xi = ((-xgrid.derivative(u, 'y') 
           + xgrid.derivative(v, 'x')
          )/f)
    return xi.rename('vorticity')


###################################################

def relative_vorticity_sigma(
    u,
    v,
    xgrid,
    z=None,
    hboundary="extend",
    hfill_value=None,
    sboundary="extend",
    sfill_value=None,
):
    """Calculate the vertical component of the relative vorticity [1/s]

    Parameters
    ----------
    u: DataArray
        xi component of velocity [m/s]
    v: DataArray
        eta component of velocity [m/s]
    xgrid: xgcm.grid
        Grid object associated with u, v
    z: DataArray
        z coordinate 
    hboundary: string, optional
        Passed to `grid` method calls; horizontal boundary selection
        for calculating horizontal derivatives of u and v.
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
        for calculating horizontal derivatives of u and v.
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

    Returns
    -------
    DataArray of vertical component of relative vorticity psi/w grids.
    Output is `[T,Z,Y,X]`.

    Notes
    -----
    relative_vorticity = v_x - u_y

    Examples
    --------
    >>> xroms.relative_vorticity(u, v, xgrid)
    """

    assert isinstance(u, xr.DataArray), "u must be DataArray"
    assert isinstance(v, xr.DataArray), "v must be DataArray"

    # assign z coordinate to v
    zcoord = [c for c in v.coords if c.startswith("z")]
    if not zcoord and z is not None: 
        z_v = gop.to_grid_point(z, xgrid, hcoord='v', vcoord='r')
        v = v.assign_coords({"z_v": z_v})
    dvdx = gop.hgrad(
        v,
        xgrid,
        which="x",
        hboundary=hboundary,
        hfill_value=hfill_value,
        sboundary=sboundary,
        sfill_value=sfill_value,
    )
    
    # assign z coordinate to u
    zcoord = [c for c in u.coords if c.startswith("z")]
    if not zcoord and z is not None: 
        z_u = gop.to_grid_point(z, xgrid, hcoord='u', vcoord='r')
        u = u.assign_coords({"z_u": z_u})
    dudy = gop.hgrad(
        u,
        xgrid,
        which="y",
        hboundary=hboundary,
        hfill_value=hfill_value,
        sboundary=sboundary,
        sfill_value=sfill_value,
    )

    var = dvdx - dudy

    var.attrs["name"] = "vorticity"
    var.attrs["long_name"] = "vertical component of vorticity"
    var.attrs["units"] = "1/s"
    var.name = var.attrs["name"]

    return var.squeeze()
    
###################################################

def ertel_pv(xgrid, u, v, w, rho, z, f, rho0=None, typ='ijk'):
    """
    
       epv    - The ertel potential vorticity with respect to property 'lambda'
    
                                           [ curl(u) + f ]
       -  epv is given by:           EPV = --------------- . del(lambda)
                                                rho
    
       -  pvi,pvj,pvk - the x, y, and z components of the potential vorticity.
    
      -  Ertel PV is calculated on horizontal rho-points, vertical w-points.
    
    Args:
        xgrid: (xgcm.grid) Grid object associated with u, v
        u: (DataArray) xi component of velocity [m/s]
        v: (DataArray) eta component of velocity [m/s]
        w: (DataArray) sigma component of velocity [m/s]
        rho: (DataArray) density
        z: (DataArray) Depth at rho points [m]. 
        f: (DataArray) Coriolis parameter 
        rho0: (float) Reference density
        typ : (string) which components of the potential vorticity to compute

    Returns:
        DataArray: The ertel potential vorticity
    
    Adapted from rob hetland.
    
    """

    assert isinstance(u, xr.DataArray), "u must be DataArray"
    assert isinstance(v, xr.DataArray), "v must be DataArray"
    assert isinstance(w, xr.DataArray), "w must be DataArray"
    assert isinstance(rho, xr.DataArray), "rho must be DataArray"
    assert isinstance(z, xr.DataArray), "z must be DataArray"
    assert isinstance(f, xr.DataArray), "f must be DataArray"
    
    dz = xgrid.diff(z,'z')
    rho0 = 1025 if rho0 is None else rho0

    if 'k' in typ:

        # Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*drho/dz
        # Compute d(v)/d(xi) at PSI-points.
        dvdxi = xgrid.derivative(v,'x')

        # Compute d(u)/d(eta) at PSI-points.
        dudeta = xgrid.derivative(u,'y')

        # Compute d(rho)/d(z) at horizontal RHO-points and vertical W-points
        drhodz = xgrid.diff(rho,'z') / dz

        #  Compute Ertel potential vorticity <k hat> at horizontal RHO-points and
        #  vertical W-points. 
        omega = dvdxi - dudeta
        omega = f + gop.to_rho(omega, xgrid)
        pvk = xgrid.interp(omega,'z') * drhodz
        del dvdxi, dudeta, drhodz, omega
    else:
        pvk = 0.

    if 'i' in typ:

        #  Ertel potential vorticity, term 2: (dw/dy - dv/dz)*(drho/dx)
        #  Compute d(w)/d(y) at horizontal RHO-points and vertical RHO-points
        dwdy = xgrid.derivative(w,'y')
        dwdy = gop.to_rho(dwdy, xgrid)

        #  Compute d(v)/d(z) at horizontal RHO-points and vertical W-points
        dz_v = xgrid.interp(dz,'y')
        dvdz = xgrid.diff(v,'z') / dz_v
        dvdz = gop.to_rho(dvdz, xgrid)

        #  Compute d(rho)/d(xi) at horizontal RHO-points and vertical RHO-points
        drhodx = xgrid.derivative(rho,'x')
        drhodx = gop.to_rho(drhodx, xgrid)

        #  Add in term 2 contribution to Ertel potential vorticity at horizontal RHO-points and
        #  vertical W-points.
        pvi = (gop.to_s_w(dwdy, xgrid) - dvdz) * gop.to_s_w(drhodx, xgrid)
        del dwdy, dz_v, dvdz, drhodx
    else:
        pvi = 0.

    if 'j' in typ:

        #  Ertel potential vorticity, term 3: (du/dz - dw/dx)*(drho/dy)
        #  Compute d(u)/d(z) at horizontal RHO-points and vertical W-points
        dz_u = xgrid.interp(dz, 'x')
        dudz = xgrid.diff(u,'z') / dz_u
        dudz = gop.to_rho(dudz, xgrid)

        #  Compute d(w)/d(x) at horizontal RHO-points and vertical RHO-points
        dwdx = xgrid.derivative(w,'x')
        dwdx = gop.to_rho(dwdx, xgrid)

        #  Compute d(rho)/d(eta) at horizontal RHO-points and vertical RHO-points
        drhodeta = xgrid.derivative(rho,'y')
        drhodeta = gop.to_rho(drhodeta, xgrid)

        #  Add in term 3 contribution to Ertel potential vorticity at horizontal RHO-points and
        #  vertical W-points..
        pvj =  (dudz - gop.to_s_w(dwdx, xgrid)) * gop.to_s_w(drhodeta, xgrid)
        del dz_u, dudz, dwdx, drhodeta

    else:
        pvj = 0.

    #
    #
    # Sum potential vorticity components, and divide by rho0
    #
    pvi = pvi / rho0
    pvj = pvj / rho0
    pvk = pvk / rho0
    pv = pvi + pvj + pvk

    pv = pv.assign_coords(coords={"lon":rho.lon})
    pv = pv.assign_coords(coords={"lat":rho.lat})
    z_w = xgrid.interp(z, 'z') #.fillna(0.)
    pv = pv.assign_coords(coords={"z_w":z_w})

    return pv.squeeze().rename('pv')

###################################################

def dtempdz(xgrid, temp, z):
    """
    Compute dT/dz at horizontal rho point/vertical w point
    ds : dataset, containing T field
    z : xarray.DataArray, z in meters at rho points
    time : int, time index 
    """

    assert isinstance(temp, xr.DataArray), "temp must be DataArray"
    assert isinstance(z, xr.DataArray), "z must be DataArray"

    # compute z coordinates at w level
    z_w = xgrid.interp(z,'z').squeeze()

    dtempdz = (xgrid.diff(temp,'z') / xgrid.diff(z,'z')).squeeze()
    if 'lon' in temp.coords: dtempdz = dtempdz.assign_coords(coords={"lon":temp.lon})
    if 'lat' in temp.coords: dtempdz = dtempdz.assign_coords(coords={"lat":temp.lat})
    dtempdz = dtempdz.assign_coords(coords={"z_w":z_w})
    return dtempdz.rename('dtdz')
    
###################################################

def richardson(xgrid, u, v, rho, z, rho0=None):
    """
    Ri is given by:      N²/((du/dz)² - (dv/dz)²)
         with N = sqrt(-g/rho0 * drho/dz)
    Ri is calculated at RHO-points and w level

    ds : dataset, which contains 3D u,v and rho fields
    z : xarray datarray, z depths at rho levels
    time : int, time index
    """

    assert isinstance(u, xr.DataArray), "u must be DataArray"
    assert isinstance(v, xr.DataArray), "v must be DataArray"
    assert isinstance(rho, xr.DataArray), "rho must be DataArray"
    assert isinstance(z, xr.DataArray), "z must be DataArray"

    g = 9.81
    rho0 = 1027 if 'rho0' is None else rho0
    
    # compute Z at w levels
    z_w = xgrid.interp(z,'z').squeeze()
    
    u = gop.to_grid_point(u, xgrid, hcoord="r", vcoord="r")
    v = gop.to_grid_point(v, xgrid, hcoord="r", vcoord="r")
    z = gop.to_grid_point(z, xgrid, hcoord="r", vcoord="r")

    # N2 = get_N2(model, ds=ds, rho=rho, z=z, g=g)
    N2 = get_N2(xgrid, rho, z, rho0=rho0)
    dudz = xgrid.diff(u,'z') / xgrid.diff(z,'z')
    dvdz = xgrid.diff(v,'z') / xgrid.diff(z,'z')

    Ri = np.log10(N2 / (dudz**2 +  dvdz**2)).squeeze()
    if 'lon' in rho.coords: Ri = Ri.assign_coords(coords={"lon":rho.lon})
    if 'lat' in rho.coords: Ri = Ri.assign_coords(coords={"lat":rho.lat})
    Ri = Ri.assign_coords(coords={"z_w":z_w})
    return Ri.rename('Ri')
    

###################################################

def get_N2(xgrid, rho, z, rho0=None):
    """ Compute square buoyancy frequency N2 
    ... doc to be improved
    """
    
    if xgrid is None: xgrid = model.xgrid
    
    rho0 = 1027 if rho0 is None else rho0
    g = 9.81
    
    N2 = -g/rho0 * xgrid.diff(rho, 'z', boundary='fill', fill_value=np.NaN) \
            / xgrid.diff(z, 'z', boundary='fill', fill_value=np.NaN)
    # cannot find a solution with xgcm, weird
    N2 = N2.fillna(N2.shift(s_w=-1))
    N2 = N2.fillna(N2.shift(s_w=1))
    return N2.rename('N2')
    
###################################################

def poisson_matrix(pm,pn):    
    """
    Initialize the elliptic equation matrix (d_xx + d_yy)
    Input :
        - pm : (ndarray) 1/dx coefficents
        - pn : (ndarray) 1/dy coefficents
    Output:
        return a compressed sparse matrix of the laplacian operator
    """
    # elliptic equation matrix:  d_xx + d_yy
    [nx,ny] = pm.shape
    
    ndim = ny*nx
    i_s  = ny
    js  = 1
    A=lil_matrix((ndim,ndim))
    ############################
    
    for i in range(nx): 
        for j in range(ny): 
                idx = i*ny + j;
                diag = 0.;
                if j>0:
                    dy2i = 0.5*(pn[i,j]+pn[i,j-1])*pn[i,j]
                    A[idx,idx-js] = dy2i;
                    diag -= dy2i;
                if i>0:
                    dx2i = 0.5*(pm[i,j]+pm[i-1,j])*pm[i,j]
                    A[idx,idx-i_s] = dx2i;
                    diag -= dx2i;
                if i<nx-1:
                    dx2i = 0.5*(pm[i,j]+pm[i+1,j])*pm[i,j]
                    A[idx,idx+i_s] = dx2i;
                    diag -= dx2i;
                if j<ny-1:
                    dy2i = 0.5*(pn[i,j]+pn[i,j+1])*pn[i,j]
                    A[idx,idx+js] = dy2i;
                    diag -= dy2i;
                A[idx,idx] = diag
        
    return A.tocsr()

###################################################

def get_streamfunction(model, pm, pn, pv,
                       tol=1e-5, solver='classical', verb=False):
    """
    !!!DOES NOT WORK YET
    
    Compute the stream function from the relative vorticity
    Invert the laplacian to solve the poisson equation Ax=b
    A is the horizontal laplacian, b is the vorticity
    Input:
        - pm : (DataArray) 1/dx metric at f grid point
        - pn : (DataArray) 1/dy metric at f grid point
        - pv : (DataArray) relative vorticity at f grid point
        - tol : (Float) Tolerance for stopping criteria
        - solver : (String) type of solver
            'classical' : multilevel solver using Classical AMG (Ruge-Stuben AMG)
            ' bicgstab' : Biconjugate Gradient Algorithm with Stabilization
        - verb : (Boolean) verbose mode
    Output:
        (DataArray) the computed streamfunction at f grid point
    """

    if np.any(np.isnan(pv)): 
        print("Can't inverse the laplacian, non compact domain, pv contains nan values")
        return None

    # pm pn at the vorticity f grid point
    pm = pm.fillna(0)
    pn = pn.fillna(0)
    
    #######################################################
    #Create matrix A
    #######################################################
    
    if verb: print('creating matrix A')
    A = poisson_matrix(pm.values,pn.values)

    #######################################################
    #Solve matrix A
    #######################################################

    if verb: print('creating matrix b')
    b = -1. * pv.values.flatten() # right hand side
    
    if solver=='classical':
        # multilevel solver using Classical AMG (Ruge-Stuben AMG)
        ml = ruge_stuben_solver(A)                # construct the multigrid hierarchy
        if verb: print("ml=",ml)                       # print hierarchy information
        x = solve(A,b,verb=verb,tol=tol)         # solve Ax=b  
    else:
        # Biconjugate Gradient Algorithm with Stabilization
        (x,flag) = bicgstab(A,b, tol=tol)
        
    if verb: print(f'residual: {norm(b - A*x):.6}') # compute norm of residual vector
    
    dims = gop.get_spatial_dims(pv)
    coords = gop.get_spatial_coords(pv)
    chi = xr.DataArray(
        data=x.reshape(pv.shape),
        dims=[dims['y'], dims['x']],
        coords={coords['lon']:pv[coords['lon']], coords['lat']:pv[coords['lat']]}
        )
    # mask_f = gop.to_psi(model.ds.mask,model.xgrid)
    # chi = chi * mask_f
    dx = 1./pm.where(pm>0,np.nan)
    dy = 1./pn.where(pm>0,np.nan)
    chi = chi / (dx*dy).sum().values
    return chi.rename('streamfct')

###################################################

def get_p(xgrid, rho, z_w, z_r, rho0=None):
    """ 
    Compute (not reduced) pressure by integration from the surface, 
    taking rho at rho points and giving results on rho points (z grid).

    Parameters
    ----------
    grid : xgcm grid
    rho  : Density (DataArray)
    z_w  : depth at vertical w points (DataArray)
    z_r  : depth at vertical rho points (DataArray)
    g    : acceleration of gravity (float)
    rho0 : mean density (float)

    """

    assert isinstance(rho, xr.DataArray), "rho must be DataArray"
    assert isinstance(z_w, xr.DataArray), "z_w must be DataArray"
    assert isinstance(z_r, xr.DataArray), "z_r must be DataArray"
    
    # useful parameters
    eps = 1.0e-10
        
    g=9.81
    rho0 = 1027 if rho0 is None else rho0
                
    GRho=g/rho0
    HalfGRho=0.5*GRho
    N = rho.sizes['s']
    OneFifth = 1.0/5.0
    OneTwelfth = 1.0/12.0

    # dz and drho on w levels
    dR = xgrid.diff(rho,'z', boundary='extend').rename('dR')
    dZ = xgrid.diff(z_r,'z', boundary='extend').rename('dZ')

    # modified dz and dr on w levels
    dZi = 2. * ( dZ * dZ.shift(s_w=1,fill_value=0) 
            / (dZ + dZ.shift(s_w=1,fill_value=0)) )
    dZi = xr.concat([dZ.isel(s_w=0), dZi.isel(s_w=slice(1,None))], dim='s_w')
    dZi = dZi.chunk(chunks={'s_w':-1})
    dRi = 2. * ( dR * dR.shift(s_w=1,fill_value=0) 
            / (dR + dR.shift(s_w=1,fill_value=0)) )
    dRi = xr.concat([dR.isel(s_w=0), dRi.isel(s_w=slice(1,None))], dim='s_w')
    dRi = dRi.chunk(chunks={'s_w':-1})

    # Pressure at the surface on rho level
    Pn = (g*z_w.isel(s_w=-1) + GRho*( rho.isel(s=-1) 
                                   + 0.5*(rho.isel(s=-1)-rho.isel(s=-2)) 
                                   * (z_w.isel(s_w=-1)-z_r.isel(s=-1)) 
                                   / (z_r.isel(s=-1)-z_r.isel(s=-2))
                                  ) * (z_w.isel(s_w=-1)-z_r.isel(s=-1))
         )

    # Pressure of each slice
    rhoz = (rho + rho.shift(s=-1))*(z_r.shift(s=-1) - z_r)
    dRz = ((xgrid.diff(dRi,'z'))*(z_r.shift(s=-1) - z_r
                                 -OneTwelfth*(xgrid.interp(dZi, 'z'))) )
    dZr = ((xgrid.diff(dZi,'z'))*(rho.shift(s=-1) - rho
                                 -OneTwelfth*(xgrid.interp(dRi,'z'))) )
    Pk = (HalfGRho*( rhoz - OneFifth*(dRz - dZr)))

    # replace pressure at the surface
    Pk = xr.concat([Pk.isel(s=slice(0,-1)), Pn], dim='s')

    # integrate from the surface to bottom
    P = (Pk
            .sortby(Pk.s, ascending=False)                
            .cumsum("s")
            .sortby(Pk.s, ascending=True)
            .assign_coords(z_r=z_r)
        )

    return P.rename('P')

###################################################

def power_spectrum(da, dims, true_phase=True, true_amplitude=True, 
                   window=None, detrend=None, spacing_tol=0.001):
    """
    Compute the spectrum of the dataarray over the dimensions dims. By default, the signal is supposed 
    to be periodic through all these dimensions.
    See the documentation https://xrft.readthedocs.io/en/latest/api.html for more details
    
    Input arguments:
        da : (DataArray) input data 
        dims : (str or list of str) dimensions of da on which to take the FFT
        true_phase : (Boolean) If set to False, standard fft algorithm is applied on signal without 
                     consideration of coordinates
        true_amplitude : (Boolean) If set to True, output is multiplied by the spacing of the transformed 
                        variables to match theoretical FT amplitude
        window : (str) Whether to apply a window to the data before the Fourier transform is taken. 
                 A window will be applied to all the dimensions in dim. Please follow scipy.signal.windows’ 
                 naming convention (ex: 'hann', 'gaussian')
        detrend : (str: 'constant', 'linear' or None) 
                  If constant, the mean across the transform dimensions will be subtracted before 
                  calculating the Fourier transform (FT). If linear, the linear least-square fit will be 
                  subtracted before the FT
    
    Returns the power spectrum of da as a DataArray. The dimensions of the FFT are modified in frequency ('y' -> 'freq_y')
    """
    
    dimcoord={'x_u':'lon_u', 'x':'lon', 'y_v':'lat_v', 'y':'lat', 's':'z', 's_w':'z_w', 't':'t'}

    # verify that dims is a list
    if not isinstance(dims,list): dims = [dims]
    
    # verify that every dimension in dims belong to da
    if not all([d in da.dims for d  in dims]):
        print('There are missing dimensions in the input DataArray')
        return None

    # liste des coordonnées de da
    dacoords = [c for c in da.coords]
    
    # retrieve coordinates of da corresponding to dims
    coords={d:da[v].reset_coords(drop=True) for d in dims for v in dacoords if v.startswith(dimcoord[d])}
      
    # makes the coordinates of da become 1D
    # mean of the coordinates over all the dimensions but one (t(t), lon(x), lat(y)...)
    coords1D={}
    for k,v in coords.items():
        if len(v.shape)>1:
            maindim = [d for d in v.dims if d.startswith(k)]
            removedims = list(set(v.dims) - set(maindim))
            v = v.mean(dim=removedims)         
        coords1D[k] = v
        
    # remove all the coordinates from da
    newda = da.reset_coords(drop=True)
    
    # assign 1D coordinates to newda
    for d in dims:
        newda = newda.assign_coords({d:coords1D[d]})
        
    # Compute FFT
    Fda = xrft.xrft.fft(newda, dim=dims,  
                        detrend=detrend, window=window, 
                        true_phase=true_phase, true_amplitude=true_amplitude,
                        spacing_tol=spacing_tol)
    
    # return power spectra
    return abs(Fda*Fda.conj())

###################################################
