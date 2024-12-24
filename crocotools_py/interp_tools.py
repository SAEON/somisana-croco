import numpy as np
from scipy.spatial import Delaunay
from scipy import interpolate
#
# Interpolation functions for get section in postprocessing.py
#
def csf(sc, theta_s, theta_b):
    """
    Allows use of theta_b > 0 (July 2009)
    is required in z_levels.py
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

    this is adapted (by J.Veitch - Feb 2022) from z_levels.m in roms_tools (by P. Penven)
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


def get_tri_coef(X, Y, newX, newY, verbose=0):
    """
    INPUTS:
    X, Y: Old x and y grid points
    newX, newY: New x and y grid points 
    
    OUTPUTS:
    elem: pointers to 2d gridded data (at lonp,latp locations) from
            which the interpolation is computed (3 for each child point)
    coef: linear interpolation coefficients
    
    Use:
    To subsequently interpolate data from Fp to Fc, the following will work:      
        Fc  = sum(coef.*Fp(elem),3)  
        This line  should come in place of all griddata calls. It avoids repeated triangulations and tsearches (that are done
        with every call to griddata) it should be much faster.
    """


    Xp = np.array([X.ravel(), Y.ravel()]).T
    Xc = np.array([newX.ravel(), newY.ravel()]).T


    #Compute Delaunay triangulation
    if verbose==1: tstart = tm.time()
    tri = Delaunay(Xp)
    if verbose==1: print('Delaunay Triangulation', tm.time()-tstart)

    #Compute enclosing simplex and barycentric coordinate (similar to tsearchn in MATLAB)
    npts = Xc.shape[0]
    p = np.zeros((npts,3))

    points = tri.points[tri.simplices[tri.find_simplex(Xc)]]
    if verbose==1: tstart = tm.time()
    for i in range(npts):

        if verbose==1: print(np.float(i)/npts)

        if tri.find_simplex(Xc[i])==-1:  #Point outside triangulation
             p[i,:] = p[i,:] * np.nan

        else:

            if verbose==1: tstart = tm.time()
            A = np.append(np.ones((3,1)),points[i] ,axis=1)
            if verbose==1: print('append A', tm.time()-tstart)

            if verbose==1: tstart = tm.time()
            B = np.append(1., Xc[i])
            if verbose==1: print('append B', tm.time()-tstart)

            if verbose==1: tstart = tm.time()
            p[i,:] = np.linalg.lstsq(A.T,B.T)[0]
            if verbose==1: print('solve', tm.time()-tstart)

    if verbose==1: print('Coef. computation 1', tm.time()-tstart)

    if verbose==1: tstart = tm.time()
    elem = np.reshape(tri.simplices[tri.find_simplex(Xc)],(newX.shape[0],newY.shape[1],3))
    coef = np.reshape(p,(newX.shape[0],newY.shape[1],3))
    if verbose==1: print('Coef. computation 2', tm.time()-tstart)

    return(elem,coef)


def horiz_interp_delaunay(lonold,latold,varold,lonnew,latnew,elem=0,coef=0):
    """
    Horizontal Interpolation from croco varoables grid to new grid
    
    INPUT:
    lonold: 2D original longitude matrix
    latold: 2D original longitude matrix
    varold: 2D original data matrix
    lonnew: 2D new longitude matrix
    latnew: 2D new longitude matrix
    elem,coef: interpolation coefficients. 

    OUTPUT:
    varnew: variable interpolated onto new grid.
    
    NOTE:
    If elem,coef is == 0, then it computes the coeficients and returns the arrays. 
    If elem,coef is != 0, and it is an input, then we do not compute it so its faster, but we do not return the arrays. 
    """
    if np.all(elem==0):
        [elem,coef] = get_tri_coef(lonold, latold,lonnew, latnew)
        coefnorm=np.sum(coef,axis=2)
        coef=coef/coefnorm[:,:,np.newaxis]
        varnew = np.sum(coef*varold.ravel()[elem],2)
        return elem,coef,varnew

    else:
        varnew = np.sum(coef*varold.ravel()[elem],2)
        return varnew

def croco_3d_interp_delaunay_spline_atrho(lon,lat,h,zeta,mask,grid,\
                                          theta_s,theta_b,hc,N,\
                                          lonnew,latnew,depthnew,
                                          *var):
    """
    Vertical interpolation of CROCO variable from sigma grid to z grid.
    lon: original 2D longitude matrix.
    lat: original 2D latitude matrix.
    h,zeta: CROCO topo (h) and free surface (zeta) 2D matrices.
    grid: string of the grid points at which the variable is poisioned ('rho', 'u', 'v' or 'psi' grid points).
    theta_s,theta_b,hc,N: parameters of croco vertical grid.
    lonnew: new 2D longitude matrix.
    latnew: new 2D latitude matrix.
    depthnew: new 1D vertical grid.
    var: 3D original croco data matrix (can be several variables as long as they are on the same grid points).
    """
    bad_value=999

    #Build croco z levels at rho points
    zvarold =  z_levels(h,zeta,theta_s,theta_b,hc,N,grid,2) # 2 should probilly not be hardcoded into the code (vtransform: 1 or 2)

    #lonnew and latnew 2D horizontal grids
    nz_new=np.size(depthnew)
    if np.ndim(lonnew)==1:
        ny_new=1
        nx_new=np.shape(lonnew)[0]
        lonnew=lonnew[np.newaxis,:]
        latnew=latnew[np.newaxis,:]
    elif np.ndim(lonnew==2):
        ny_new=np.shape(lonnew)[0]
        nx_new=np.shape(lonnew)[1]

    # Horizontal interpolation of the topo
    elem,coef,toponew=horiz_interp_delaunay(lon,lat,h,lonnew,latnew)

    # Horizontal interpolation of the mask nut it is double check with zlev>0
    # It masks points that are really close to the coast
    masknew=horiz_interp_delaunay(lon,lat,mask,lonnew,latnew,elem,coef)
    masknew[np.round(masknew,5)<1]=0
    masknew[np.round(masknew,5)>=1]=1

    # Horizontal interpolation of zlevels and variables
    # add a point at the surface to fall in the range of interpolation    
    varnew_oldz = np.zeros((np.shape(var)[0],N+1,ny_new, nx_new))
    
    #depth of sigma levels at lonnew,latnew
    zvarnew_oldz= np.zeros((N+1,ny_new, nx_new))

    #Vertical interpolation on new vertical grid
    varnew_newz=np.zeros((np.shape(var)[0],np.size(depthnew),ny_new,nx_new))+bad_value

    #multiplying by masknew results in setting to zero values very close to the coasts
    for nvar in range(np.shape(var)[0]):
        tmp=var[nvar]
        for k in range(N):
            varnew_oldz[nvar,k,:,:]=horiz_interp_delaunay(lon,lat,tmp[k],
                                                          lonnew,latnew,elem,coef)*masknew
            zvarnew_oldz[k,:,:] =horiz_interp_delaunay(lon,lat,zvarold[k,:,:],
                                                            lonnew,latnew,elem,coef)*masknew

        varnew_oldz[nvar,N] = varnew_oldz[nvar,N-1]

        for i in range(ny_new):
            for j in range(nx_new):
                nznew_above_zero = np.sum(depthnew>zvarnew_oldz[0,i,j])
                k_start = nz_new - nznew_above_zero
                if nznew_above_zero!=0:
                    f_section = interpolate.interp1d(zvarnew_oldz[:,i,j],varnew_oldz[nvar,:,i,j], kind='cubic')
                    varnew_newz[nvar,k_start:,i,j] = f_section(depthnew[k_start:])


        varnew_newz = np.ma.masked_where((varnew_newz==bad_value),varnew_newz)

    return toponew,masknew,varnew_newz
