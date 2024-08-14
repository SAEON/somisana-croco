from progressbar import progressbar
import numpy as np

import warnings
warnings.filterwarnings("ignore")
def scoord2z(point_type, zeta, topo,theta_s, theta_b,N,hc,scoord='new2008',Dcrit=0.2):
    '''
    scoord2z finds z at either rho or w points (positive up, zero at rest surface)

    Inputs:
      point_type        'r' or 'w'
      zeta               sea surface height
      topo              array of depths (e.g., from grd file)
      theta_s           surface focusing parameter
      theta_b           bottom focusing parameter
      N                 number of vertical rho-points
      hc                critical depth
      scoord            'new2008' :new scoord 2008  or 'old1994' for Song scoord
    
    Outputs:
      z                 sigma coordinates
      Cs                Cs parameter
    '''
    def CSF(sc,theta_s,theta_b):
        '''
        Allows use of theta_b > 0 (July 2009)
        '''
        one64 = np.float64(1)
        if theta_s > 0.:
            csrf = ((one64-np.cosh(theta_s*sc))
                       /(np.cosh(theta_s)-one64))
        else:
            csrf = -sc**2
        sc1 = csrf+one64
        if theta_b > 0.:
            Cs = ((np.exp(theta_b*sc1)-one64)
                /(np.exp(theta_b)-one64)-one64)
        else:
            Cs = csrf
        return Cs

    N = np.float64(N)
    if isinstance(zeta,float):
        zeta = np.ones(topo.shape)*zeta
    if scoord not in 'new2008':
        cff1 = 1./np.sinh(theta_s)
        cff2 = 0.5/np.tanh(0.5*theta_s)
    sc_w = (np.arange(N+1,dtype=np.float64)-N)/N
    sc_r = ((np.arange(1,N+1,dtype=np.float64))-N-0.5)/N

    if 'w' in point_type:
        sc = sc_w
        N += 1. # add a level
    else:
        sc = sc_r
        
    if len(np.array(zeta).shape)>2: # case zeta is 3-D (in time)
        z  = np.empty((int(zeta.shape[0]),) + (int(N),) + topo.shape, dtype=np.float64)
    else:
        z  = np.empty((int(N),) + topo.shape, dtype=np.float64)

    if scoord in 'new2008':
        Cs = CSF(sc,theta_s,theta_b)
    elif scoord in 'old1994':
        Cs = (1.-theta_b)*cff1*np.sinh(theta_s*sc)+ \
           theta_b*(cff2*np.tanh(theta_s*(sc+0.5))-0.5)

    if scoord in 'new2008':
        hinv = 1. / (abs(topo) + hc)
        cff = hc * sc
        cff1 = Cs
            
        if len(np.array(zeta).shape)>2:
            for t in range(zeta.shape[0]):
                zeta[t][zeta[t]<(Dcrit-topo)] = Dcrit-topo[zeta[t]<(Dcrit-topo)]
                if 'w' in point_type:
                    z[t,0] = -topo
                    start = 1
                else:
                    start = 0
                for k in np.arange(start,N, dtype=int):
                    z[t,k] = zeta[t]+(zeta[t]+topo)* \
                             (cff[k]+cff1[k]*abs(topo))*hinv
        else:
            zeta[zeta<(Dcrit-topo)] = Dcrit-topo[zeta<(Dcrit-topo)]
            for k in np.arange(N, dtype=int):
                z[k] = zeta+(zeta+topo)* \
                       (cff[k]+cff1[k]*abs(topo))*hinv

    elif scoord in 'old1994':
        topo[topo==0] = 1.e-2
        hinv = 1./topo
        cff = hc*(sc-Cs)
        cff1 = Cs
        cff2 = sc + 1

        if len(np.array(zeta).shape)>2:
            for t in range(zeta.shape[0]):
                zeta[t][zeta[t]<(Dcrit-topo)] = Dcrit-topo[zeta[t]<(Dcrit-topo)]
                for k in np.arange(N,dtype=int) + 1:
                    z0 = cff[k-1]+cff1[k-1]*topo
                    z[t,k-1, :] = z0+zeta[t,:]*(1.+z0*hinv)
        else:
            zeta[zeta<(Dcrit-topo)] = Dcrit-topo[zeta<(Dcrit-topo)]
            for k in np.arange(N,dtype=int) + 1:
                z0 = cff[k-1]+cff1[k-1]*topo
                z[k-1, :] = z0+zeta*(1.+z0*hinv)
    else:
        raise Exception("Unknown scoord, should be 'new2008' or 'old1994'")

    if sc_r is None:
        sc_r = sc_r
    return z.squeeze(), np.float32(Cs)

####################

def ztosigma(vin,Z,zcroco):
    '''
    This fonction perform the z to sigma transformation for
    3D (Z,Y,X) or 4D (T,Z,Y,X) variables

    Input:
      Vin       Input variables to put on sigma grid (3 or 4D)
      Z         Input depth values (1D)
      zcroco    CROCO vertical level (3D)

    output:
       vout     Vin projection on zcroco
    '''
# Do a vertical interpolation from z levels to sigma CROCO levels
    if len(zcroco.shape)>3:
        [T,N,M,L]=np.shape(zcroco)
        four_dim=True
    else:
        [N,M,L]=np.shape(zcroco)
        four_dim=False

    [Nz]=np.shape(Z)
#
# Find the grid position of the nearest vertical levels
#
    i1=np.arange(0,L)
    j1=np.arange(0,M)
    if four_dim:
        t1=np.arange(0,T)
        [jmat,tmat,imat]=np.meshgrid(j1,t1,i1)
        VAR=np.reshape(vin,T*Nz*M*L)
        vout=np.zeros((T,N,M,L))
    else:
        [imat,jmat]=np.meshgrid(i1,j1)
        VAR=np.reshape(vin,Nz*M*L)
        vout=np.zeros((N,M,L))

    for ks in progressbar(range(N),' Sigma layer : ', 40):
        if four_dim:
            sigmalev=zcroco[:,ks,:,:]
            thezlevs=np.zeros((T,M,L),dtype=int)-1
        else:
            sigmalev=zcroco[ks,:,:]
            thezlevs=np.zeros((M,L),dtype=int)-1

        for kz in range(Nz):
            thezlevs[np.where(sigmalev>Z[kz])]=thezlevs[np.where(sigmalev>Z[kz])]+1

        if four_dim:
            pos= L*M*Nz*tmat+ L*M*thezlevs + L*jmat + imat
        else:
            pos= L*M*thezlevs + L*jmat + imat

        z1=Z[thezlevs]
        z2=Z[thezlevs+1]
        v1=VAR[pos]
        v2=VAR[pos+L*M]

        if four_dim:
            vout[:,ks,:,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
        else:
            vout[ks,:,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
    return vout

####################

#def vinterp(var, depths, z_r, z_w=None, mask=None,imin=0,jmin=0,kmin=1,interp_sfc=1,interp_bot=0,below=None,bounded=False,**kwargs):

#    '''
#    Interpolate a variable on a level of constant value. This level can be depths, isopycnal,... it only depends
#    on you choice for depths, z_r parameter.
#    ex:
#        1- If you want to have the temperature at 500-m depth
#            t500 = vinterp(temp,[-500],z_r) with z_r CROCO vertical coordinate
#        2- If you want to have the vorticity on the 27.7 isopycnal
#            vort27 = vinterp(vort,[27.7],rhop) with rhop the surface-referenced potential density 
#
#    This function use fortran subroutine to run faster
#
#    Inputs:
#      var           Variable to interpolate
#      depths        Levels on perform the interpolation (dim = k)
#      z_r           Vertical coordinates of the desired reference field (dim: N,M,L)
#      z_w           (Default None) CROCO z_w coordinates, if you use a field (temp,salt,...) leave empty
#      mask          (Default None) 2D mask
#      imin          (Default 0) minimum index in x-dimension
#      jmin          (Default 0) minimum index in y-dimension
#      kmin          (Default 1) minimum index in z-dimension 
#      interp_sfc    (Default 1) if 1 no interpolation below ground
#      interp_bot    (Default 0) if 1 data will be interpolated below ground (have priority over interp_sfc)
#      below         (Default None) Define a below value. When leave empty it is replace by 10000
#      bounded       (Default False) if True means that data above surface take surface value and data below bottom take bottom value
#
#    Outputs
#      vnew          Values interpolated on the desired levels (dim: k,M,L)
#    '''
#
#    if mask is None:  mask = np.ones((z_r.shape[-2],z_r.shape[-1])); mask[z_r[-1,:,:]==0] = 0
#
#    if z_w is None:
#        print('no z_w specified')
#        z_w=np.zeros((z_r.shape[0]+1,z_r.shape[1],z_r.shape[2]))
#        z_w[1:-1:,:,:] = 0.5*(z_r[1:,:,:] + z_r[:-1,:,:])
#        z_w[0,:,:] = z_r[0,:,:] - (z_r[1,:,:]-z_r[0,:,:])
#        z_w[-1,:,:] = z_r[-1,:,:] + (z_r[-1,:,:]-z_r[-2,:,:])
#
#    if np.ndim(depths)==1: newz = np.zeros((len(depths),var.shape[-2],var.shape[-1])) + depths
#    else: newz = depths
#
#    if bounded:
#        vnew=toolsf.sigma_to_z_intr_bounded(z_r, z_w,mask,var,newz,imin,jmin,kmin,9999.)
#    elif interp_bot==1:
#        if below is None: below=10000.
#        vnew=toolsf.sigma_to_z_intr_bot(z_r, z_w,mask,var,newz,below,imin,jmin,kmin,9999.)
#    elif interp_sfc==1:
#        vnew=toolsf.sigma_to_z_intr_sfc(z_r, z_w,mask,var,newz,imin,jmin,kmin,9999.)
#    else:
#        vnew=toolsf.sigma_to_z_intr(z_r, z_w,mask,var,newz,imin,jmin,kmin,9999.)
#
#    vnew[np.abs(vnew)==9999.]=np.nan
#
#    return np.squeeze(vnew)
#
########
#def vinterp_2d(var, depths, z_r, z_w=None, mask=None,imin=0,jmin=0,kmin=1, below=None,**kwargs):
#    ''' 
#    Same as vinterp but for 2D var
#
#    Inputs:
#      var           Variable to interpolate
#      depths        Levels on perform the interpolation (dim = k)
#      z_r           Vertical coordinates of the desired reference field (dim: N,L)
#      z_w           (Default None) CROCO z_w coordinates, if you use a field (temp,salt,...) leave empty
#      mask          (Default None) 2D mask
#      imin          (Default 0) minimum index in x-dimension
#      jmin          (Default 0) minimum index in y-dimension
#      kmin          (Default 1) minimum index in z-dimension 
#      below         (Default None) Define a below value. When leave empty it is replace by 0
#
#    Outputs
#      vnew          Values interpolated on the desired levels (dim: k,L)
#    '''
#    
#  
#    if mask is None:  mask = np.ones((var.shape[1])); mask[z_r[-1,:]==0] = 0
#
#    if z_w is None:
#        print('no z_w specified')
#        z_w=np.zeros((z_r.shape[0],z_r.shape[1]+1))
#        z_w[1:-1,:] = 0.5*(z_r[1:,:] + z_r[:-1,:])
#        z_w[0,:] = z_r[0,:] - (z_r[1,:]-z_r[0,:])
#        z_w[-1,:] = z_r[-1,:] + (z_r[-1,:]-z_r[-2,:])
#
#    if np.ndim(depths)==1: newz = np.zeros((len(depths),z_r.shape[1])) + depths
#    else: newz = depths
#
#    if below is None: below=0.
#    vnew=toolsf.sigma_to_z_intr_bot_2d(z_r, z_w,mask,var,newz,below,imin,kmin,9999.)
#
#    vnew[np.abs(vnew)==9999.]=np.nan
#
#    return np.squeeze(vnew)



####################

def vintegr(var,zw,zr,z01,z02):
    '''
    Vertically integrate a CROCO variable (var) from a constant depth
    z01 (ex z01=-4000 m) to a constant depth z02 (ex z02=-2000m).

    If z01 = NaN : perform the integration from the bottom to z02.
    If z02 = NaN : perform the integration from z01 to the surface.
    If they are both NaNs perform the integration from the bottom
    to the surface.

    Input :
     var         CROCO variable at RHO-points (3D matrix)
     zw          Depth of the W-points (3D matrix)
     zr          Depth of the RHO-points (3D matrix)
     z01         lower limit of integration (scalar)
     z02         upper limit of integration (scalar)

    Outputs :
    V            intgrated value (2D matrix)
    h0           layer thickness (2D matrix)

    Pierrick Penven 2005
    '''

    if z02 <= z01:
        print('vintegr2:  z02 <= z01')
        sys.exit()

    [Np,M,L]=np.shape(zw)
    N=Np-1
    mask=np.zeros([M,L]) + 1.

    i1=np.arange(0,L)
    j1=np.arange(0,M)
    [i2,j2]=np.meshgrid(i1,j1)

    za=np.reshape(zw,Np*M*L)
    vara=np.reshape(var,N*M*L)
#
# Get the matrix of the variables above of z01 and below z02
# 
    if (np.isfinite(z01) & np.isfinite(z02)):
        isgood=np.int_(zw[1:,:,:]>z01) * np.int_(zw[1:,:,:]<z02)
    elif np.isfinite(z01):
        isgood=np.int_(zw[1:,:,:]>z01)
    elif np.isfinite(z02):
        isgood=np.int_(zw[1:,:,:]<z02)
    else:
        isgood=np.int_(var==var)
#
    if np.isfinite(z01):
#
# Find the bottom limit of the corresponding grid cells
#
        a=np.int_(zw<z01)
        levs=np.sum(a,axis=0)-1
        mask=np.zeros([M,L]) + 1.
        mask[np.where(levs<0)]=np.nan
        mask[np.where( (levs<0) | (levs==N) )]=np.nan
        levs[np.where(levs==N)]=1

        pos = L*M*levs + L*j2 + i2
        z1=mask*za[pos]
#
# Compute the vertical size of the partial step to add at the bottom
#
        dzbot=z1-z01
        dzbot[np.where(np.isnan(mask))]=0.
#
# Get the value of the variable in the partial step to add at the bottom
#
        Vbot=vara[pos]
    else:
        dzbot=0
        Vbot=0

    if np.isfinite(z02):
#
# Find the top positions
#
        a=np.int_(zw<z02)
        levs=np.sum(a,axis=0)-1
        mask=np.zeros([M,L]) + 1.
        mask[np.where( (levs<0) | (levs==N) )]=np.nan
        levs[np.where(levs==N)]=1

        pos = L*M*levs + L*j2 + i2
        z2=mask*za[pos]
#
# Compute the vertical size of the partial step to add at the top
#
        dztop=z02-z2
        dztop[np.where(np.isnan(mask))]=0.
#
# Get the value of the variable in the partial step to add at the bottom
#
        Vtop=vara[pos]
    else:
        dztop=0
        Vtop=0
#
# Perform the vertical integration
    dz=zw[1:,:,:]-zw[:-1,:,:]
    V=np.sum(dz*isgood*var,axis=0) + dzbot*Vbot + dztop*Vtop
#
# Get the depth
#
    h0=np.sum(dz*isgood,axis=0) + dzbot + dztop

    V[np.where(h0==0)]=0
    h0[np.where(h0==0)]=0

    return V,h0


def vintegr4D(var,zw,zr,z01,z02):
    '''
    Do the samething than vintegr but for 4D vars
    zw and zr need to be 4D too
    '''
    if z02 <= z01:
        print('vintegr2:  z02 <= z01')
        sys.exit()

    [T,Np,M,L]=np.shape(zw)
    N=Np-1
    mask=np.zeros([T,M,L]) + 1.

    i1=np.arange(0,L)
    j1=np.arange(0,M)
    t1=np.arange(0,T)
    [j2,t2,i2]=np.meshgrid(j1,t1,i1)
#    [i2,j2]=np.meshgrid(i1,j1)

    za=np.reshape(zw,T*Np*M*L)
    vara=np.reshape(var,T*N*M*L)
#
# Get the matrix of the variables above of z01 and below z02
#
    if (np.isfinite(z01) & np.isfinite(z02)):
        isgood=np.int_(zw[:,1:,:,:]>z01) * np.int_(zw[:,1:,:,:]<z02)
    elif np.isfinite(z01):
        isgood=np.int_(zw[:,1:,:,:]>z01)
    elif np.isfinite(z02):
        isgood=np.int_(zw[:,1:,:,:]<z02)
    else:
        isgood=np.int_(var==var)
#
    if np.isfinite(z01):
#
# Find the bottom limit of the corresponding grid cells
#
        a=np.int_(zw<z01)
        levs=np.sum(a,axis=1)-1
        mask=np.zeros([T,M,L]) + 1.
        mask[np.where(levs<0)]=np.nan
        mask[np.where( (levs<0) | (levs==N) )]=np.nan
        levs[np.where(levs==N)]=1

#        pos = L*M*levs + L*j2 + i2
        pos= L*M*T*levs + L*M*j2 + L*j2 + i2
        z1=mask*za[pos]
#
# Compute the vertical size of the partial step to add at the bottom
#
        dzbot=z1-z01
        dzbot[np.where(np.isnan(mask))]=0.
#
# Get the value of the variable in the partial step to add at the bottom
#
        Vbot=vara[pos]
    else:
        dzbot=0
        Vbot=0

    if np.isfinite(z02):
#
# Find the top positions
#
        a=np.int_(zw<z02)
        levs=np.sum(a,axis=1)-1
        mask=np.zeros([T,M,L]) + 1.
        mask[np.where( (levs<0) | (levs==N) )]=np.nan
        levs[np.where(levs==N)]=1

#        pos = L*M*levs + L*j2 + i2
        pos= L*M*T*levs + L*M*j2 + L*j2 + i2
        z2=mask*za[pos]
#
# Compute the vertical size of the partial step to add at the top
#
        dztop=z02-z2
        dztop[np.where(np.isnan(mask))]=0.
#
# Get the value of the variable in the partial step to add at the bottom
#
        Vtop=vara[pos]
    else:
        dztop=0
        Vtop=0
#
# Perform the vertical integration
    dz=zw[:,1:,:,:]-zw[:,:-1,:,:]
    V=np.sum(dz*isgood*var,axis=1) + dzbot*Vbot + dztop*Vtop
#
# Get the depth
#
    h0=np.sum(dz*isgood,axis=1) + dzbot + dztop

    V[np.where(h0==0)]=0
    h0[np.where(h0==0)]=0

    return V,h0

