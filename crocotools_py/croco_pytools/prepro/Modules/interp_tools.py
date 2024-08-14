import numpy as np
from scipy.interpolate import griddata
import scipy.interpolate as itp
import os
import sys
import time
from scipy.spatial import Delaunay
import Cgrid_transformation_tools as grd_tools
from sigmagrid_tools import ztosigma
from progressbar import progressbar
import xarray as xr
import pyinterp
import pyinterp.backends.xarray as pyxr
import pyinterp.fill as pyfi
#



def make_xarray(data,lon2D,lat2D):
    data_array= xr.DataArray(data,
               dims=["lat", "lon"],
               coords=[lat2D[:,0].astype(np.float64),lon2D[0,:].astype(np.float64)])
    data_array.lon.attrs['units']='degrees_east'
    data_array.lat.attrs['units']='degrees_north'
    return data_array

###############
# Get interpolation weight
###############
def get_tri_coef(X, Y, newX, newY, verbose=0):

    """
    Inputs:
        origin lon and lat 2d arrays (X,Y)
        child lon and lat 2d arrays (newX,newY)
    Ouputs:
        elem - pointers to 2d gridded data (at lonp,latp locations) from
            which the interpolation is computed (3 for each child point)
        coef - linear interpolation coefficients
    Use:
        To subsequently interpolate data from Fp to Fc, the following
        will work:      Fc  = sum(coef.*Fp(elem),3);  This line  should come in place of all
        griddata calls. Since it avoids repeated triangulations and tsearches (that are done
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

    points = tri.points[tri.vertices[tri.find_simplex(Xc)]]
    if verbose==1: tstart = tm.time()
    for i in progressbar(range(npts),'  Get_tri_coef: ', 40):

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
    elem = np.reshape(tri.vertices[tri.find_simplex(Xc)],(newX.shape[0],newY.shape[1],3))
    coef = np.reshape(p,(newX.shape[0],newY.shape[1],3))
    if verbose==1: print('Coef. computation 2', tm.time()-tstart)

    return(elem,coef)

#######################################

def get_delaunay_bry(lon_bry,lat_bry,inputfile,bdy):
    '''
    This function computes the delaunay matrices for the interpolations 
    at the boundaies

    Input:
      lon_bry      Longitudes of the boundary (vector).
      lat_bry      Latitudes of the boundary (vector).
      inputfile    netcdf structure poiting to the input file
      bdy          which boundary is done

    Output:
      LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
      LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
      LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    ''' 
    comp_delaunay=1

#
# get grid positions
#
    LonU_bry,LatU_bry   = eval(''.join(("inputfile.lonU"+bdy))),eval(''.join(("inputfile.latU"+bdy)))
    LonV_bry,LatV_bry   = eval(''.join(("inputfile.lonV"+bdy))),eval(''.join(("inputfile.latV"+bdy)))
    LonT_bry,LatT_bry   = eval(''.join(("inputfile.lonT"+bdy))),eval(''.join(("inputfile.latT"+bdy))) 

#
# Get the 2D interpolation coefficients
#

    if comp_delaunay==1:

      print('\nCompute Delaunay triangulation from T-points to CROCO rho_points...')
      print('-------------------------------------------------------------------')
      [elemT_bry,coefT_bry] = get_tri_coef(LonT_bry,LatT_bry,lon_bry,lat_bry)
      coefnorm=np.sum(coefT_bry,axis=2)
      coefT_bry=coefT_bry/coefnorm[:,:,np.newaxis]

      print('\nCompute Delaunay triangulation from U-points to CROCO rho_points...')
      print('-------------------------------------------------------------------')
      [elemU_bry,coefU_bry] = get_tri_coef(LonU_bry,LatU_bry,lon_bry,lat_bry)
      coefnorm=np.sum(coefU_bry,axis=2)
      coefU_bry=coefU_bry/coefnorm[:,:,np.newaxis]

      print('\nCompute Delaunay triangulation from V-points to CROCO rho_points...')
      print('-------------------------------------------------------------------')
      [elemV_bry,coefV_bry] = get_tri_coef(LonV_bry,LatV_bry,lon_bry,lat_bry)
      coefnorm=np.sum(coefV_bry,axis=2)
      coefV_bry=coefV_bry/coefnorm[:,:,np.newaxis]

      # Save the Delaunay triangulation matrices
      np.savez('coeffs_bry'+bdy+'.npz',\
                coefT=coefT_bry,elemT=elemT_bry,\
                coefU=coefU_bry,elemU=elemU_bry,\
                coefV=coefV_bry,elemV=elemV_bry)

    else:
#
# Load the Delaunay triangulation matrices
#
      print('Load Delaunay triangulation...')
      data=np.load('coeffs_bry'+bdy+'.npz')
      coefT_bry = data['coefT_bry']
      elemT_bry = data['elemT_bry']
      coefU_bry = data['coefU_bry']
      elemU_bry = data['elemU_bry']
      coefV_bry = data['coefV_bry']
      elemV_bry = data['elemV_bry']

    print('Delaunay triangulation done')

    return(elemT_bry,coefT_bry,elemU_bry,coefU_bry,elemV_bry,coefV_bry)



####################################################################
def add2layers(vin):
    '''
    Add a layer below the bottom and above the surface to avoid
    vertical extrapolations when doing a vertical interpolation

    Input:
      Vin    3 or 4D Variable 

    Output:
      vout   Vin with 2 new layers above and below
    '''
    if len(np.shape(vin))==3:
        [Nz,M,L]=np.shape(vin)
        vout=np.zeros((Nz+2,M,L))

        vout[1:-1,:,:]=vin
        vout[0,:,:]=vin[0,:]
        vout[-1,:,:]=vin[-1,:]

    elif len(np.shape(vin))==4:
        [T,Nz,M,L]=np.shape(vin)
        vout=np.zeros((T,Nz+2,M,L))

        vout[:,1:-1,:,:]=vin
        vout[:,0,:,:]=vin[:,0,:]
        vout[:,-1,:,:]=vin[:,-1,:]
    return vout

######################################################################

def interp_tracers(inputfile,vname,k,crocogrd,dtmin,dtmax,prev=0,nxt=0,bdy=""):
    '''
    Remove the missing values from a gridded 2D field
    and do an horizontal interpolation using pyinterp tools

    Inputs:
      inputfile     Input class containing useful tools (see input_class.py)
      vname         Variable's name to interpolate
      k             Depth index ( -1 when no depth component)
      crocogrd      CROCO grid class
      dtmin         Starting index of the time series
      dtmax         Ending index of the time series
      prev          If 1 duplicates the first index of the time series
      nxt           If 1 duplicates the last index of the time series
      bdy           Specify which boundary you are on. Leave empty if not

    Output:
      vout          Interpolation of the field time series on 2D CROCO grid
      Nzgood        Number of good points on the input level k        
    '''

# 0: Read input data informations
    nc        = inputfile.ncglo
    varinp    = inputfile.var
    if vname in [ 'u','ubar']:
        Lon,Lat = eval(''.join(("inputfile.lonU"+bdy))),eval(''.join(("inputfile.latU"+bdy)))
    elif  vname in [ 'v','vbar']:
        Lon,Lat = eval(''.join(("inputfile.lonV"+bdy))),eval(''.join(("inputfile.latV"+bdy)))
    else:
        Lon,Lat = eval(''.join(("inputfile.lonT"+bdy))),eval(''.join(("inputfile.latT"+bdy)))

# 1: Read data
    if dtmin != dtmax:
        l=np.arange(dtmin,dtmax+1)
    else:
        l=dtmin

    Vin = inputfile.var_periodicity(vname,l,k,bdy=bdy)

    if dtmin==dtmax:
        Vin=Vin[np.newaxis,:]

    if prev == 1: # create overlap before when no data avail
        Vtmp=np.zeros([Vin.shape[0]+1,Vin.shape[1],Vin.shape[2]])
        Vtmp[1:,:]=np.copy(Vin)
        Vtmp[0,:]=np.copy(Vin[0,:])
        Vin=Vtmp
        del Vtmp
    if nxt==1: # create overlap after when no data avail
        Vtmp=np.zeros([Vin.shape[0]+1,Vin.shape[1],Vin.shape[2]])
        Vtmp[:-1,:]=np.copy(Vin)
        Vtmp[-1,:]=np.copy(Vin[-1,:])
        Vin=Vtmp
        del Vtmp

# 2: Remove bad values (using nearest values)
    if "_FillValue" not in inputfile.ncglo[vname].encoding:# If no FillValue in netcdf, assume 0 as value for the mask  
        Vin[Vin==0]=np.nan 

    if bdy != "":
        if bdy == 'S':
            bound='south'
        elif bdy == 'W':
            bound='west'
        elif bdy == 'E':
            bound='east'
        elif bdy == 'N':
            bound='north'
        crocolon=eval(''.join(('crocogrd.lon_',bound)))
        crocolat=eval(''.join(('crocogrd.lat_',bound)))
    else:
        crocolon=crocogrd.lon
        crocolat=crocogrd.lat

    Vout=np.zeros([Vin.shape[0],crocolon.shape[0],crocolon.shape[1]])
   
    for tt in range(Vin.shape[0]):
        igood = np.where(np.isnan(Vin[tt,:])==False)
        ibad  = np.where(np.isnan(Vin[tt,:]))
 
        NzGood=np.size(igood)
        Nbad=np.size(ibad)

        if NzGood==0:
#            print('\nWarning: no good data')
            Vin[:]=np.nan
        elif NzGood<10:
#            print('\nWarning: less than 10 good values')
            Vinfilled=np.copy(Vin[tt,:])
            Vinfilled[:]=np.nanmean(Vinfilled[igood])
            val_interpolator=pyinterp.backends.xarray.Grid2D(make_xarray(Vinfilled,Lon,Lat))           
        else:
            spline = itp.NearestNDInterpolator((Lon[igood].ravel(),Lat[igood].ravel()),Vin[tt,igood[0],igood[1]])
            Vinfilled =np.copy(np.squeeze(Vin[tt,:]))
            Vinfilled[ibad] = spline(Lon[ibad[0],ibad[1]],Lat[ibad[0],ibad[1]])
            val_interpolator=pyinterp.backends.xarray.Grid2D(make_xarray(Vinfilled,Lon,Lat))
# 3: 2D interpolation
        if NzGood ==0:
            Vout[tt,:]=np.nan
        else:
            Vout[tt,:]  = val_interpolator.bivariate(coords=dict(lon=crocolon.flatten(),lat=crocolat.flatten()),num_threads=1).reshape(crocolon.shape)
    return Vout,NzGood


#######################################################

def interp(inputfile,vname,Nzgoodmin,z_rho,crocogrd,dtmin,dtmax,prev=0,nxt=0,bdy=""):
    '''
    Do a full interpolation of a 3-4d variable from z-grid data to a CROCO sigma grid
    1 - Horizontal interpolation on each z-levels
    2 - Vertical Interpolation from z to CROCO sigma levels

    Inputs:
      inputfile     Input class containing useful tools (see input_class.py)
      vname         Variable's name to interpolate
      Nzgoodmin     Number of value to consider a z-level fine to be used
      z_rho         CROCO vertical levels
      crocogrd      CROCO grid class
      dtmin         Starting index of the time series
      dtmax         Ending index of the time series
      prev          If 1 duplicates the first index of the time series
      nxt           If 1 duplicates the last index of the time series
      bdy           Specify which boundary you are on. Leave empty if not

    Output:
      vout          3-4D interpolation of vname variable
    '''

    if np.ndim(z_rho)==3:
        z_rho=z_rho[np.newaxis,:]
    [T,N,M,L]=np.shape(z_rho)
    depth= inputfile.depth
    [Nz]=np.shape(depth)

    comp_horzinterp = 1 # if 1 compute horizontal interpolations - 0 use saved matrices (for debugging)
    if comp_horzinterp==1:
        print('Horizontal interpolation over z levels')
        t4d=np.zeros((T,Nz,M,L))
        kgood=-1
        for k in progressbar(range(Nz),vname+': ', 40):#range(Nz)
            (t3d,Nzgood) = interp_tracers(inputfile,vname,k,crocogrd,dtmin,dtmax,prev,nxt,bdy=bdy)
            if Nzgood>Nzgoodmin:
                kgood=kgood+1
                t4d[:,kgood,:,:]=t3d
                
        t4d=t4d[:,0:kgood,:,:]
        depth=depth[0:kgood]
        np.savez('t4d.npz',t4d=t4d,depth=depth)

    else:
        print('Load matrix...')
        data=np.load('t4d.npz')
        t3d = data['t4d']
        depth = data['depth']

    [Nz]=np.shape(depth)
    if depth[0]<0:
        Z=depth
    else:
        Z=-depth
    
#  Vertical interpolation
    print('Vertical interpolation')

# Add a layer below the bottom and above the surface to avoid vertical extrapolations
# and flip the matrices upside down (Z[Nz]=surface)

    Z=np.flipud(np.concatenate(([100.],Z,[-10000.])))
    [Nz]=np.shape(Z)
    t4d=np.flip(add2layers(t4d),axis=1)

# Do the vertical interpolations
    vout=ztosigma(t4d,Z,z_rho)
    return vout

####################################################################


def interp_uv(inputfile,Nzgoodmin,z_rho,cosa,sina,\
        crocogrd,dtmin,dtmax,prev=0,nxt=0,bdy=""):
    '''
    Same as interp but for horizontal velocities u,v

    Inputs:
      inputfile     Input class containing useful tools (see input_class.py)
      Nzgoodmin     Number of value to consider a z-level fine to be used
      z_rho         CROCO vertical levels
      cosa          Cosine value of grid angle
      sina          Sinus value of grid angle
      crocogrd      CROCO grid class
      dtmin         Starting index of the time series
      dtmax         Ending index of the time series
      prev          If 1 duplicates the first index of the time series
      nxt           If 1 duplicates the last index of the time series
      bdy           Specify which boundary you are on. Leave empty if not

    Outputs:
       uout         U velocity on 4D (Time,CROCO-Ugrid)
       vout         V velocity on 4D (Time,CROCO-Vgrid)
       ubar         Integrated U velocity on (Time,CROCO-Ugrid)
       vbar         Integrated V velocity on (Time,CROCO-Vgrid)
    '''

    if np.ndim(z_rho)==3:
        z_rho=z_rho[np.newaxis,:]
    [T,N,M,L]=np.shape(z_rho)

    cosa3d=np.tile(cosa, (T, 1, 1))
    sina3d=np.tile(sina, (T, 1, 1))

    depth=inputfile.depth
    dz=np.gradient(depth)
    [Nz]=np.shape(depth)
    comp_horzinterp = 1 # if 1 compute horizontal interpolations - 0 use saved matrices (for debugging)
    if comp_horzinterp==1:
        print('Horizontal interpolation of u and v over z levels')
        u4d=np.zeros((T,Nz,M,L-1))
        v4d=np.zeros((T,Nz,M-1,L))
        ubar=np.zeros((T,M,L-1))
        vbar=np.zeros((T,M-1,L))
        zu  =ubar
        zv  =vbar
        kgood=-1
        for k in progressbar(range(Nz),' uv : ', 40):
            (u3d,Nzgood_u) = interp_tracers(inputfile,'u',k,crocogrd,dtmin,dtmax,prev,nxt,bdy=bdy)
            (v3d,Nzgood_v) = interp_tracers(inputfile,'v',k,crocogrd,dtmin,dtmax,prev,nxt,bdy=bdy)
            Nzgood=np.min((Nzgood_u,Nzgood_v))
            if Nzgood>Nzgoodmin:
                kgood=kgood+1
#
# Rotation and put to u-points and v-points 
#
                u4d[:,kgood,:,:]=grd_tools.rho2u(u3d*cosa3d+v3d*sina3d)
                v4d[:,kgood,:,:]=grd_tools.rho2v(v3d*cosa3d-u3d*sina3d)

                ubar = ubar + grd_tools.rho2u((u3d*dz[kgood])*cosa3d+(v3d*dz[kgood])*sina3d)
                zu   = zu   + dz[kgood]*np.ones(ubar.shape)
                vbar = vbar + grd_tools.rho2v((v3d*dz[kgood])*cosa3d-(u3d*dz[kgood])*sina3d)
                zv   = zv   + dz[kgood]*np.ones(vbar.shape)


        u4d=u4d[:,0:kgood,:,:]
        v4d=v4d[:,0:kgood,:,:]
        ubar=ubar/zu
        vbar=vbar/zv
        depth=depth[0:kgood]
        np.savez('u4d.npz',u4d=u4d,v4d=v4d,depth=depth,ubar=ubar,vbar=vbar)

    else:
        print('Load matrices...')
        data=np.load('u4d.npz')
        u4d = data['u4d']
        v4d = data['v4d']
        depth = data['depth']
        ubar = data['ubar']
        vbar = data['vbar']

    [Nz]=np.shape(depth)
    if depth[0]<0:
        Z=depth
    else:
        Z=-depth

#--------------------------------------------------
#  Vertical interpolation
#----------------------------------------------------
#
    print('Vertical interpolation')
#
# Add a layer below the bottom and above the surface to avoid vertical extrapolations
# and flip the matrices upside down (Z[Nz]=surface)
#
    Z=np.flipud(np.concatenate(([100.],Z,[-10000.])))
    [Nz]=np.shape(Z)
    u4d=np.flip(add2layers(u4d),axis=1)
    v4d=np.flip(add2layers(v4d),axis=1)
#
# Do the vertical interpolations 
#
    uout=ztosigma(u4d,Z,grd_tools.rho2u(z_rho))
    vout=ztosigma(v4d,Z,grd_tools.rho2v(z_rho))

    return uout,vout,ubar,vbar




def interp_tides(inputfile,vname,k,crocogrd,dtmin,dtmax,input_type,prev=0,nxt=0,bdy=""):

    '''
    Compute Extrapolate data with pyinterp gauss_seidel and perform horizontal interpolation for tides

    Inputs:
      inputfile     Input class containing useful tools (see input_class.py)
      vname         Variable's name to interpolate
      Nzgoodmin     Number of value to consider a z-level fine to be used
      z_rho         CROCO vertical levels
      coef          Coef from Delaunay triangulation
      elem          Elem from Delaunay triangulation
      dtmin         Starting index of the time series
      dtmax         Ending index of the time series
      prev          If 1 duplicates the first index of the time series
      nxt           If 1 duplicates the last index of the time series
      bdy           Specify which boundary you are on. Leave empty if not

    Output:
      vout          4D interpolation of vname variable (Time+ CROCO-grid)
    '''
# 0: Read input data informations
    nc        = inputfile.ncglo
    varinp    = inputfile.var

    if vname == 'H':
        var=['ssh']
    elif vname == 'cur':
        var=['u','v']

    for field in var:
        if field == 'u':
            Lon,Lat = eval(''.join(("inputfile.lonU"+bdy))),eval(''.join(("inputfile.latU"+bdy)))
        elif  field == 'v':
            Lon,Lat = eval(''.join(("inputfile.lonV"+bdy))),eval(''.join(("inputfile.latV"+bdy)))
        else:
            Lon,Lat = eval(''.join(("inputfile.lonT"+bdy))),eval(''.join(("inputfile.latT"+bdy)))

# 1: Read data
        if dtmin != dtmax:
            l=np.arange(dtmin,dtmax+1)
        else:
            l=dtmin

        part1 = inputfile.var_periodicity(field+'_part1',l,k,bdy=bdy)
        part2 = inputfile.var_periodicity(field+'_part2',l,k,bdy=bdy)
        
        if input_type=='Amp_Phase':
            rpart,impart = part1*np.cos(part2),part1*np.sin(part2)
        else:
            rpart,impart = part1,part2

        if dtmin==dtmax:
            rpart=rpart[np.newaxis,:]
            impart=impart[np.newaxis,:]

# 2: Remove bad values (using nearest values)
#        if "_FillValue" not in inputfile.ncglo[field+'_part1'].encoding:# If no FillValue in netcdf, assume 0 as value for the mask  
#            rpart[rpart==0]=np.nan
#            impart[impart==0]=np.nan

        if bdy != "":
            if bdy == 'S':
                bound='south'
            elif bdy == 'W':
                bound='west'
            elif bdy == 'E':
                bound='east'
            elif bdy == 'N':
                bound='north'
            crocolon=eval(''.join(('crocogrd.lon_',bound)))
            crocolat=eval(''.join(('crocogrd.lat_',bound)))
        else:
            crocolon=crocogrd.lon
            crocolat=crocogrd.lat

        Rout=np.zeros([rpart.shape[0],crocolon.shape[0],crocolon.shape[1]])
        Iout=np.zeros([impart.shape[0],crocolon.shape[0],crocolon.shape[1]])

        for tt in range(rpart.shape[0]):
            # assume rpart and impart have the same mask
            igood = np.where(np.isnan(rpart[tt,:])==False)
            ibad  = np.where(np.isnan(rpart[tt,:]))

            NzGood=np.size(igood)
            Nbad=np.size(ibad)

            if NzGood==0:
                rpart[:]=np.nan
                impart[:]=np.nan
            elif NzGood<10:
                rpartfilled=np.copy(rpart[tt,:])
                rpartfilled[:]=np.nanmean(rpartfilled[igood])
                r_val_interpolator=pyinterp.backends.xarray.Grid2D(make_xarray(rpartfilled,Lon,Lat))
                impartfilled=np.copy(impart[tt,:])
                impartfilled[:]=np.nanmean(impartfilled[igood])
                i_val_interpolator=pyinterp.backends.xarray.Grid2D(make_xarray(impartfilled,Lon,Lat))
            else:
                r_spline = itp.NearestNDInterpolator((Lon[igood].ravel(),Lat[igood].ravel()),rpart[tt,igood[0],igood[1]])
                rpartfilled =np.copy(np.squeeze(rpart[tt,:]))
                rpartfilled[ibad] = r_spline(Lon[ibad[0],ibad[1]],Lat[ibad[0],ibad[1]])
                r_val_interpolator=pyinterp.backends.xarray.Grid2D(make_xarray(rpartfilled,Lon,Lat))
                i_spline = itp.NearestNDInterpolator((Lon[igood].ravel(),Lat[igood].ravel()),impart[tt,igood[0],igood[1]])
                impartfilled =np.copy(np.squeeze(impart[tt,:]))
                impartfilled[ibad] = i_spline(Lon[ibad[0],ibad[1]],Lat[ibad[0],ibad[1]])
                i_val_interpolator=pyinterp.backends.xarray.Grid2D(make_xarray(impartfilled,Lon,Lat))
# 3: 2D interpolation
            if NzGood ==0:
                Rout[tt,:]=np.nan
                Iout[tt,:]=np.nan
            else:
                Rout[tt,:]  = r_val_interpolator.bivariate(coords=dict(lon=crocolon.flatten(),lat=crocolat.flatten()),num_threads=1).reshape(crocolon.shape)
                Iout[tt,:]  = i_val_interpolator.bivariate(coords=dict(lon=crocolon.flatten(),lat=crocolat.flatten()),num_threads=1).reshape(crocolon.shape)

        Vout=(Rout+1j*Iout)

 
        if field =='u':
            Uout=np.copy(Vout)


    if vname == 'H':
        return Vout,NzGood

    elif vname == 'cur':
        return Uout,Vout,NzGood

