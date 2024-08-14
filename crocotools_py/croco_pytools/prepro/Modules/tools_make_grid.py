import numpy as np
import netCDF4 as netcdf
import xarray as xr
import topo_reader
import toolsf
import netCDF4 as netcdf
from datetime import datetime
import scipy.interpolate as itp
from scipy.interpolate import griddata
import regionmask
import geopandas as gp
import pandas as pd
from shapely.geometry import Polygon
from itertools import product
import sys
#
def topo_periodicity(topo_file, geolim):
    '''
    topo_periodicity checks whether domain is inside the topo file.
    If so, check if there is a need for create a periodicity between
    the last and first longitude points ( for global data).
    It is returning lon/lat/topo adapted to the desired domain
    geolim = [lonmin,lonmax,latmin,latmax]
    '''

    topo_type = topo_reader.lookvar(topo_file)

    try:
        print('Reading topography file:', topo_file)
#        nc = netcdf.Dataset(topo_file)
        nc = xr.open_dataset(topo_file)
    except:
        sys.exit(''.join(('ERROR: \n', 'Topo file -> ',topo_file ,' does not exist... ')))
    
    topo_lon = eval(''.join(("nc."+topo_type['lon']+'.values')))
    topo_lat = eval(''.join(("nc."+topo_type['lat']+'.values')))
    if topo_lon.ndim==2: # gebco is a bit different
        topo_lon = np.linspace(topo_lon[0,0],
                               topo_lon[0,-1], num=topo_lon.shape[1])
        topo_lat = np.linspace(topo_lat[0,0],
                               topo_lat[-1,0], num=topo_lat.shape[0])
        gebco = True
    else:
        gebco = False

    if topo_type['zaxis']== 'down':
        topo_fact=-1
    else:
        topo_fact=1

    for i in range(1,topo_lon.shape[0]): # Fix etopo5 discontinuity
        if topo_lon[i]<topo_lon[i-1]:    # between 180/-180 in the
            topo_lon[i]=topo_lon[i]+360  # middle

####
    jmin=indx_bound(topo_lat, geolim[2])
    jmax=indx_bound(topo_lat, geolim[-1])
    if -1<jmin and jmin<topo_lat.shape[0] and \
       -1<jmax and jmax<topo_lat.shape[0] :
        if jmin > 0 :
            jmin=jmin-1
        jmax=jmax+2
    else:
        print('North-south extents of the dataset ',topo_lat[0],topo_lat[-1],' are not sufficient to cover the entire model grid.')
        exit()

    imin=indx_bound(topo_lon, geolim[0])
    imax=indx_bound(topo_lon, geolim[1])

    # Do not include border to check periodicity
    if -1<imin and imin<topo_lon.shape[0] and \
       -1<imax and imax<topo_lon.shape[0] :
        if imin>0:
            imin=imin-1
        imax=imax+1
        shft_west=0 ; shft_east=0
        print('Single region dataset imin/imax=',imin,imax, )
    else:
        ######
        ptest=topo_lon[-1]-topo_lon[0]-360
        dx=(topo_lon[-1]-topo_lon[0])/(topo_lon.shape[0]-1)
        epsil=0.01*abs(dx)
        if abs(ptest) < epsil :
            period=topo_lon.shape[0]-1
        elif abs(ptest+dx) < epsil :
            period=topo_lon.shape[0]
        else:
            period=0

        if period>0:
            print('Identified periodicity domain in data of ', period,' points out of', topo_lon.shape[0])
        else :
            print('ERROR: The data does not cover the entire grid. Change your grid definition')
            exit()
        ##
        shft_west=0
        if imin==-1 :
            shft_west=-1
            imin=indx_bound(topo_lon, geolim[0]+360)
            if imin == topo_lon.shape[0]: imin = topo_lon.shape[0]-1
        elif imin==topo_lon.shape[0] :
            shft_west=+1
            imin=indx_bound(topo_lon, geolim[0]-360)
            if imin == -1: imin = topo_lon.shape[0]-1
        ##
        shft_east=0
        if imax == -1:
            shft_east=-1
            imax=indx_bound(topo_lon, geolim[1]+360)
            if imax == topo_lon.shape[0]: imax = 0
        elif imax == topo_lon.shape[0]:
            shft_east=+1
            imax=indx_bound(topo_lon, geolim[1]-360)
            if imax == -1: imax = 0

        # Problem if index not in [0,lon.shape[0]-1]
        if -1<imin and imin<topo_lon.shape[0] and \
           -1<imax and imax<topo_lon.shape[0] :
            if imin>0:
                imin=imin-1
            imax=imax+1
        else:
            print('ERROR: Data longitude covers 360 degrees, but still cannot find  starting and ending indices.')
            exit()
    
    print('Bounding indices of the relevant part to be extracted from the entire dataset:\n', \
          'imin,imax =', imin,imax,'out of', topo_lon.shape[0],'jmin,jmax =',jmin,jmax, 'out of',topo_lat.shape[0])
    ny_lat=jmax-jmin+1
    start2=jmin ; end2=start2+ny_lat; count2=ny_lat
    lat_tmp=np.zeros([ny_lat])
    for j in range(0,ny_lat):
        lat_tmp[j]=topo_lat[j+jmin]
 
    #####

    if imin < imax :
        nx_lon=imax-imin+1
        start1=imin ; end1=start1+nx_lon ; count1=nx_lon
        if gebco:
            topo = topo_fact*eval(''.join(("nc."+topo_type['topo']+'.values')))
            topo = np.reshape(topo, (topo_lat.size, topo_lon.size))
            topo = topo[start2:end2, start1:end1]
        else:
            topo = topo_fact*eval(''.join(("nc."+topo_type['topo']+'[start2:end2, start1:end1]'+'.values')))
        nc.close()

        ishft=imin
        lon_tmp=np.zeros([nx_lon])
        if shft_west>0 and shft_east>0:
            for i in range(0,nx_lon):
                lon_tmp[i]=topo_lon[i+ishft] +360
        elif shft_west<0 and shft_east<0:
            for i in range(0,nx_lon):
                 lon_tmp[i]=topo_lon[i+ishft]-360
        elif shft_west== 0 and shft_east==0:
            for i in range(0,nx_lon) :
                lon_tmp[i]=topo_lon[i+ishft]
        else:
            print('Error in shifting algoritm')
            exit()

    elif imin>imax:
        print('Reading topography in two separate parts adjacent through 360-degree periodicity,' )
        print('first..., ')
        nx_lon=imax+period-imin+1
        htopo = np.zeros([count2,nx_lon])
        xtmp  = np.zeros([nx_lon])
        start1=0 ; end1=start1+imax+1; count1=imax+1

        if gebco:
            topo = topo_fact*eval(''.join(("nc."+topo_type['topo']+'.values')))
            topo = np.reshape(topo, (topo_lat.size, topo_lon.size))
            topo = topo[start2:end2, start1:end1]
        else:
            topo = topo_fact*eval(''.join(("nc."+topo_type['topo']+'[start2:end2, start1:end1]'+'.values')))
        for j in range(0,count2):
            for i in range(0,count1):
                htopo[j,nx_lon-imax+i-1]=topo[j,i]
        del topo

        ishft=nx_lon-count1
        if shft_east>0:
            for i in range(0,count1):
                xtmp[i+ishft]=topo_lon[i] +360
        elif shft_east<0:
            for i in range(0,count1):
                xtmp[i+ishft]=topo_lon[i] -360
        else:
            for i in range(0,count1):
                xtmp[i+ishft]=topo_lon[i]

        print('second...')
        start1=imin ; count1=period-imin; end1=start1+count1
        if gebco:
            topo = topo_fact*eval(''.join(("nc."+topo_type['topo']+'.values')))
            topo = np.reshape(topo, (topo_lat.size, topo_lon.size))
            topo = topo[start2:end2, start1:end1]
        else:
            topo = topo_fact*eval(''.join(("nc."+topo_type['topo']+'[start2:end2, start1:end1]'+'.values')))
        nc.close()

        for j in range(0,count2):
            for i in range(0,count1):
                htopo[j,i]=topo[j,i]
        del topo

        ishft=imin
        if shft_west>0:
            for i in range(0,count1):
                xtmp[i]=topo_lon[i+ishft] +360
        elif shft_west<0 :
            for i in range(0,count1):
                xtmp[i]=topo_lon[i+ishft] -360
        else:
            for i in range(0,count1):
                xtmp[i]=topo_lon[i+ishft]
        lon_tmp=np.zeros([xtmp.shape[0]])
        for i in range(0,nx_lon):
            lon_tmp[i]=xtmp[i]

        topo=np.copy(htopo)

    del topo_lon,topo_lat
   
    topo[np.isnan(topo)]=np.nanmax(topo.ravel())
    
    topo_lon=np.copy(lon_tmp)
    topo_lat=np.copy(lat_tmp)

    return topo_lon,topo_lat,topo




class GetTopo():
    """
    GetTopo object.  At present this class will identify and read nc files from:
      ETOPO
      GEBCO
      Romstools (etopo2.nc; http://www.romsagrif.org)
    
    The default topography is an Etopo5 file (etopo5.nc).
    A file selector is available to point to choose alternative topo products.
    """
    def nearest1d(self, point, array):
        '''
        Return index of nearest point in array
        '''
        return np.argmin(np.abs(array - point))

    def topo(self,outputs, topo_file,shpfile,smooth=None,hmin=None,hmax=None,sgl_connect=None,prt_grd=None,coef=None):

        if smooth is not None:
            rd=smooth.smthr
        else:
            rd = 1 

        topo_type=topo_reader.lookvar(topo_file)
        if 'srtm' in topo_type.keys():
            srtm_file='/'.join(topo_file.split('/')[:-1])
            topo=toolsf.srtopo(srtm_file,outputs.lon_rho,outputs.lat_rho,outputs.pm,outputs.pn,rd)
        else:          
            lonmin,lonmax,latmin,latmax=toolsf.roms_grid_geo_bounds(outputs.lon_rho,outputs.lat_rho,rd)
            topo_lon,topo_lat,topo=topo_periodicity(topo_file,[lonmin,lonmax,latmin,latmax])
            print('Interpolating topography to CROCO grid')
            topo=toolsf.compute_hraw(topo_lon,topo_lat,topo.T,outputs.lon_rho,outputs.lat_rho,outputs.pm,outputs.pn,rd)
            print('Finished interpolating')

        if smooth is not None:
            outputs.hraw = topo
            if hmin is not None and hmax is not None: #Means you are in zoom AGRIF        
                GetMask.mask(GetMask,outputs,shpfile,hmin=hmin,sgl_connect=sgl_connect,prt_grd=prt_grd,ref_coef=coef)   
                topo=eval(''.join(("toolsf.",smooth.smooth,"(topo,hmin,hmax, \
                                    smooth.rfact,outputs.mask_rho)")))
            else:
                if prt_grd is not None: #compute approximatively the ref coeff
                    coef=int(np.nanmean(outputs.pm)/np.nanmean(prt_grd.pm))
                    print('ratio between prt and chld grid is approx:', coef)
           
                GetMask.mask(GetMask,outputs,shpfile,hmin=smooth.depthmin,sgl_connect=sgl_connect,prt_grd=prt_grd)
                topo=eval(''.join(("toolsf.",smooth.smooth,"(topo,smooth.depthmin,smooth.depthmax, \
                                    smooth.rfact,outputs.mask_rho)")))
            outputs.h=topo
            return outputs
        else:
            outputs.hraw = topo
            return outputs

    def match_topo(self,prt,outputs,WESN):
        WEST=True if 'West' in WESN[0] else False
        EAST=True if 'East' in WESN[0] else False
        SOUTH=True if 'South' in WESN[0] else False
        NORTH=True if 'North' in WESN[0] else False
        mask_tmp=np.copy(outputs.mask_rho)
        topo=toolsf.r2r_match_topo(WEST,EAST,SOUTH,NORTH,WESN[1],outputs.lon_rho.T,outputs.lat_rho.T,outputs.h.T,mask_tmp.T,\
                                  prt.lon_rho.T,prt.lat_rho.T,prt.h.T)
        outputs.h=topo.T
        return outputs

class GetMask():
    def outline(lon, lat):
        '''
        Return lon, lat of perimeter around the grid
        '''
        def func(var):
            return np.hstack([var[:, 0], var[-1, 1:-1],
                              var[::-1, -1], var[0, ::-1][1:]])
        return func(lon), func(lat)

    def process_mask(maskin):
 
        print('Processing mask to close narrow bay and narrow land (1 point wide)')
        maskout=np.copy(maskin)
        [M,L]=maskout.shape
        neibmask=0.*maskout
        neibmask[2:-2,2:-2]=maskout[1:-3,2:-2]+maskout[3:-1,2:-2]+\
                            maskout[2:-2,1:-3]+maskout[2:-2,3:-1]
        cpt=0
        
        while np.sum(((neibmask[2:-2,2:-2]>=3) & (maskout[2:-2,2:-2]==0)) |
                      ((neibmask[2:-2,2:-2]<=1) & (maskout[2:-2,2:-2]==1)))>0:
            cpt+=1
            maskout[(neibmask>=3) & (maskout==0)]=1
            maskout[(neibmask<=1) & (maskout==1)]=0
 
            neibmask[2:-2,2:-2]=maskout[1:-3,2:-2]+maskout[3:-1,2:-2]+\
                                maskout[2:-2,1:-3]+maskout[2:-2,3:-1]
            if cpt>15:
                print('Process mask did not converge inside the domain')
                break

        maskout[:,0:2]=maskin[:,0:2]
        maskout[:,-2:]=maskin[:,-2:]
        maskout[0:2,:]=maskin[0:2,:]
        maskout[-2:,:]=maskin[-2:,:]

        return maskout


#     def mask(self, outputs,gfile,sgl_connect=None):
    def mask(self, outputs,gfile,hmin=None,sgl_connect=None,prt_grd=None,ref_coef=None):
        # handle basic mask
        if hmin is None or hmin>0: # Case hmin >0
            llcrnrlon = outputs.lon_rho.ravel().min()
            urcrnrlon = outputs.lon_rho.ravel().max()
            llcrnrlat = outputs.lat_rho.ravel().min()
            urcrnrlat = outputs.lat_rho.ravel().max()
            
            geoshp=gp.read_file(gfile,bbox=(llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat))
            if urcrnrlon>180: # 0-360 format
                # select -180-urcrnrlon part that was not previously included
                geoshp_bis = gp.read_file(gfile,
                               bbox=(-180,llcrnrlat,urcrnrlon-360,urcrnrlat))
                # put it in 0-360 format
                geoshp_bis = gp.GeoDataFrame(geometry=geoshp_bis.translate(xoff=360))
                # Merge both side
                geoshp = gp.GeoDataFrame( pd.concat( [geoshp,geoshp_bis], ignore_index=True) )

            ny,nx = outputs.lon_rho.shape
            outputs.mask_rho = np.zeros([ny,nx])
            n2max = 50000 # max point by chunk

            if geoshp.shape[0]>3000:
                nchunk = int(np.max([np.sqrt((ny)*(nx)/n2max),1]))
            else :
              nchunk = 1

            if nchunk>1:
                print(f"Chunk format (y,x):({nchunk},{nchunk})")
 
            for i,j in product(list(range(nchunk)),list(range(nchunk))):
                if nchunk>1:
                    print(f"Doing chunk ({i+1},{j+1})")
                dx1=2; dx2=2; dy1=2; dy2=2 #chunk overlap
                if i==0: dx1=0
                if i==nchunk-1: dx2=0
                if j==0: dy1=0
                if j==nchunk-1: dy2=0

                nx1i = int(i*(nx)/nchunk-2*dx1)
                nx2i = int((i+1)*(nx)/nchunk+2*dx2)
                ny1i = int(j*(ny)/nchunk-2*dy1)
                ny2i = int((j+1)*(ny)/nchunk+2*dy2)
    
                llcrnrlon = np.nanmin(outputs.lon_rho[ny1i:ny2i,nx1i:nx2i])
                urcrnrlon = np.nanmax(outputs.lon_rho[ny1i:ny2i,nx1i:nx2i])
                llcrnrlat = np.nanmin(outputs.lat_rho[ny1i:ny2i,nx1i:nx2i])
                urcrnrlat = np.nanmax(outputs.lat_rho[ny1i:ny2i,nx1i:nx2i])

                gs=geoshp.clip((llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat))
                try: # Error if no polygon in the chunck
                    rmask = regionmask.mask_geopandas(
                              gs.geometry,outputs.lon_rho[ny1i:ny2i,nx1i:nx2i],outputs.lat_rho[ny1i:ny2i,nx1i:nx2i])
                    outputs.mask_rho[ny1i+dy1:ny2i-dy2,nx1i+dx1:nx2i-dx2]\
                        [np.isnan(rmask[dy1:ny2i-ny1i-dy2,dx1:nx2i-nx1i-dx2])] = 1
                except:# Put 1 everywhere in the chunck when no polygon
                    outputs.mask_rho[ny1i+dy1:ny2i-dy2,nx1i+dx1:nx2i-dx2] = 1
                    continue
 
        else: # case hmin<=0
            outputs.mask_rho=np.ones(outputs.hraw.shape)
            if hmin<np.nanmin(outputs.hraw.ravel()):
                hmin=np.ceil(np.nanmin(outputs.hraw.ravel()))
            outputs.mask_rho[outputs.hraw<=hmin]=0       

#        ##
        if prt_grd is not None: # Means we are in a zoom and we use prt grid mask at boundaries
            print('Matching Parent and Child mask close to boundary')
            spline = itp.NearestNDInterpolator((prt_grd.lon_rho.ravel(),prt_grd.lat_rho.ravel()),prt_grd.mask_rho.ravel())
            maskr_coarse = spline((outputs.lon_rho,outputs.lat_rho)) # prt mask on chd grid

            outputs.mask_rho[0:2,:]  = maskr_coarse[0:2,:]
            outputs.mask_rho[-2:,:] = maskr_coarse[-2:,:]
            outputs.mask_rho[:,0:2]  = maskr_coarse[:,0:2]
            outputs.mask_rho[:,-2:] = maskr_coarse[:,-2:]
            if ref_coef is not None:# check whether there is a coef
                [M,L]=outputs.mask_rho.shape
                [imat,jmat]=np.meshgrid(np.arange(L),np.arange(M))
                dist=0*outputs.mask_rho+np.inf
                for j in range(M):
                    if outputs.mask_rho[j,0]==1:
                        dist=np.minimum(dist,np.sqrt((imat)**2+(jmat-j)**2))
                    if outputs.mask_rho[j,-1]==1:
                        dist=np.minimum(dist,np.sqrt((imat-L+1)**2+(jmat-j)**2))
                for i in range(L):
                    if outputs.mask_rho[0,i]==1:
                        dist=np.minimum(dist,np.sqrt((imat-i)**2+(jmat)**2))
                    if outputs.mask_rho[-1,i]==1:
                        dist=np.minimum(dist,np.sqrt((imat-i)**2+(jmat-M+1)**2))

#                # Put the parent mask close to the boundaries
                nmsk=1+2*ref_coef
                outputs.mask_rho[dist<=nmsk]=maskr_coarse[dist<=nmsk]
                test=np.copy(outputs.mask_rho)
                outputs.mask_rho=self.process_mask(outputs.mask_rho)

        if sgl_connect is not None:
            if sgl_connect[0]:
                if outputs.mask_rho[sgl_connect[2],sgl_connect[1]]<0.5 :
                    print('ERROR: selected point i =', sgl_connect[1], 'j =', sgl_connect[2], 'is on land. Try another point.')
                    exit()
                outputs.mask_rho=toolsf.single_connect(sgl_connect[1],sgl_connect[2],outputs.mask_rho.T).T
        return outputs

class EasyGrid():
    """
    EasyGrid object. Implements both the easygrid computation, and
    the picture acquisition.
    """
    '''
    def easygrid(self, inputs, outputs):
    
        lon_rho = np.linspace(inputs.tra_lon-size_x/111.,
                              inputs.tra_lon+size_x/111.,inputs.nx)
        lat_rho = np.linspace(inputs.tra_lat-size_y/111.,
                              inputs.tra_lat+size_y/111.,inputs.ny)
                            
        lon_rho, lat_rho = np.meshgrid(lon_rho,lat_rho)
        
        outputs.lon_rho = lon_rho
        outputs.lat_rho = lat_rho
        
        return outputs
    '''
   
    def easygrid(self, inputs, outputs):
        """
        Easy grid makes a rectangular, orthogonal grid with minimal gridsize variation
        It uses a Mercator projection around the equator and then rotates the sphere around its three axes to position the grid wherever it is desired.
   
        Inputs:
          nx:      Number of grid point in the x direction
          ny:      Number of grid point in the y direction
          size_x:  Domain size in x-direction
          size_y:  Domain size in y-direction
          tra_lon: Desired longitude of grid center
          tra_lat: Desired latitude of grid center
          rot:     Rotation of grid direction (0: x direction is west-east)

        Example:  > lon_rho,lat_rho,pm,pn,ang,lon_psi,lat_psi = easy_grid(30,20,4e6,3e6,-180,59,0)

        Translated to Python from EGRID (c) (Matlab) by Jeroen Molemaker, UCLA, 2008
        """
        def tra_sphere(lon1, lat1, tra):
            """
            Translate sphere about its y-axis
            Part of easy grid
            (c) 2008, Jeroen Molemaker, UCLA
            """

            n, m = lon1.shape
            tra = np.deg2rad(tra) # translation in latitude direction

            # translate into x,y,z
            # conventions:  (lon,lat) = (0,0)  corresponds to (x,y,z) = ( 0,-r, 0)
            #               (lon,lat) = (0,90) corresponds to (x,y,z) = ( 0, 0, r)
            x1 = np.sin(lon1) * np.cos(lat1)
            y1 = np.cos(lon1) * np.cos(lat1)
            z1 = np.sin(lat1)

            """
            We will rotate these points around the small circle defined by
            the intersection of the sphere and the plane that
            is orthogonal to the line through (lon,lat) (90,0) and (-90,0)

            The rotation is in that plane around its intersection with
            aforementioned line.

            Since the plane is orthogonal to the x-axis (in my definition at least),
            rotations in the plane of the small circle maintain constant x and are around
            (x,y,z) = (x1,0,0)
            """

            rp1 = np.sqrt(y1 ** 2 + z1 ** 2)

            ap1 = 0.5 * np.pi * np.ones((n, m))
            tmpi = np.abs(y1) > 1.e-5 
            ap1[tmpi] = np.arctan(np.abs(z1[tmpi] /
                                         y1[tmpi]))
            tmpi = y1 < 0.
            ap1[tmpi] = np.pi - ap1[tmpi]
            tmpi = z1 < 0.
            ap1[tmpi] = -ap1[tmpi]

            ap2 = ap1 + tra
            x2  = x1.copy()
            y2  = rp1 * np.cos(ap2)
            z2  = rp1 * np.sin(ap2)

            # transformation from (x,y,z) to (lat,lon)
            lon2 = 0.5 * np.pi * np.ones((n, m))
            tmpi = np.abs(y2) > 1.e-5
            lon2[tmpi] = np.arctan(np.abs(x2[tmpi] /
                                          y2[tmpi]))
            tmpi = y2 < 0.
            lon2[tmpi] = np.pi - lon2[tmpi]
            tmpi = x2 < 0.
            lon2[tmpi] = -lon2[tmpi]

            pr2 = np.sqrt(x2 ** 2 + y2 ** 2)
            lat2 = 0.5 * np.pi * np.ones((n, m))
            tmpi = np.abs(pr2) > 1.e-5
            lat2[tmpi] = np.arctan(np.abs(z2[tmpi] /
                                         pr2[tmpi]))
            tmpi = z2 < 0.
            lat2[tmpi] = -lat2[tmpi]

            return lon2, lat2
        
        
        def rot_sphere(lon1, lat1, rot):
            """
            Rotate sphere around its y-axis
            Part of Easy Grid
            (c) 2008, Jeroen Molemaker, UCLA
            """
            n, m = lon1.shape
            rot = np.deg2rad(rot)
            eps=1e-20
            # translate into x,y,z
            # conventions:  (lon,lat) = (0,0)  corresponds to (x,y,z) = ( 0,-r, 0)
            #                   (lon,lat) = (0,90) corresponds to (x,y,z) = ( 0, 0, r)
            x1 = np.sin(lon1) * np.cos(lat1)
            y1 = np.cos(lon1) * np.cos(lat1)
            z1 = np.sin(lat1)
            """
            We will rotate these points around the small circle defined by
            the intersection of the sphere and the plane that
            is orthogonal to the line through (lon,lat) (0,0) and (180,0)

            The rotation is in that plane around its intersection with
            aforementioned line.

            Since the plane is orthogonal to the y-axis (in my definition at least),
            rotations in the plane of the small circle maintain constant y and are around
            (x,y,z) = (0,y1,0)
            """
            rp1 = np.sqrt(x1 ** 2 + z1 ** 2)

            ap1 = 0.5 * np.pi * np.ones((n, m))
            ap1[np.abs(x1) > eps] = np.arctan(np.abs(z1[np.abs(x1) > eps]/x1[np.abs(x1) > eps]))
            ap1[x1 < 0.] = np.pi - ap1[x1 < 0.]
            ap1[z1 < 0.] = -ap1[z1 < 0.]

            ap2 = ap1 + rot
            x2 = rp1 * np.cos(ap2)
            y2 = y1.copy()
            z2 = rp1 * np.sin(ap2)

            lon2 = 0.5 * np.pi * np.ones((n, m))
            lon2[np.abs(y2) > eps] = np.arctan(np.abs(x2[np.abs(y2) > eps]/ y2[np.abs(y2) > eps]))
            lon2[y2 < 0.] = np.pi - lon2[y2 < 0.]
            lon2[x2 < 0.] = -lon2[x2 < 0.]

            pr2 = np.hypot(x2, y2)
            lat2 = 0.5 * np.pi * np.ones((n, m))
            lat2[np.abs(pr2) > eps] = np.arctan(np.abs(z2[np.abs(pr2) > eps]/pr2[np.abs(pr2) > eps]))
            lat2[z2 < 0] = -lat2[z2 < 0]

            return lon2, lat2
        
        
        def gc_dist(lon1, lat1, lon2, lat2):
            '''
            Use Haversine formula to calculate distance
            between one point and another
            '''
            dlat = lat2 - lat1
            dlon = lon2 - lon1
            dang = 2. * np.arcsin(np.sqrt(np.power(np.sin(0.5 * dlat),2) + \
                        np.cos(lat2) * np.cos(lat1) * np.power(np.sin(0.5 * dlon),2)))
            return r_earth * dang # distance


        r_earth = 6371315. # Mean earth radius in metres (from scalars.h)

        size_x = inputs.size_x * 1000. # convert to meters
        size_y = inputs.size_y * 1000. # convert to meters

        # Mercator projection around the equator
        if size_y > size_x:
            length = np.float64(size_y)
            nl = np.float64(inputs.ny)
            width = np.float64(size_x)
            nw = np.float64(inputs.nx)
        else:
            length = np.float64(size_x)
            nl = np.float64(inputs.nx)
            width = np.float64(size_y)
            nw = np.float64(inputs.ny)
    
        dlon = length / r_earth
        lon1d = (dlon * np.arange(-0.5, nl + 1.5, 1.) / nl ) - ( 0.5 * dlon)

        mul = 1.
        dlat = width / r_earth

        for it in range(100):
            y1 = np.log(np.tan((0.25 * np.pi) - (0.25 * dlat)))
            y2 = np.log(np.tan((0.25 * np.pi) + (0.25 * dlat)))
            y = ((y2 - y1) * np.arange(-0.5, nw + 1.5, 1.) / nw ) + y1
            
            lat1d = np.arctan(np.sinh(y))
            dlat_cen = 0.5 * (lat1d[np.int32(np.round(0.5 * nw))] -  
                              lat1d[np.int32(np.round(0.5 * nw) - 2)])
            dlon_cen = dlon / nl
            mul = (dlat_cen / dlon_cen) * (length/width) * (nw / nl)
            dlat /= mul

        lon1de = (dlon * np.arange(-1., nl + 2., 1.) / nl ) - (0.5 * dlon)
        ye = ((y2-y1) * np.arange(-1., nw + 2., 1.) / nw ) + y1
        lat1de = np.arctan(np.sinh(ye))
        lat1de /= mul
        
        lon1, lat1 = np.meshgrid(lon1d, lat1d)
        lone, late = np.meshgrid(lon1de, lat1de)
        lonu = 0.5 * (lon1[:, :-1] + lon1[:, 1:])
        latu = 0.5 * (lat1[:, :-1] + lat1[:, 1:])
        lonv = 0.5 * (lon1[:-1] + lon1[1:])
        latv = 0.5 * (lat1[:-1] + lat1[1:])
        
        if size_y > size_x:
            lon1, lat1 = rot_sphere(lon1, lat1, 90.)
            lonu, latu = rot_sphere(lonu, latu, 90.)
            lonv, latv = rot_sphere(lonv, latv, 90.)
            lone, late = rot_sphere(lone, late, 90.)
            
            lon1 = lon1[::-1].T
            lat1 = lat1[::-1].T
            lone = lone[::-1].T
            late = late[::-1].T

            lonu_tmp = lonv[::-1].T
            latu_tmp = latv[::-1].T
            lonv = lonu[::-1].T
            latv = latu[::-1].T
            lonu = lonu_tmp
            latu = latu_tmp

        lon2, lat2 = rot_sphere(lon1, lat1, inputs.rot)
        lonu, latu = rot_sphere(lonu, latu, inputs.rot)
        lonv, latv = rot_sphere(lonv, latv, inputs.rot)
        lone, late = rot_sphere(lone, late, inputs.rot)

        lon3, lat3 = tra_sphere(lon2, lat2, inputs.tra_lat)
        lonu, latu = tra_sphere(lonu, latu, inputs.tra_lat)
        lonv, latv = tra_sphere(lonv, latv, inputs.tra_lat)
        lone, late = tra_sphere(lone, late, inputs.tra_lat)

        lon4 = lon3 + np.deg2rad(inputs.tra_lon)
        lonu += np.deg2rad(inputs.tra_lon)
        lonv += np.deg2rad(inputs.tra_lon)
        lone += np.deg2rad(inputs.tra_lon)
        
        pi2 = 2 * np.pi 
        lon4[lon4 < -np.pi] += pi2
        lonu[lonu < -np.pi] += pi2
        lonv[lonv < -np.pi] += pi2
        lone[lone < -np.pi] += pi2
        lat4 = lat3.copy()

        # Compute pm and pn
        pmu = gc_dist(lonu[:, :-1], latu[:, :-1],
                      lonu[:, 1:], latu[:, 1:])
        pm = np.zeros_like(lon4)
        pm[:, 1:-1] = pmu
        pm[:, 0] = pm[:, 1]
        pm[:, -1] = pm[:, -2]
        pm = 1. / pm
        
        pnv = gc_dist(lonv[:-1], latv[:-1],
                      lonv[1:], latv[1:])
        pn = np.zeros_like(lon4)
        pn[1:-1] = pnv
        pn[0] = pn[1]
        pn[-1] = pn[-2]
        pn = 1. / pn

        # Compute angles of local grid positive x-axis relative to east
        dellat = latu[:, 1:] - latu[:, :-1]
        dellon = lonu[:, 1:] - lonu[:, :-1]
        dellon[dellon > np.pi] -= pi2
        dellon[dellon < -np.pi] += pi2
        dellon = dellon * np.cos(0.5 * (latu[:, 1:] + latu[:, :-1]) )

        ang = np.zeros_like(lon4)
        ang_s = np.arctan(dellat / (dellon + 1.e-16))
        deli = np.logical_and(dellon < 0., dellat < 0.)
        ang_s[deli] -= np.pi
        
        deli = np.logical_and(dellon < 0., dellat >= 0.)
        ang_s[deli] += np.pi
        ang_s[ang_s > np.pi] -= np.pi
        ang_s[ang_s < -np.pi] += np.pi

        ang[:, 1:-1] = ang_s
        ang[:, 0] = ang[:, 1]
        ang[:, -1] = ang[:, -2]

        outputs.lon_rho = np.rad2deg(lon4)
        outputs.lat_rho = np.rad2deg(lat4)
        outputs.lon_u = np.rad2deg(lonu)
        outputs.lat_u = np.rad2deg(latu)
        outputs.lon_v = np.rad2deg(lonv)
        outputs.lat_v = np.rad2deg(latv)
        outputs.pm = pm
        outputs.pn = pn
        outputs.angle = ang
        outputs.lon_psi = np.rad2deg(lone)
        outputs.lat_psi = np.rad2deg(late)

        return outputs

    def AGRIFgrid(self,prt_grd,inputs,outputs):
        def gc_dist(lon1, lat1, lon2, lat2):
            r_earth=6371315. # Mean earth radius in metres (from scalars.h
            '''
            Use Haversine formula to calculate distance
            between one point and another
            '''
            dlat = lat2 - lat1
            dlon = lon2 - lon1
            dang = 2. * np.arcsin(np.sqrt(np.power(np.sin(0.5 * dlat),2) + \
                        np.cos(lat2) * np.cos(lat1) * np.power(np.sin(0.5 * dlon),2)))
            return r_earth * dang # distance

        [Mp,Lp]=prt_grd.lon_rho.shape
        
        igrdp=np.arange(0,Lp-1);jgrdp=np.arange(0,Mp-1)
        igrdr=np.arange(0,Lp);jgrdr=np.arange(0,Mp)
        igrdu=np.arange(0,Lp-1);jgrdu=np.arange(0,Mp)
        igrdv=np.arange(0,Lp);jgrdv=np.arange(0,Mp-1)
       
        [iprtgrd_p,jprtgrd_p]=np.meshgrid(igrdp,jgrdp)
        [iprtgrd_r,jprtgrd_r]=np.meshgrid(igrdr,jgrdr)
        [iprtgrd_u,jprtgrd_u]=np.meshgrid(igrdp,jgrdr)
        [iprtgrd_v,jprtgrd_v]=np.meshgrid(igrdr,jgrdp)

        bbound_east=1
        bbound_west=1
        bbound_south=1
        bbound_north=1

        while (bbound_east | bbound_west | bbound_south | bbound_north):
            
            ipch=np.arange(inputs.imin,inputs.imax+0.5/inputs.coef,1/inputs.coef)
            jpch=np.arange(inputs.jmin,inputs.jmax+0.5/inputs.coef,1/inputs.coef)
            ################
            irch=np.arange(inputs.imin+0.5-0.5/inputs.coef,inputs.imax+0.5+0.75/inputs.coef,1/inputs.coef)
            jrch=np.arange(inputs.jmin+0.5-0.5/inputs.coef,inputs.jmax+0.5+0.75/inputs.coef,1/inputs.coef)

            [ichildgrd_p,jchildgrd_p]=np.meshgrid(ipch,jpch)
            [ichildgrd_r,jchildgrd_r]=np.meshgrid(irch,jrch)
            [ichildgrd_u,jchildgrd_u]=np.meshgrid(ipch,jrch)
            [ichildgrd_v,jchildgrd_v]=np.meshgrid(irch,jpch)
         
            spline_lonp    = itp.RectBivariateSpline(jgrdp,igrdp,prt_grd.lon_psi)
            spline_latp    = itp.RectBivariateSpline(jgrdp,igrdp,prt_grd.lat_psi)
            outputs.lon_psi = spline_lonp(jchildgrd_p,ichildgrd_p,grid=False)
            outputs.lat_psi = spline_latp(jchildgrd_p,ichildgrd_p,grid=False)
            ###############
            spline_lonr    = itp.RectBivariateSpline(jgrdr,igrdr,prt_grd.lon_rho)
            spline_latr    = itp.RectBivariateSpline(jgrdr,igrdr,prt_grd.lat_rho)
            outputs.lon_rho = spline_lonr(jchildgrd_r,ichildgrd_r,grid=False)
            outputs.lat_rho = spline_latr(jchildgrd_r,ichildgrd_r,grid=False)
            ###############
            spline_lonu    = itp.RectBivariateSpline(jgrdu,igrdu,prt_grd.lon_u)
            spline_latu    = itp.RectBivariateSpline(jgrdu,igrdu,prt_grd.lat_u)
            outputs.lon_u = spline_lonr(jchildgrd_u,ichildgrd_u,grid=False)
            outputs.lat_u = spline_latr(jchildgrd_u,ichildgrd_u,grid=False)
            ###############
            spline_lonv    = itp.RectBivariateSpline(jgrdv,igrdv,prt_grd.lon_v)
            spline_latv    = itp.RectBivariateSpline(jgrdv,igrdv,prt_grd.lat_v)
            outputs.lon_v = spline_lonr(jchildgrd_v,ichildgrd_v,grid=False)
            outputs.lat_v = spline_latr(jchildgrd_v,ichildgrd_v,grid=False)
            #################

            maskr_coarse = griddata((jprtgrd_r.ravel(),iprtgrd_r.ravel()),prt_grd.mask_rho.ravel(),\
                                    (jchildgrd_r.ravel(),ichildgrd_r.ravel()),method='nearest')

            maskr_coarse = np.reshape(maskr_coarse,[jchildgrd_r.shape[0],jchildgrd_r.shape[1]])
         
            eastchk = abs(maskr_coarse[:,-2]-maskr_coarse[:,-1]);
            westchk = abs(maskr_coarse[:,0]-maskr_coarse[:,1]);
            southchk = abs(maskr_coarse[0,:]-maskr_coarse[1,:]);
            northchk = abs(maskr_coarse[-2,:]-maskr_coarse[-1,:]);

            if sum(eastchk)!=0:
                inputs.imax=inputs.imax+1
                bbound_east=1
                print('==> East limits displacement +1')
            else:
                bbound_east=0
     
            if sum(westchk)!=0:
                inputs.imin=inputs.imin-1
                bbound_west=1
                print('==> West limits displacement -1')
            else:
                bbound_west=0

            if sum(southchk)!=0:
                inputs.jmin=inputs.jmin-1
                bbound_south=1
                print('==> South limits displacement -1')
            else:
                bbound_south=0

            if sum(northchk)!=0:
                inputs.jmax=inputs.jmax+1
                bbound_north=1
                print('==> North limits displacement +1')
            else:
                bbound_north=0



        # Compute pm and pn
        pmu = gc_dist(np.deg2rad(outputs.lon_u[:, :-1]), np.deg2rad(outputs.lat_u[:, :-1]),
                      np.deg2rad(outputs.lon_u[:, 1:]), np.deg2rad(outputs.lat_u[:, 1:]))
        pm = np.zeros_like(outputs.lon_rho)
        pm[:, 1:-1] = pmu
        pm[:, 0] = pm[:, 1]
        pm[:, -1] = pm[:, -2]
        pm = 1. / pm

        pnv = gc_dist(np.deg2rad(outputs.lon_v[:-1]),np.deg2rad(outputs.lat_v[:-1]),
                      np.deg2rad(outputs.lon_v[1:]), np.deg2rad(outputs.lat_v[1:]))
        pn = np.zeros_like(outputs.lon_rho)
        pn[1:-1] = pnv
        pn[0] = pn[1]
        pn[-1] = pn[-2]
        pn = 1. / pn

        # Compute angles of local grid positive x-axis relative to east
        pi2 = 2 * np.pi
        dellat = np.deg2rad(outputs.lat_u[:, 1:]) - np.deg2rad(outputs.lat_u[:, :-1])
        dellon = np.deg2rad(outputs.lon_u[:, 1:]) - np.deg2rad(outputs.lon_u[:, :-1])
        dellon[dellon > np.pi] -= pi2
        dellon[dellon < -np.pi] += pi2
        dellon = dellon * np.cos(0.5 * (np.deg2rad(outputs.lat_u[:, 1:]) +np.deg2rad(outputs.lat_u[:, :-1])) )

        ang = np.zeros_like(outputs.lon_rho)
        ang_s = np.arctan(dellat / (dellon + 1.e-16))
        deli = np.logical_and(dellon < 0., dellat < 0.)
        ang_s[deli] -= np.pi

        deli = np.logical_and(dellon < 0., dellat >= 0.)
        ang_s[deli] += np.pi
        ang_s[ang_s > np.pi] -= np.pi
        ang_s[ang_s < -np.pi] += np.pi

        ang[:, 1:-1] = ang_s
        ang[:, 0] = ang[:, 1]
        ang[:, -1] = ang[:, -2]

        outputs.pm = pm
        outputs.pn = pn
        outputs.angle = ang

        return outputs

def indx_bound(x, x0):
    """
    Conversion of fortran tools indx_bound
    """
    n=x.shape[0]
    if x0 < x[0] :
        i=-1                     # if x0 is outside the full range
    elif x0 > x[-1] :            # of x(1) ... x(n), then return
        i=n                      # i=0 or i=n.
    else:
        i=int( ( x[-1]-x0 +(n-1)*(x0-x[0]) )/(x[-1]-x[0]) )
        if x[i+1]<x0 :
            while x[i+1] <x0 :  # This algorithm computes "i" as
                i=i+1           # linear interpolation between x(1)
                                # and x(n) which should yield the
        elif x[i] > x0 :        # correct value for "i" right a way
            while x[i] > x0 :   # because array elements x(i) are
                i=i-1           # equidistantly spaced.  The while
                                # loops are here merely to address
                                # possible roundoff errors.

        if x[i+1]-x0 < 0 or x0-x[i] < 0 :
            print('### ERROR: indx_bound :: ',x[i], x0, x[i+1], x0-x[i], x[i+1]-x0)
            exit()
    indx_bound=i
    return indx_bound





##################################

#############################
### Class for normal mode ###
class inputs():
    '''
    Inputs to locate grid
    '''
    def __init__(self,tra_lon,tra_lat,size_x,size_y,nx,ny,rot):

        self.tra_lon = tra_lon
        self.tra_lat = tra_lat
        self.size_x  = size_x
        self.size_y  = size_y
        self.rot     = rot
        self.nx      = nx
        self.ny      = ny

class inputs_smth():
    '''
    Inputs for smoothing
    '''
    def __init__(self,depthmin,depthmax,smthr,rfact,smooth):
        self.depthmin  = depthmin
        self.depthmax  = depthmax
        self.smthr     = smthr
        self.rfact     = rfact
        self.smooth    = smooth
#######################
### Read parent grid ##

class topo_prt():

    def __init__(self,prt_file):
       nc=netcdf.Dataset(prt_file,'r')
       self.lon_rho  = nc.variables['lon_rho'][:]
       self.lat_rho  = nc.variables['lat_rho'][:]
       self.lon_psi  = nc.variables['lon_psi'][:]
       self.lat_psi  = nc.variables['lat_psi'][:]
       self.lon_u    = nc.variables['lon_u'][:]
       self.lat_u    = nc.variables['lat_u'][:]
       self.lon_v    = nc.variables['lon_v'][:]
       self.lat_v    = nc.variables['lat_v'][:]
       self.mask_rho = nc.variables['mask_rho'][:]
       self.h        = nc.variables['h'][:]
       self.pm       = nc.variables['pm'][:]
       self.pn       = nc.variables['pn'][:]
       nc.close()
