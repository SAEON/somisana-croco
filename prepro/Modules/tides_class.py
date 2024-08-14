import numpy as np
import netCDF4 as netcdf
import cftime
import xarray as xr
import tides_reader as dico
import os
import sys
import glob

class getdata():
    def __init__(self,inputdata,inputfile,crocogrd,file_format,tideslist,currentu=None,currentv=None):
 
        self.var=dico.lookvar(inputdata) # Dictionary to find the names of the input variables
        if file_format =='Amp_Phase':
            part1_suf='a'
            part2_suf='p'
        elif file_format =='Re_Im':
            part1_suf='r'
            part2_suf='i'
        else:
            sys.exit('input_type must be Amp_Phase or Re_Im depending on what is in your input file(s)')

        xr_list_ssh=[]
        for inpt in inputfile:
            xr_list_ssh+=[xr.open_dataset(inpt)]
        dataxr_ssh=xr.concat(xr_list_ssh,dim='ntides',data_vars='different')
        
        self.ncglo={'ssh_part1':eval(''.join(("dataxr_ssh."+self.var['H_'+part1_suf]))),\
                    'ssh_part2':eval(''.join(("dataxr_ssh."+self.var['H_'+part2_suf])))
                   }
        
        # Handling current variables
        if currentu is not None: 
            xr_list_u=[]
            xr_list_v=[]
            for inpt in currentu:
                xr_list_u+=[xr.open_dataset(inpt)]
            dataxr_u=xr.concat(xr_list_u,dim='ntides',data_vars='different')
            for inpt in currentv:
                xr_list_v+=[xr.open_dataset(inpt)]
            dataxr_v=xr.concat(xr_list_v,dim='ntides',data_vars='different')
     
            self.ncglo['u_part1']=eval(''.join(("dataxr_u."+self.var['U_'+part1_suf])))
            self.ncglo['u_part2']=eval(''.join(("dataxr_u."+self.var['U_'+part2_suf])))
            self.ncglo['v_part1']=eval(''.join(("dataxr_v."+self.var['V_'+part1_suf])))
            self.ncglo['v_part2']=eval(''.join(("dataxr_v."+self.var['V_'+part2_suf])))
            
            if 'tpxo' in inputdata and \
               'transport' in self.ncglo['u_part1'].long_name:
            # set transport to velocity           
                print('Looking for grid file')
                grd_file=glob.glob(
                             currentu[0].replace(currentu[0].split('/')[-1],
                             'grid*'))
                # check transport units to change to m2/s
                if ('cm' in self.ncglo['u_part1'].attrs['units'] or 
                    'centimeter' in self.ncglo['u_part1'].attrs['units']):
                   coef = 100
                elif ('mm' in self.ncglo['u_part1'].attrs['units'] or
                    'millimeter' in self.ncglo['u_part1'].attrs['units']):
                   coef = 1000
                else:
                   coef = 1

                try:
#                   if len(grd_file)>0 and os.path.exists(grd_file[0]):
                    print('Found '+grd_file[0]+' for grid file')
                    xar=xr.open_dataset(grd_file[0])
                    hu=eval(''.join((f"xar.{self.var['topou']}")))
                    hv=eval(''.join((f"xar.{self.var['topov']}")))
                except:
                    try:
                        hu=eval(''.join((f"dataxr_u.{self.var['topou']}")))
                        hv=eval(''.join((f"dataxr_v.{self.var['topov']}")))
                    except:
                        try:
                            hu=eval(''.join((f"dataxr_u.{self.var['topor']}")))
                            hv=eval(''.join((f"dataxr_v.{self.var['topor']}")))
                        except:
                            try:
                                hu=eval(''.join((f"dataxr_ssh.{self.var['topor']}")))
                                hv=eval(''.join((f"dataxr_ssh.{self.var['topor']}")))
                            except:
                                sys.exit('Need either a grid file or topo in tpxo file to correct transport.') 
                self.ncglo['u_part1']=self.ncglo['u_part1']/coef**2/hu.values
                self.ncglo['u_part1'].attrs['units']='meter/sec'
                self.ncglo['u_part2']=self.ncglo['u_part2']/coef**2/hu.values
                self.ncglo['u_part2'].attrs['units']='meter/sec'
                self.ncglo['v_part1']=self.ncglo['v_part1']/coef**2/hv.values
                self.ncglo['v_part1'].attrs['units']='meter/sec'
                self.ncglo['v_part2']=self.ncglo['v_part2']/coef**2/hv.values
                self.ncglo['v_part2'].attrs['units']='meter/sec'                    
        # If all waves are in one file, reads period along record dimension 
        if len(inputfile)==1 and self.ncglo['ssh_part1'].ndim ==3:
            rcd_dim=self.ncglo['ssh_part1'].dims[0]
            self.ntides = np.array((eval(''.join(("dataxr_ssh."+rcd_dim+".data")))))
       
        # If len(inputfile)==1 and wave dimension is 2-d, need to add a dummy dimension   
        for key in self.ncglo:
            if self.ncglo[key].ndim == 2:
                self.ncglo[key]=self.ncglo[key].expand_dims(dim={'ntides':1},axis=0)       
        # Take care of unit
        for key in self.ncglo:
            try:         
                uni = self.ncglo[key].units
            except:
                print('No units in input file.')
                if 'part2' in key:
                    if file_format =='Amp_Phase':
                        print('looking at Phase. Assume it is in degree')
                        uni='degrees'
                    else:
                        if 'u' in self.ncglo[key] or 'v' in self.ncglo[key]: 
                            print('looking at velocity Imaginary part. Assume it is in meter.second-1')
                            uni='meter.second-1'
                        else:
                            print('looking at Imaginary part. Assume it is in meter')
                            uni='meter'
                elif 'part1' in key:
                    if file_format =='Amp_Phase':
                        if 'u' in self.ncglo[key] or 'v' in self.ncglo[key]:
                            print('looking at a velocity amplitude. Assume it is in meter.second-1')
                            uni='meter.second-1'
                        elif 'ssh' in self.ncglo[key]:
                            print('looking at an Elevation. Assume it is in meter')
                            uni='meter'
                    else:
                        if 'u' in self.ncglo[key] or 'v' in self.ncglo[key]:
                            print('looking at velocity Real part. Assume it is in meter.second-1')
                            uni='meter.second-1'
                        else:
                            print('looking at Real part. Assume it is in meter')
                            uni='meter'

            if uni.lower() in ['millimeters','millimeter','mm']:
                print('converting %s to meter' % uni)
                self.ncglo[key]=self.ncglo[key]/1000
            elif uni.lower() in ['centimeters','centimeter','cm']:
                print('converting %s to meter' % uni)
                self.ncglo[key]=self.ncglo[key]/100
            elif uni.lower() in ['degrees','degree']:
                print('converting %s to radian' % uni)
                self.ncglo[key]=self.ncglo[key]*np.pi/180
            elif uni.lower() in ['centimeter.second-1', 'cm.s-1','centimeters/second','cm/s','centimeter/sec','centimeters/sec']:
                print('converting %s to meter.second-1' % uni)
                self.ncglo[key]=self.ncglo[key]/100
            elif uni.lower() in ['millimeter.second-1', 'mm.s-1','millimeters/second','mm/s','milliimeter/sec','millimeters/sec']:
                print('converting %s to meter.second-1' % uni)
                self.ncglo[key]=self.ncglo[key]/1000
         

        [self.lonT ,self.latT ,self.idmin ,self.idmax ,self.jdmin ,self.jdmax ,self.period ]  = self.handle_periodicity(crocogrd,dataxr_ssh,'r')

        if currentu is not None:
            [self.lonU ,self.latU ,self.idminU ,self.idmaxU ,self.jdminU ,self.jdmaxU ,self.periodU ]  = self.handle_periodicity(crocogrd,dataxr_u,'u')
            [self.lonV ,self.latV ,self.idminV ,self.idmaxV ,self.jdminV ,self.jdmaxV ,self.periodV ]  = self.handle_periodicity(crocogrd,dataxr_v,'v')
       
        # Checking that var in format [ntides,lat,lon] 
        for key in self.ncglo:
            if 'u' in self.ncglo[key]:
                nx=eval(''.join(("dataxr_u."+self.var['lonu']))).shape[0]
                ny=eval(''.join(("dataxr_u."+self.var['latu']))).shape[-1]
                if (self.ncglo[key][:].shape[1] != ny and self.ncglo[key][:].shape[2] != nx ) and (self.ncglo[key][:].shape[1] == nx and self.ncglo[key][:].shape[2] == ny):
                    print('%s in format [ntides,Lon,Lat], switching Lon/Lat axes' % key)
                    dim=self.ncglo[key].dims
                    self.ncglo[key]=self.ncglo[key].transpose(dim[0],dim[2],dim[1])
            elif 'v' in self.ncglo[key]:
                nx=eval(''.join(("dataxr_v."+self.var['lonv']))).shape[0]
                ny=eval(''.join(("dataxr_v."+self.var['latv']))).shape[-1]
                if (self.ncglo[key][:].shape[1] != ny and self.ncglo[key][:].shape[2] != nx ) and (self.ncglo[key][:].shape[1] == nx and self.ncglo[key][:].shape[2] == ny):
                    print('%s in format [ntides,Lon,Lat], switching Lon/Lat axes' % key)
                    dim=self.ncglo[key].dims
                    self.ncglo[key]=self.ncglo[key].transpose(dim[0],dim[2],dim[1])
            else:
                nx=eval(''.join(("dataxr_ssh."+self.var['lonr']))).shape[0]
                ny=eval(''.join(("dataxr_ssh."+self.var['latr']))).shape[-1]
                if (self.ncglo[key][:].shape[1] != ny and self.ncglo[key][:].shape[2] != nx ) and (self.ncglo[key][:].shape[1] == nx and self.ncglo[key][:].shape[2] == ny):
                    print('%s in format [ntides,Lon,Lat], switching Lon/Lat axes' % key)
                    dim=self.ncglo[key].dims
                    self.ncglo[key]=self.ncglo[key].transpose(dim[0],dim[2],dim[1])       
        
    def ap2ep(self,Au, PHIu, Av, PHIv):
        '''
        This function compute Semi Major axis along with eccentricity and inclination and phase angle for tide ellipse
        Input:
          Au       Amplitude U tide
          PHIu     Phase u tide
          Av       Amplitude V tide
          PHIv     Phase V tide

        Output:
          SEMA     Semi  Major Axis, or maximum speed
          ECC      Eccentricity
          INC      Ellipse inclination
          PHA      Phase angle
        '''
        # Assume the input phase lags are in degrees and convert them in radians.
        PHIu = PHIu/180*np.pi
        PHIv = PHIv/180*np.pi

        # Make complex amplitudes for u and v
        u = Au*np.exp(-1j*PHIu)
        v = Av*np.exp(-1j*PHIv)

        # Calculate complex radius of anticlockwise and clockwise circles:
        wp = (u+1j*v)/2      # for anticlockwise circles
        wm = np.conj(u-1j*v)/2  # for clockwise circles

        # and their amplitudes and angles
        Wp = np.abs(wp)
        Wm = np.abs(wm)
        THETAp = np.angle(wp)
        THETAm = np.angle(wm)

        # calculate e-parameters (ellipse parameters)
        SEMA = Wp+Wm              # Semi  Major Axis, or maximum speed
        SEMI = Wp-Wm              # Semin Minor Axis, or minimum speed
        ECC = SEMI/SEMA          # Eccentricity

        PHA = (THETAm-THETAp)/2   # Phase angle, the time (in angle) when 
                                 # the velocity reaches the maximum
        INC = (THETAm+THETAp)/2   # Inclination, the angle between the 
                                 # semi major axis and x-axis (or u-axis).

        # convert to degrees for output
        PHA = PHA/np.pi*180
        INC = INC/np.pi*180
        THETAp = THETAp/np.pi*180
        THETAm = THETAm/np.pi*180

        # flip THETAp and THETAm, PHA, and INC in the range of 
        # [-pi, 0) to [pi, 2*pi), which at least is my convention.
        id = THETAp < 0;   THETAp[id] = THETAp[id]+360
        id = THETAm < 0;   THETAm[id] = THETAm[id]+360
        id = PHA < 0;      PHA[id] = PHA[id]+360
        id = INC < 0;      INC[id] = INC[id]+360

        return SEMA,  ECC, INC, PHA

    #####################################
    def indx_bound(self,x, x0):
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
                sys.exit()
        indx_bound=i
        return indx_bound



    def handle_periodicity(self,crocogrd,data,grid):
        '''
        handle_periodicity checks whether domain is inside the croco file.
        If so, check if there is a need for create a periodicity between
        the last and first longitude points ( for global data).
        It is returning lon/lat/topo adapted to the desired domain
        geolim = [lonmin,lonmax,latmin,latmax]
 
        input: inputfile : data file
               crocogrd  : croco grid lon/lat
               grid      : which grid (rho: 'r', u: 'u', v: 'v')
 
        output: lon/lat of the grid
                imin/imax index min/max xaxis
                jmin/jmax index min/max yaxis
        '''

        print('Reading coordinate file for grid: '+grid )
        print('-----------------------------------')

        lon=eval(''.join(("data."+self.var['lon'+grid])))
        lat=eval(''.join(("data."+self.var['lat'+grid])))

        if len(lon.shape)==2: # Some datasets are a bit different
            if lon[0,1]-lon[0,0]==0: # axes are swicthed
                lon = lon[:,0]
                lat = lat[0,:]
            else:
                lon = lon[0,:]
                lat = lat[:,0]
                
        for i in range(1,lon.shape[0]): # Fix discontinuity
            if lon[i]<lon[i-1]:        # between 180/-180 in the
                lon[i]=lon[i]+360      # middle
        ####
        geolim=[crocogrd.lonmin(),crocogrd.lonmax(),crocogrd.latmin(),crocogrd.latmax()]

        jmin=self.indx_bound(lat.data, geolim[2])
        jmax=self.indx_bound(lat.data, geolim[-1])
        
        if -1<jmin and jmin<lat.shape[0] and \
           -1<jmax and jmax<lat.shape[0] :
            if jmin > 1 :
                jmin=jmin-1
            jmax=jmax+2
        else:
            print('North-south extents of the dataset ',lat.data[0],lat.data[-1],' are not sufficient to cover the entire model grid.')
            sys.exit()
        ####
        imin=self.indx_bound(lon, geolim[0])
        imax=self.indx_bound(lon, geolim[1])

        if -1<imin and imin<lon.shape[0] and \
           -1<imax and imax<lon.shape[0] :
            if imin > 0:
                imin=imin-1
            imax=imax+1
            shft_west=0 ; shft_east=0 ; period=0
            print('Single region dataset imin/imax=',imin,imax)
        else:
        ######
            ptest=lon[-1]-lon[0]-360
            dx=(lon[-1]-lon[0])/(lon.shape[0]-1)
            epsil=0.01*abs(dx)
            if abs(ptest) < epsil :
                period=lon.shape[0]-1
            elif abs(ptest+dx) < epsil :
                period=lon.shape[0]
            else:
                period=0

            if period>0:
                print('Identified periodicity domain in data of ', period,' points out of', lon.shape[0])
            else :
                print('ERROR: The data does not cover the entire grid. Change your grid definition')
                sys.exit()
        ##
            shft_west=0
            if imin==-1 :
                shft_west=-1
                imin=self.indx_bound(lon, geolim[0]+360)
                if imin == lon.shape[0]: imin = lon.shape[0]-1
            elif imin==lon.shape[0] :
                shft_west=+1
                imin=self.indx_bound(lon, geolim[0]-360)
                if imin == -1: imin = lon.shape[0]-1
        ##
            shft_east=0
            if imax == -1:
                shft_east=-1
                imax=self.indx_bound(lon, geolim[1]+360)
                if imax == lon.shape[0]: imax = 0
            elif imax == lon.shape[0]:
                shft_east=+1
                imax=self.indx_bound(lon, geolim[1]-360)
                if imax == -1: imax = 0
    
            if -1<imin and imin<lon.shape[0] and \
               -1<imax and imax<lon.shape[0] :
                if imin>0:
                    imin=imin-1
                imax=imax+1
            else:
                print('ERROR: Data longitude covers 360 degrees, but still cannot find  starting and ending indices.')
                sys.exit()

        print('Bounding indices of the relevant part to be extracted from the entire dataset:\n', 
              'imin,imax =', imin,imax,'out of', lon.shape[0],'jmin,jmax =',jmin,jmax, 'out of',lat.shape[0])
        ny_lat=jmax-jmin+1
        start2=jmin ; end2=start2+ny_lat; count2=ny_lat
        lat_tmp=np.zeros([ny_lat])
        for j in range(0,ny_lat):
            lat_tmp[j]=lat[j+jmin]
        #####
        if imin < imax :
            nx_lon=imax-imin+1
            start1=imin ; end1=start1+nx_lon ; count1=nx_lon

            ishft=imin
            lon_tmp=np.zeros([nx_lon])
            if shft_west>0 and shft_east>0:
                for i in range(0,nx_lon):
                    lon_tmp[i]=lon[i+ishft] +360
            elif shft_west<0 and shft_east<0:
                for i in range(0,nx_lon):
                     lon_tmp[i]=lon[i+ishft]-360
            elif shft_west== 0 and shft_east==0:
                for i in range(0,nx_lon) :
                    lon_tmp[i]=lon[i+ishft]
            else:
                print('Error in shifting algoritm')
                sys.exit()
            (lon,lat)=np.meshgrid(lon_tmp,lat_tmp)
        ###
        elif imin>imax:
            print('Reading topography in two separate parts adjacent through 360-degree periodicity\n First...' )
            nx_lon=imax+period-imin+1
            xtmp  = np.zeros([nx_lon])
            start1=0 ; end1=start1+nx_lon; count1=imax
            ishft=nx_lon-count1-1
            if shft_east>0:
                for i in range(0,count1):
                    xtmp[i+ishft]=lon[i] +360
            elif shft_east<0:
                for i in range(0,count1):
                    xtmp[i+ishft]=lon[i] -360
            else:
                for i in range(0,count1):
                    xtmp[i+ishft]=lon[i]

            print('Second...')
            start1=imin ; count1=period-imin; end1=start1+count1
            ishft=imin
            if shft_west>0:
                for i in range(0,count1):
                    xtmp[i]=lon[i+ishft] +360
            elif shft_west<0 :
                for i in range(0,count1):
                    xtmp[i]=lon[i+ishft] -360
            else:
                for i in range(0,count1):
                    xtmp[i]=lon[i+ishft]
            lon_tmp=np.zeros([xtmp.shape[0]])
            for i in range(0,nx_lon):
                lon_tmp[i]=xtmp[i]

            del lon,lat
            (lon,lat)=np.meshgrid(lon_tmp,lat_tmp)
    
        return lon,lat,imin,imax,jmin,jmax,period

    #############################   
    def var_periodicity(self,vname,l,k,bdy=""):
        '''
        handle periodicity for tracers. Limits (imin,imax,jmin,jmax) are fixed before
        by runing 
        '''
        if vname in ['u','u_part1','u_part2']:
            imin=eval(''.join(("self.idminU"+bdy))) ; imax=eval(''.join(("self.idmaxU"+bdy)))
            jmin=eval(''.join(("self.jdminU"+bdy))) ; jmax=eval(''.join(("self.jdmaxU"+bdy)))
            period=eval(''.join(("self.periodU"+bdy))); grdid='u'
        elif vname in ['v','v_part1','v_part2']:
            imin=eval(''.join(("self.idminV"+bdy))) ; imax=eval(''.join(("self.idmaxV"+bdy)))
            jmin=eval(''.join(("self.jdminV"+bdy))) ; jmax=eval(''.join(("self.jdmaxV"+bdy)))
            period=eval(''.join(("self.periodV"+bdy)));grdid='v'
        else:
            imin=eval(''.join(("self.idmin"+bdy))) ; imax=eval(''.join(("self.idmax"+bdy)))
            jmin=eval(''.join(("self.jdmin"+bdy))) ; jmax=eval(''.join(("self.jdmax"+bdy)))
            period=eval(''.join(("self.period"+bdy)));grdid='r'
        
        try:
            mintime=min(l);maxtime=max(l)+1
        except:
            mintime=l;maxtime=l+1

        ny_lat=jmax-jmin+1
        start2=jmin ; end2=start2+ny_lat; count2=ny_lat

        if imin < imax :
            nx_lon=imax-imin+1
            start1=imin ; end1=start1+nx_lon ; count1=nx_lon
            if k==-1:
                field=np.array(np.squeeze(self.ncglo[vname][mintime:maxtime,start2:end2,start1:end1]))
            else:
                field=np.array(np.squeeze(self.ncglo[vname][mintime:maxtime,k,start2:end2,start1:end1]))

        elif imin>imax:    
            nx_lon=imax+period-imin+1
            try:
                lent=l.shape[0]
                ftmp = np.zeros([lent,ny_lat,nx_lon])
            except:
                ftmp = np.zeros([ny_lat,nx_lon])
            # First
            start1=0 ; end1=start1+nx_lon; count1=imax 
            if k==-1:
                field=np.array(np.squeeze(self.ncglo[vname][mintime:maxtime,start2:end2,start1:end1]))
            else:
                field=np.array(np.squeeze(self.ncglo[vname][mintime:maxtime,k,start2:end2,start1:end1]))

            if len(ftmp.shape)<3:
                for j in range(0,count2):
                    for i in range(0,count1):
                        ftmp[j,nx_lon-imax+i-1]=field[j,i]
            else:
                for j in range(0,count2):
                    for i in range(0,count1):
                        ftmp[:,j,nx_lon-imax+i-1]=field[:,j,i]

            del field

            # Second
            start1=imin ; count1=period-imin; end1=start1+count1
            if k==-1:
                field=np.array(np.squeeze(self.ncglo[vname][mintime:maxtime,start2:end2,start1:end1]))
            else:
                field=np.array(np.squeeze(self.ncglo[vname][mintime:maxtime,k,start2:end2,start1:end1]))

            if len(ftmp.shape)<3:
                for j in range(0,count2):
                    for i in range(0,count1):
                        ftmp[j,i]=field[j,i]
            else:
                for j in range(0,count2):
                    for i in range(0,count1):
                        ftmp[:,j,i]=field[:,j,i]

            del field

            field=np.copy(ftmp)

        return field

    def egbert_correction(self,input_wav,date_ini):
        '''
        Correct phases and amplitudes for real time runs
        Use parts of pos-processing code from Egbert's & Erofeeva's (OSU) 
        TPXO model. Their routines have been adapted from code by Richard Ray 
        (@?) and David Cartwright.

        Inputs:
          input_wav        Tidal wave on which to apply correction
          date_ini         Date of initialiation
        '''
        class phase_mkb():
            '''
            Astronomical arguments, obtained with Richard Ray's
            "arguments" and "astrol", for Jan 1, 1992, 00:00 Greenwich time
            Corrected July 12, 2000
            '''

            val=np.array([1.731557546,0.000000000,0.173003674,1.558553872,\
                6.050721243,6.110181633,3.487600001,5.877717569,\
                4.086699633,3.463115091,5.427136701,0.553986502,\
                0.052841931,2.137025284,2.436575100,1.929046130,\
                5.254133027,1.756042456,1.964021610,3.487600001,\
                3.463115091,1.731557546,1.499093481,5.194672637,\
                6.926230184,1.904561220,0.000000000,4.551627762,\
                3.809122439,0.,3.913707,5.738991])
            cbt=0
            fwav=['m2','s2','k1','o1','n2','p1','k2','q1','_2n2','mu2','nu2','l2',\
                  't2','j1','m1','oo1','rho1','mf','mm','ssa','m4','ms4','mn4','m6',\
                  'm8','mk3','s6','_2sm2','_2mk3','s1','_2q1','m3']

            for wav in fwav:
                exec(f'{wav}={val[cbt]}')
                cbt+=1
        
        # Time in Modified Julian Date
        # To compute lunar position time should be in days relatively Jan 1 2000 12:00:00
#        timetemp=cftime.date2num(date_ini,'days since 2000-01-01:12:00:00')
        timetemp=cftime.date2num(date_ini,'days since 1858-11-17:00:00:00')-51544.4993
        # mean longitude of lunar perigee
        P =  83.3535 +  0.11140353 * timetemp;
        P = np.mod(P,360.0)
#        if P<0:
#            P+360
        P=np.deg2rad(P)

        # mean longitude of ascending lunar node
        N = 125.0445 -  0.05295377 * timetemp
        N = np.mod(N,360.0)
#        if N<0:
#            N=N+360 
        N=np.deg2rad(N)

        # nodal corrections: pf = amplitude scaling factor [], 
        #                    pu = phase correction [rad]
        sinn = np.sin(N)
        cosn = np.cos(N)
        sin2n = np.sin(2*N)
        cos2n = np.cos(2*N)
        sin3n = np.sin(3*N)
        tmp1  = 1.36*np.cos(P)+.267*np.cos((P-N))
        tmp2  = 0.64*np.sin(P)+.135*np.sin((P-N))
        temp1 = 1.-0.25*np.cos(2*P)-0.11*np.cos((2*P-N))-0.04*cosn
        temp2 =    0.25*np.sin(2*P)+0.11*np.sin((2*P-N))+0.04*sinn
        ##
        if input_wav.lower() in ['sa','ssa','msf','alpha1','tau1','pi1','p1','s1','psi1','phi1','theta1','m2a','m2b','lambda2','la2','t2','s2','r2','s3','s4','s5','s6','s7','s8']:
            pf = 1.
        elif input_wav.lower() in ['2q1','sigma1','q1','rho1']:
            pf = np.sqrt((1.+.188*cosn)**2+(.188*sinn)**2)
        elif input_wav.lower() in ['2n2','mu2','n2','nu2','m2','2sm2']:
            pf = np.sqrt((1.-.03731*cosn+.00052*cos2n)**2 +(.03731*sinn-.00052*sin2n)**2)
        elif input_wav.lower() == 'mm':
            pf= 1 - 0.130*cosn   
        elif input_wav.lower() == 'mf':
            pf = 1.043 + 0.414*cosn
        elif input_wav.lower() == 'mt':
            pf = np.sqrt((one+.203*cosn+.040*cos2n)**2 + (.203*sinn+.040*sin2n)**2)
        elif input_wav.lower() == 'chi1':
            pf = np.sqrt((1.+.221*cosn)**2+(.221*sinn)**2)
        elif input_wav.lower() == 'k1':
            pf = np.sqrt((1.+.1158*cosn-.0029*cos2n)**2 + (.1554*sinn-.0029*sin2n)**2)
        elif input_wav.lower() == 'j1':
            pf = np.sqrt((1.+.169*cosn)**2+(.227*sinn)**2)
        elif input_wav.lower() == 'oo1':
            pf = np.sqrt((1.0+0.640*cosn+0.134*cos2n)**2 + (0.640*sinn+0.134*sin2n)**2 )
        elif input_wav.lower() == 'k2':
            pf = np.sqrt((1.+.2852*cosn+.0324*cos2n)**2 + (.3108*sinn+.0324*sin2n)**2)
        elif input_wav.lower() in ['mns2','mn4','m4','ms4']:
            pf = (np.sqrt((1.-.03731*cosn+.00052*cos2n)**2 + (.03731*sinn-.00052*sin2n)**2))**2
        elif input_wav.lower() == 'm6':
            pf = (np.sqrt((1.-.03731*cosn+.00052*cos2n)**2 + (.03731*sinn-.00052*sin2n)**2))**3
        elif input_wav.lower() == 'mk4':
            pf =np.sqrt((1.-.03731*cosn+.00052*cos2n)**2+(.03731*sinn-.00052*sin2n)**2)* \
                np.sqrt((1.+.2852*cosn+.0324*cos2n)**2 + (.3108*sinn+.0324*sin2n)**2)
        elif input_wav.lower() == 'mk3':
            pf = np.sqrt((1.+.1158*cosn-.0029*cos2n)**2 + (.1554*sinn-.0029*sin2n)**2) \
                *np.sqrt((1.-.03731*cosn+.00052*cos2n)**2 + (.03731*sinn-.00052*sin2n)**2)
        elif input_wav.lower() == 'eta2':
            pf = np.sqrt((1.+.436*cosn)**2+(.436*sinn)**2)
        elif input_wav.lower() == 'o1':
            pf = np.sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 + (0.189*sinn-0.0058*sin2n)**2)
        elif input_wav.lower() == 'm1':
            pf = np.sqrt(tmp1**2 + tmp2**2)
        elif input_wav.lower() == 'l2':
            pf = np.sqrt(temp1**2 + temp2**2)    
        else:
            print('No amplitude correction for wave %s' % input_wav)
            pf = 1
        
        ##
        if input_wav.lower() in ['sa','ssa','mm','msf','alpha1','tau1','pi1','p1','s1','psi1','phi1','theta1','m2a','m2b','lambda2','la2','t2','s2','r2','s3','s4','s5','s6','s7','s8']:
            pu=0
        elif input_wav.lower() in ['2q1','sigma1','q1','rho1']:
            pu = np.arctan(.189*sinn/(1.+.189*cosn))
        elif input_wav.lower() in ['2n2','mu2','n2','nu2','m2','2sm2','ms4']:
            pu = np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))
        elif input_wav.lower() in ['mns2','mn4','m4']:
            pu = 2*np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))
        elif input_wav.lower() == 'm3':
            pu = 1.5*np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))
        elif input_wav.lower() == 'mk3':
            pu = np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))+ \
                 np.arctan((-.1554*sinn+.0029*sin2n)/ (1.+.1158*cosn-.0029*cos2n))
        elif input_wav.lower() == 'mk4':
           pu = np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))+ \
                np.arctan(-(.3108*sinn+.0324*sin2n)/(1.+.2852*cosn+.0324*cos2n))
        elif input_wav.lower() == 'm6':
           pu = 3* np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))
        elif input_wav.lower() == 'k2':
           pu = np.arctan(-(.3108*sinn+.0324*sin2n)/(1.+.2852*cosn+.0324*cos2n))
        elif input_wav.lower() == 'k1':
           pu = np.arctan((-.1554*sinn+.0029*sin2n)/(1.+.1158*cosn-.0029*cos2n))
        elif input_wav.lower() == 'eta2':
           pu = np.arctan(-.436*sinn/(1.+.436*cosn))
        elif input_wav.lower() == 'j1':
           pu = np.arctan(-.227*sinn/(1.+.169*cosn))
        elif input_wav.lower() == 'oo1':
           pu = np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))
        elif input_wav.lower() == 'chi1':
           pu = np.arctan(-.221*sinn/(1.+.221*cosn))
        elif input_wav.lower() == 'm1':
           pu = np.arctan2(tmp2,tmp1) 
        elif input_wav.lower() == 'l2':
           pu = np.arctan(-temp2/temp1)
        elif input_wav.lower() == 'o1':
           pu = np.deg2rad(10.8*sinn - 1.3*sin2n + 0.2*sin3n)
        elif input_wav.lower() == 'mf':
           pu = np.deg2rad(-23.7*sinn + 2.7*sin2n - 0.4*sin3n)
        elif input_wav.lower() == 'mt':
           pu = np.arctan(-(.203*sinn+.040*sin2n)/(one+.203*cosn+.040*cos2n))
        else:
            print('No phase correction for wave %s' % input_wav)        
            pu=0
        
        try:
            mkb=eval(''.join(('phase_mkb.',input_wav.lower())))
        except:
            try:
                mkb=eval(''.join(('phase_mkb._',input_wav.lower())))
            except:
                print('No phase_mkB defined for wave %s' % input_wav)
                mkb=0

        return pf,pu,mkb

    class pot_tide():
        '''
        Includes the direct astronomical contribution from the sun and moon 
        (factor amp, from Schwiderski, 1978) and the contributions 
        from solid Earth body tide (factor elas, from Wahr, 1981).
        
        More informations about potential tides (or equilibrium tides) can be found here:
        - Global Ocean Tides. Part I , Schwiderski, 1978
        - Sea-level Science,p.43-44, David Pugh and Philip Woodworth
            
        '''

        amp=np.array([0.242334,0.112743,0.141565,0.100661,\
                      0.046397,0.046848,0.030684,0.019273,\
                      0.006141,0.007408,0.008811,0.006931,\
                      0.006608,0.007915,0.007915,0.004338,\
                      0.003661,0.042041,0.022191,0.019567,\
                      0.,0.,0.,0.,\
                      0.,0.,0.,0.,\
                      0.,7.6464e-04,0.002565,0.003192])

        elas=np.array([0.693,0.693,0.736,0.695,\
                       0.693,0.706,0.693,0.695,\
                       0.693,0.693,0.693,0.693,\
                       0.693,0.695,0.695,0.695,\
                       0.695,0.693,0.693,0.693,\
                       0.693,0.693,0.693,0.693,\
                       0.693,0.693,0.693,0.693,\
                       0.693,0.693,0.693,0.802])

        cbt=0
        fwav=['m2','s2','k1','o1','n2','p1','k2','q1','_2n2','mu2','nu2','l2',\
              't2','j1','m1','oo1','rho1','mf','mm','ssa','m4','ms4','mn4','m6',\
              'm8','mk3','s6','_2sm2','_2mk3','s1','_2q1','m3']

        for wav in fwav:
            exec(f'{wav}=[{amp[cbt]},{elas[cbt]}]')
            cbt+=1



