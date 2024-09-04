import numpy as np
import cftime
import pylab as plt
from dateutil.relativedelta import relativedelta
import xarray as xr
import _pickle as pickle
import scipy.signal as ss
from scipy.interpolate import interp1d
import matplotlib.pyplot as py
import sys,os
import copy
import warnings
warnings.filterwarnings("ignore")
##############################
def lonlat_to_m(lon,lat,plon,plat):
    '''
    This routine compute the distance in meter (on a sphere) between a grid and a specified point following Haversine formula

    Inputs:
      lon    Grid longitude
      lat    Grid latitude
      plon   Longitude of the reference point
      plat   Latitude of the reference point

    Outputs:
      dx     Grid distance to the reference point
    '''
    lon  = lon*2.*np.pi/360.
    lat  = lat*2.*np.pi/360.
    plon = plon*2*np.pi/360.
    plat = plat*2*np.pi/360
    if isinstance(lon,float):
        dx = np.arccos(np.sin(lat)*np.sin(plat) + np.cos(lat)*np.cos(plat)*np.cos(lon-plon))*6371000.
    else:
        dx = np.arccos(np.sin(lat[:])*np.sin(plat) + np.cos(lat[:])*np.cos(plat)*np.cos(lon[:]-plon))*6371000.
    return dx

def patch_adjust(lor, lar, lon, lat, mask, iz, jz,value=0):
    '''
    Recursive function allowing to find the closest masked/non-masked point depending on
    the integer 'value'

    Inputs:
      lor    Longitude of the reference point 
      lar    Latitude of the reference point
      lon    Grid longitude
      lat    Grid latitude
      mask   Grid mask
      iz     Reference point index on the grid in x direction
      jz     Reference point index on the grid in y direction
      value  if 0 looks for the closest non-masked point (if 1 closest masked point)
    Outputs:
      jz     Index of the closest (non)masked point in y direction
      iz     Index of the closest (non)masked point in x direction
    '''
    if value == 1:
        print('Look for the closest masked point ...')
    else:
        print('Look for the closest non-masked point ...')
    lono=np.copy(lon)
    lato=np.copy(lat)
    lon=np.ma.masked_where(mask==value,lon)
    lat=np.ma.masked_where(mask==value,lat)   
    
    dist = lonlat_to_m(lon,lat,lor,lar)      
    [jz,iz]=np.where(dist == np.nanmin(dist))
    if jz.shape[0]>1:
        jz=jz[0];iz=iz[0]
    if value == 1:
        jz,iz=patch_adjust(lono[int(jz),int(iz)], lato[int(jz),int(iz)], lono, lato, mask, int(iz), int(jz))
    return int(jz), int(iz)


def convtime(s):
    '''
    Try to convert a string to a date time
    
    Inputs:
      s    String to convert
    
    Outputs:
      datetime
    '''
    try:
        return plt.datetime.datetime.strptime(s,'%d/%m/%Y %H:%M:%S')
    except:
        try:
              return plt.datetime.datetime.strptime(s,'%Y/%m/%d %H:%M:%S')
        except:
            print('Could not find out date format from date %s' % s)
            sys.exit()


def nine_points_max_iter(Vin,Maskin):
    '''
    Iterative procedure to link all point of one river (in the coastal area)
    '''
    Vout=0*Vin
    res=100
    n=0
    while res>0:
        n=+1
        a=np.stack( (Vin[:-2,2:  ],Vin[1:-1,2:  ],Vin[2:,2:  ],\
                 Vin[:-2,1:-1],Vin[1:-1,1:-1],Vin[2:,1:-1],\
                 Vin[:-2,:-2 ],Vin[1:-1,:-2 ],Vin[2:,:-2 ] ),axis=2)
        Vin[1:-1,1:-1]=Maskin[1:-1,1:-1]*np.max(a,axis=2)
        res=np.max(Vin-Vout)
        Vout=Vin+0
        if n>15:
            print('nine_points_max_iter did not converge')
            break
     
    
#  print('nine_points_max_iter did converge after ',n,' iterations')
    Vout=Vout[np.where(Vout!=0)]

    return(Vout)



def read_river_netcdf(files,tstart,tend,time_units,Qmin=None,lon_pos=None,lat_pos=None):
    '''
    Read rivers in a netcdf files depending on what is given in Qmin,lon_pos,lat_pos.
    
    Inputs:
      files       Netcdf file to read
      tstart      Starting date
      tend        Ending date
      time_units  Reference time expressed as 'days since YYYY-MM-DD HH:mm:ss' 
      Qmin        Minimal flow condition to consider a river. Only used with 2-d
                  map of river flow (ex: https://cds.climate.copernicus.eu/cdsapp!/dataset/cems-glofas-historical)
      lon_pos     Longitude of the desired river. If None will look for a longitude var
      lat_pos     Latitude of the desired river. If None will look for a latitude var

    Ouputs:
      Depending on the input you can have different variables but the main are:
        Q         Flow of the rivers on the desired time period
        time      Date of each records
        lon       Longitude of rivers
        lat       Latitude of rivers
    '''
    data=xr.open_dataset(files)
    print("Asuming time variable is \'time\'")
    try:
        time=plt.date2num(data.time.values)
    except:
        try:
            time=plt.date2num(data.TIME.values)
        except:
            print('Could not find time variable... exit')
            sys.exit()
    # -- Selection pour la periode consideree
    rstrr = plt.date2num(tstart)
    rendr = plt.date2num(tend)
    interm = np.where((time >=rstrr) & (time<=rendr))
    if len(interm[0]) ==0:
        print('No data found between ',tstart, ' and ',tend, 'for %s' %files)
        sys.exit()

    # Check if data for previous month 
    prev_day = plt.date2num(tstart+relativedelta(months=-1))
    next_day = plt.date2num(tstart)
    ind_prev = np.where((time>=prev_day) & (time<next_day))
    if len(ind_prev[0])!=0:
        interm=(np.concatenate((ind_prev[0],interm[0])),)
    del prev_day,next_day,ind_prev    
    # Check if data for next month
    prev_day = plt.date2num(tend+relativedelta(hours=12))
    next_day = plt.date2num(tend+relativedelta(months=1,hours=12))
    ind_next = np.where((time>=prev_day) & (time<next_day))
    if len(ind_next[0])!=0:
        interm=(np.concatenate((interm[0],ind_next[0])),)
    del prev_day,next_day,ind_next
    dbt_var=[]
    lat_var=[];flip_lon=0
    lon_var=[];flip_lat=0
    data_units=['m3 s-1','m**3 s**-1','m3/s']
    for ii in data.data_vars: # Loop to find variables to use
        try:
            if data[ii].units in data_units:
                dbt_var.append(ii)
        except:
            pass
    for ii in data.coords:
        try:
            if data[ii].attrs['units'] in ['degrees_north','degree_north'] and lat_pos is None:
                lat_var.append(ii)
                if data[ii][1]-data[ii][0]<0: # if lat is decreasing
                    flip_lat=1
                     
            elif data[ii].attrs['units'] in ['degrees_east','degree_east'] and lon_pos is None:
                lon_var.append(ii)
                if data[ii][1]-data[ii][0]<0:
                    flip_lon=1
        except:
            pass
#########
    if len(dbt_var)>1:
        string=''
        for ii,jj in enumerate(dbt_var):
            string+='%s: %i, '%(jj,ii)
            choice=input('Several flow detected which one to choose?\n%s ?'%string[:-2])
            dbt_var=dbt_var[choice]
    elif len(dbt_var)==1:
        dbt_var=dbt_var[0]
    else:
        print('Did not found any variables with units', data_units,'... exit')
        sys.exit()
##############
    if len(lon_var)>1:  
        string=''
        for ii,jj in enumerate(lon_var):
            string+='%s: %i, '%(jj,ii)
            choice=input('Several longitude detected which one to choose?\n%s ?'%string[:-2])
            lon_var=lon_var[choice]
    elif len(lon_var)==1:
        lon_var=lon_var[0]
    elif len(lon_var)==0 and lon_pos is None:
        print('Did not found any lon variables with units:\'degrees_east\' ... exit')
        sys.exit()
#############
    if len(lat_var)>1:
        string=''
        for ii,jj in enumerate(lat_var):
            string+='%s: %i, '%(jj,ii)
            choice=input('Several longitude detected which one to choose?\n%s ?'%string[:-2])
            lat_var=lat_var[choice]
    elif len(lat_var)==1:
        lat_var=lat_var[0]
    elif len(lat_var)==0 and lat_pos is None:
        print('Did not found any lat variables with units:\'degrees_west\' ... exit')
        sys.exit()
#############

    if lon_pos is not None and lat_pos is not None: # Reading netcdf with lon lat specified
        print('Reading data for river %s' %files)
        Q=eval(''.join(('data.',dbt_var,'[interm[0],:].values')))
        return plt.num2date(time[interm]),Q,np.array([lon_pos]),np.array([lat_pos])

    elif Qmin is not None: # Reading netcdf with 2D data
        lon=np.squeeze(np.unique(eval(''.join(('data.',lon_var,'.values')))))
        lat=np.squeeze(np.unique(eval(''.join(('data.',lat_var,'.values')))))
        (Lon,Lat)=np.meshgrid(lon,lat)

        def_Fillval=False
        if "_FillValue" not in eval(''.join(('data.',dbt_var))).encoding:
            def_Fillval=True

        Q=eval(''.join(('data.',dbt_var,'[interm[0],:].values')))
        if flip_lat==1:
            Q=Q[:,::-1,:]
        if flip_lon==1:
            Q=Q[:,:,::-1]

        
        Qmean=np.squeeze(np.nanmean(Q,axis=0))

        # compute landmask
        mask = 0.*Qmean
        if def_Fillval:# If no FillValue in netcdf, assume 0 as value for the mask  
            Qmean[Qmean==0]=np.nan
        mask[np.where(np.isnan(Qmean)==False)]=1.
        mask[np.where(np.isnan(Qmean))]=0.

        Qmean=mask*Qmean
        lon_all=Lon[np.where(Qmean>=Qmin)]
        lat_all=Lat[np.where(Qmean>=Qmin)]
    
        # Get the coast mask
        (mx,my)=np.gradient(mask)
        coast=mask*np.sqrt(mx*mx+my*my)
        coast[np.where(coast!=0.)]=1.
        for i in range(1):
            # Broaden the coast mask to be sure to have all the rivers...
            coast[1:-1,1:-1]=ss.convolve(coast[1:-1,1:-1],np.ones((3,3)), mode='same')
            coast[np.where(coast>=1)]=1
            coast=coast*mask
  
        # Get the river mouths position
        Rmouth=Qmean*coast
        # Select the rivers with a mean discharge above a given treshold (Qmin)
        # each pixel with a value > Qmin is considered a specific river
        Rmouth[np.where(Rmouth<Qmin)]=0
        Rmouth[np.where(np.isnan(Rmouth))]=0
        Rmouth[np.where(Rmouth>=Qmin)]=1
        Rmouth[0,:]=0
        Rmouth[-1,:]=0
        Rmouth[:,0]=0
        Rmouth[:,-1]=0
        [Ny,Nx]=np.shape(Rmouth)
        # For each group of connected points selected as river mouths keep only one point using an iterative procedure 
        i=np.arange(Nx)
        j=np.arange(Ny)

        [I,J]=np.meshgrid(i,j)

        Imouth=I*Rmouth
        Jmouth=J*Rmouth

        i_mou=nine_points_max_iter(Imouth,Rmouth)
        j_mou=nine_points_max_iter(Jmouth,Rmouth)
        
        ij_mou=i_mou + 1e5*j_mou
        ij_mou=np.unique(ij_mou)

        j_mou=np.int32(np.floor(ij_mou*1e-5))
        i_mou=np.int32(ij_mou - 1e5*j_mou)
        # For each river mouth, get the connected location where the discharge is the largest
        Nrivers=np.size(i_mou)
        if Nrivers ==0:
            print('No rivers found leaving...')
            sys.exit()   
        else:
            print('Found %i rivers' %Nrivers)
               
        for i in range(Nrivers):
            Qtmp=Qmean[j_mou[i]-1:j_mou[i]+2,i_mou[i]-1:i_mou[i]+2]
            (jtmp,itmp)=np.where(Qtmp==np.nanmax(Qtmp))
            jtmp=jtmp[0]
            itmp=itmp[0]
            j_mou[i]=j_mou[i]+jtmp-1
            i_mou[i]=i_mou[i]+itmp-1

        Q_mou=Q[:,j_mou,i_mou]
        lon_mou=Lon[j_mou,i_mou]
        lat_mou=Lat[j_mou,i_mou]
        return Nrivers,plt.num2date(time[interm]),Q_mou,lon_mou,lat_mou 

    else: # Reading data with multiple rivers
        lon=eval(''.join(('data.',lon_var,'.values')))
        lat=eval(''.join(('data.',lat_var,'.values')))
        Q=eval(''.join(('data.',dbt_var,'[interm[0],:].values')))
        print('Found %i rivers in %s' %(Q.shape[1],files))
        
        return Q.shape[1], plt.num2date(time[interm]), Q, lon, lat            



def read_river(list_river_files,lon_inp,lat_inp,rstr,rend,time_units):
    ''' 
    Read river flows files. It handle netcdf,.txt or .dat files
    
    Inputs:
      list_river_files     List of all river files
      lon_inp              Longitude of the river. To read netcdf value put None
      lat_inp              Latitude of the river OR minimal river flow to consider.
                           To read netcdf value put None
      rstr                 Starting date
      rend                 Ending date
      time_units           Reference time expressed as 'days since YYYY-MM-DD HH:mm:ss'

    Outputs:
      river                Dictionnary with all river intels. Each river dictionnary has
                           for key 'time','flow','longitude','latitude'
    '''
    # -- Obs
    river=dict()
    riv_cnt=0;Nriv=1
    # -- Boucle sur les fichiers
    for ific,fic in enumerate(list_river_files):
        # - Compteur
        if '.nc' in fic:
            if np.isnan(lon_inp[ific])  and np.isnan(lat_inp[ific]):
            # Case of a netcdf with multiple station
                Nriv,time,flw,lon,lat = read_river_netcdf(fic,rstr,rend,time_units)
            elif np.isnan(lon_inp[ific]):
            # case of a 2D map 
                print('Script will look for rivers with flow >',lat_inp[ific],'m3/s')
                Nriv,time,flw,lon,lat = read_river_netcdf(fic,rstr,rend,time_units,Qmin=lat_inp[ific])
            else:
            # Case of a netcdf for 1 river with lon/lat specified
                time,flw,lon,lat = read_river_netcdf(fic,rstr,rend,time_units,lon_pos=lon_inp[ific],lat_pos=lat_inp[ific])

            for k in range(Nriv): 
                river['river_'+str(riv_cnt)]=dict()
                river['river_'+str(riv_cnt)]['time']=cftime.date2num(time,time_units)
                river['river_'+str(riv_cnt)]['flow']=np.squeeze(flw[:,k])
                river['river_'+str(riv_cnt)]['longitude']=np.squeeze(lon[k])
                river['river_'+str(riv_cnt)]['latitude']=np.squeeze(lat[k])
                riv_cnt+=1
            Nriv=1
        elif '.dat' or '.txt' in fic: 
            # -- Lecture du fichier
            file_data = open(fic).read()
            lline = file_data.split('\n')
            # ---- First scan ----
            detect=0
            for nline,ll in enumerate(lline):
                if ('$' not in ll) and (detect==0):
                    try:
                        d = convtime(ll[:19])
                        detect=1
                        first_data=nline
                    except:
                        continue
            # ---- Lecture des donnees
#            data=np.genfromtxt(file_data,dtype=None,converters={0:convtime},delimiter=[19,20],skip_header=first_data)
            data=np.genfromtxt(fic,dtype=str,skip_header=first_data,delimiter=[19,20],comments='#')
   
            time = np.array([convtime(x[0]) for x in data])
            flw  = np.array([float(x[1]) for x in data])

            # -- Selection pour la periode consideree
            tr = [cftime.date2num(x,time_units) for x in time]
            tr = np.array(tr)
            rstrr = cftime.date2num(rstr,time_units)
            rendr = cftime.date2num(rend,time_units)
            interm = np.where((tr >=rstrr) & (tr<=rendr))
            if len(interm) ==0:
                print('No data found between ',rstr, ' and ',rend, 'for %s' %fic)
                sys.exit()

            # Check if data for previous month
            prev_day = rstr+relativedelta(months=-1)
            next_day = rstr
            ind_prev = np.where((tr>=cftime.date2num(prev_day,time_units)) & (tr<cftime.date2num(next_day,time_units)))
            if len(ind_prev[0])!=0:
                interm=(np.concatenate((ind_prev[0],interm[0])),)

           # Check if data for next day
            prev_day = rend+relativedelta(hours=12)
            next_day = rend+relativedelta(months=1,hours=12)
            ind_prev = np.where((tr>=cftime.date2num(prev_day,time_units)) & (tr<cftime.date2num(next_day,time_units)))
            if len(ind_prev[0])!=0:
                interm=(np.concatenate((interm[0],ind_prev[0])),)

            flw  = flw[interm]
            time = time[interm]

            river[os.path.basename(fic)]=dict()
            river[os.path.basename(fic)]['time']=cftime.date2num(time,time_units)
            river[os.path.basename(fic)]['flow']=flw
            river[os.path.basename(fic)]['longitude']=lon_inp[ific]
            river[os.path.basename(fic)]['latitude']=lat_inp[ific]

#            # --- Sauvegarde des fichiers lus ---
#            pickle.dump(river, open('total_river_'+plt.datetime.datetime.strftime(rdeb,'%Y%m%d')+'.p','wb'))
#        else:
#            river=pickle.load(open('total_river_'+plt.datetime.datetime.strftime(rdeb,'%Y%m%d')+'.p','rb'))

    return river


def fill_period(river,rstr,rend,time_units,output_type):
    # Function to put all rivers on a daily scale
    l=river.keys()
    rstr_sec=cftime.date2num(rstr,time_units)*86400+43200 #start at 12:00:00 on first day
    rend_sec=cftime.date2num(rend,time_units)*86400
    if output_type.upper()=="HOURLY":
        new_time=np.arange(rstr_sec-46800,rend_sec+50400,3600)
    elif output_type.upper()=='DAILY':
        new_time=np.arange(rstr_sec-86400,rend_sec+129600,86400)
    elif output_type.upper()=='MONTHLY':
        new_time=np.array([])
        start_mth = rstr+relativedelta(hours=12)
        end_mth   = start_mth+relativedelta(months=1,days=-1)
        while (cftime.date2num(end_mth,time_units) <= cftime.date2num(rend,time_units)):
            mid_mth=(cftime.date2num(start_mth,time_units)+cftime.date2num(end_mth,time_units))/2
            if len(new_time)==0:
                prev_mid_mth=cftime.date2num(rstr,time_units)-(mid_mth-cftime.date2num(rstr,time_units))
                new_time = np.append(new_time,prev_mid_mth)
            new_time = np.append(new_time,mid_mth)
            #
            start_mth = start_mth+relativedelta(months=1)
            end_mth   = start_mth+relativedelta(months=1,days=-1)
        next_mid_mth=(cftime.date2num(start_mth,time_units)+cftime.date2num(end_mth,time_units))/2
        new_time = np.append(new_time,next_mid_mth)
        new_time = new_time*86400 # convert in second)
    ##
    for ir,r in enumerate(l):
        t=np.squeeze(river[r]['time'][:]*86400)
        f=np.squeeze(river[r]['flow'][:])
        if output_type.upper()=='DAILY' and np.nanmean(np.diff(t))<86400:
            # Need to average hourly data on daily
            f_cor=np.zeros(len(new_time))
            for nn in range(1,len(new_time)-1):
                index=np.where((t>new_time[nn]-43200) & (t<new_time[nn]+43200))
                f_cor[nn]=np.nanmean(f[index])

            # look for previous data
            ind_prev=np.where((t>new_time[0]-43200) & (t<new_time[0]+43200))
            if len(ind_prev[0])==0:
                f_cor[0]=f_cor[1]
            else:
                f_cor[0]=np.nanmean(f[ind_prev])
            del ind_prev
            # look for next data
            ind_next=np.where((t>new_time[-1]-43200) & (t<new_time[-1]+43200))
            if len(ind_next[0])==0:
                f_cor[-1]=f_cor[-2]
            else:
                f_cor[0]=np.nanmean(f[ind_next])
            del ind_next
            
        elif output_type.upper()=='HOURLY' and np.nanmean(np.diff(t))>3600:
            print('Output asked on HOURLY frequency but %s seems to have lower frequency. Performing linear interpolation' % r)
            interp=interp1d(t,f,bounds_error=False,fill_value=(f[0],f[-1]))
            f_cor=interp(new_time)
        elif output_type.upper()=='MONTHLY':
            f_cor=np.zeros(len(new_time))
            start_mth = rstr
            end_mth   = start_mth+relativedelta(months=1,hours=-12)
            # check if values for previous month
            prev_day = cftime.date2num(start_mth+relativedelta(months=-1),time_units)*86400
            next_day = cftime.date2num(start_mth,time_units)*86400
            ind_prev = np.where((t>=prev_day) & (t<next_day))
            if len(ind_prev[0])==0:
               prev=0
            else:
               prev=1
               f_cor[0]=np.nanmean(f[ind_prev])
            
            del prev_day,next_day,ind_prev
            cpt=1
            while (cftime.date2num(end_mth,time_units)<=cftime.date2num(rend,time_units)):
                index=np.where((t>=cftime.date2num(start_mth,time_units)*86400) \
                           & (t<cftime.date2num(end_mth,time_units)*86400 ))
                f_cor[cpt]=np.nanmean(f[index])
                cpt+=1          
                start_mth = start_mth+relativedelta(months=1)
                end_mth   = start_mth+relativedelta(months=1,hours=-12)
            # check if values for next month
            prev_day = cftime.date2num(start_mth,time_units)*86400
            next_day = cftime.date2num(end_mth,time_units)*86400
            ind_next = np.where((t>=prev_day) & (t<next_day))
            if len(ind_next[0])==0:
               next=0
            else:
               next=1
               f_cor[-1]=np.nanmean(f[ind_next])
            # If no value for previous or next month, fill with first and last month
            if prev==0:
                f_cor[0]=f[1]
            if next==0:
                f_cor[-1]=f_cor[-2]
        else:
            interp=interp1d(t,f,bounds_error=False,fill_value=(f[0],f[-1]))
            f_cor=interp(new_time)
        
        river[r]['time']=new_time/86400
        river[r]['flow']=f_cor

    return river

def get_river_index(river, croco_class):
    '''
    Routine to compute river first guess position on croco grid
  
    Inputs:
      river           River dictionnary
      croco_class     CROCO class including croco grid
    
    Outputs:
      river           River dictionnary with new keys ii,jj indicatingfirst guess position                      on croco grid
    '''    
    # Keep rivers in domain and exclude others
    river_out=dict()
    geolim=[croco_class.lonmin(),croco_class.lonmax(),croco_class.latmin(),croco_class.latmax()]
    for ll in river.keys():
        # Eventually need to adapt lon/lat format (-180,180)/(0,360)      
        if (river[ll]['longitude']>=(geolim[0]-360)) & (river[ll]['longitude']<=(geolim[1]-360)):
            river[ll]['longitude']+=360
        elif (river[ll]['longitude']>=geolim[0]+360) & (river[ll]['longitude']<=geolim[1]+360) :
            river[ll]['longitude']-=360

        lor = river[ll]['longitude']
        lar = river[ll]['latitude']
        if (lor >= geolim[0] ) & (lor <= geolim[1]) &\
           (lar >= geolim[-2]) & (lar <= geolim[-1]):
            river_out[ll]=river[ll]

            dist = lonlat_to_m(croco_class.lon,croco_class.lat,lor,lar)
            [jpos,ipos]=np.where(dist == np.min(dist))
            jpos,ipos=int(jpos),int(ipos) 
                                         # Init Kernel of shape
            kernel=np.arange(9)%2        #  0 1 0
            kernel=kernel.reshape((3,3)) #  1 0 1
                                         #  0 1 0
            check= ss.convolve(croco_class.maskr[jpos-1:jpos+2,ipos-1:ipos+2], kernel, mode='same')[1,1] # check if point is close to mask
            
            if croco_class.maskr[jpos,ipos] == 0:
                jpos,ipos = patch_adjust(lor, lar, croco_class.lon,croco_class.lat, croco_class.maskr, ipos, jpos)
            elif croco_class.maskr[jpos,ipos] == 1 and check==4: #mean only surrounded with water, need at least 1 land point around
                print('Moving point, only surounded with water')
                jpos,ipos = patch_adjust(lor, lar, croco_class.lon,croco_class.lat, croco_class.maskr, ipos, jpos,value=1)
            
            river_out[ll]['ii']=float(ipos)+1 # +1 To follow fortran indexing.....
            river_out[ll]['jj']=float(jpos)+1
            
            print('Choose a default value for runoff direction')
            kernel=np.zeros((3,3,4))
            
            kernel[:,:,0]=np.array(([0,1,0],[0,0,0],[0,0,0])) # look if land up
            kernel[:,:,1]=np.array(([0,0,0],[1,0,0],[0,0,0])) # look if land right
            kernel[:,:,2]=np.array(([0,0,0],[0,0,1],[0,0,0])) # look if land left
            kernel[:,:,3]=np.array(([0,0,0],[0,0,0],[0,1,0])) # look if land down
            mask_inv=np.zeros((3,3))
            mask_inv[croco_class.maskr[jpos-1:jpos+2,ipos-1:ipos+2] ==0 ] =1
            land=np.array([ss.convolve(mask_inv, kernel[:,:,i], mode='same')[1,1] for i in range(kernel.shape[-1]) ])
            if (land == np.array([0,1,0,0])).all(): # Land only at the right 
                river_out[ll]['dsrc']=0
                river_out[ll]['qbardir']=-1
            elif (land == np.array([0,0,1,0])).all(): # Land only at the left
                river_out[ll]['dsrc']=0
                river_out[ll]['qbardir']=1
                river_out[ll]['ii']=float(ipos)
            elif (land == np.array([1,0,0,0])).all(): # Land only at the top
                river_out[ll]['dsrc']=1
                river_out[ll]['qbardir']=-1
            elif (land == np.array([0,0,0,1])).all(): # Land only at the bottom
                river_out[ll]['dsrc']=1
                river_out[ll]['qbardir']=1
                river_out[ll]['jj']=float(jpos)
            elif (land == np.array([1,0,0,1])).all(): # Land at top and bottom
                river_out[ll]['dsrc']=1
                river_out[ll]['qbardir']=1
            elif (land == np.array([0,1,1,0])).all(): # Land at left and right
                river_out[ll]['dsrc']=0
                river_out[ll]['qbardir']=1
            elif (land == np.array([1,1,0,0])).all() or (land == np.array([0,1,0,1])).all(): # Land at the top/bottom and right
                river_out[ll]['dsrc']=0
                river_out[ll]['qbardir']=-1
            elif (land == np.array([0,0,1,1])).all() or (land == np.array([1,0,1,0])).all(): # Land at the top/bottom and left
                river_out[ll]['dsrc']=0
                river_out[ll]['qbardir']=1
                river_out[ll]['ii']=float(ipos)
            elif (land == np.array([1,1,1,0])).all(): # Bay oriented southward
                river_out[ll]['dsrc']=1
                river_out[ll]['qbardir']=-1
            elif (land == np.array([0,1,1,1])).all(): # Bay oriented northward
                river_out[ll]['dsrc']=1
                river_out[ll]['qbardir']=1
                river_out[ll]['jj']=float(jpos)
            elif (land == np.array([1,1,0,1])).all(): # Bay oriented eastward
                river_out[ll]['dsrc']=0
                river_out[ll]['qbardir']=1
                river_out[ll]['ii']=float(ipos)
            elif (land == np.array([1,0,1,1])).all(): # Bay oriented westward
                river_out[ll]['dsrc']=0
                river_out[ll]['qbardir']=1
                


     
        else:
            print('Remove %s , river mouth (%.02f,%.02f) not in the domain:%.02f,%.02f/%.02f,%.02f ' % (ll,lor,lar\
  ,geolim[0],geolim[1],geolim[-2],geolim[-1]))

    return river_out



def locate_runoff(direc,j,i,mask,masku,maskv):
    '''
     Function derived from locate_runoff.m used for make_runoff.m in ROMSTOOLS
    '''

    j=j-1
    i=i-1


    if mask[int(j),int(i)]==1:
        print( 'River positionned in sea')
        insea=1
    else:
        print( 'River positionned in land')
        insea=0
    if direc[0]==0:
        if insea==1:
            if direc[1]==-1:      #est - west => TESTED
                while masku[int(j),int(i)]==1:
                    i=i+1
                    # disp(['i:',num2str(i)])
                    # disp(['MASKU:',num2str(masku(j,i))])

                # disp(['--'])
                # disp(['i:',num2str(i)])
                # disp(['MASKU:',num2str(masku(j,i))])
            elif direc[1]==1: # west - est => TESTED
                while masku[int(j),int(i)]==1:
                    i=i-1
                    # disp(['MASKU:',num2str(masku(j,i))])

                # disp(['--'])
                # disp(['MASKU:',num2str(masku(j,i))])

        else: #inland
            if direc[1]==1:      #west-est  => TESTED
                while masku[int(j),int(i)]!=1:
                    i=i+1
                    # disp(['MASKU:',num2str(masku(j,i))])

                # disp(['--'])
                i=i-1
                # disp(['MASK:',num2str(masku(j,i))])
            elif direc[1]==-1: #est-west => TESTED
                while masku[int(j),int(i)]!=1:
                    i=i-1
                    # disp(['MASKU:',num2str(masku(j,i))])

                # disp(['--'])
                i=i+1
    else: # direc(k,1)=1
        if insea==1:
            if direc[1]==-1:     # nord - sud  => TESTED
                while maskv[int(j),int(i)]==1:
                    j=j+1
                    # disp(['MASKV:',num2str(maskv(j,i))])

                # disp(['--'])
                # disp(['MASKV:',num2str(maskv(j,i))])
            elif direc[1]==1: # sud - nord => TESTED
                while maskv[int(j),int(i)]==1:
                    j=j-1
                    # disp(['MASKV:',num2str(maskv(j,i))])

                # disp(['--'])
                # disp(['MASKV:',num2str(maskv(j,i))])

        else: # inland  
            if direc[1]==1:      #sud-nord  => TESTED          
                while maskv[int(j),int(i)]!=1:
                    j=j+1
                    # disp(['MASKV:',num2str(maskv(j,i))])

                # disp(['--'])
                j=j-1
                # disp(['MASKV:',num2str(maskv(j,i))])
            elif direc[1]==-1: #nord-sud => TESTED
                while maskv[int(j),int(i)]!=1:
                    j=j-1
                    # disp(['MASKV:',num2str(maskv(j,i))])

                # disp(['--'])
                j=j+1
                # disp(['MASKV:',num2str(maskv(j,i))])
    j2=j
    i2=i
   # Convert from j,i for array 0:.... (Python) to j,i for array 1:... (ROMS)
    j2=j2+1
    i2=i2+1

    return int(j2), int(i2)

def write_croco_in(river, croco_dir,riverfile,filw='for_croco_in.txt'):
    '''
    Write a text file to help fill in croco.in
    
    Inputs:
      river       River dictionnary
      croco_dir   Directory to put the text file
      riverfile   Netcdf file with interranual forcing
      filw        Name of the text file created ( default: for_croco_in.txt)
    '''

    fw = open(croco_dir+filw, "w")
    fw.write(' \n')
    fw.write('Line to enter in the croco.in file in the psource_ncfile section :\n')
    fw.write('-----------------------------------------------------------------\n')
    fw.write('psource_ncfile:   Nsrc  Isrc  Jsrc  Dsrc qbardir  Lsrc  Tsrc   runoff file name\n')
    fw.write('                           %s\n'%(croco_dir+riverfile))
    fw.write('                   '+str(len(river))+'\n')
    for il,ll in enumerate(river.keys()):
        fw.write('                       %5i %5i   %1i      %2.1d    T T  14.0  5.0\n' % (river[ll]['ii'],river[ll]['jj'],river[ll]['dsrc'],river[ll]['qbardir']))
    fw.write(' \n')
    fw.write('Line to enter in the croco.in file in the psource section :\n')
    fw.write('-----------------------------------------------------------------\n')
    fw.write('psource:   Nsrc  Isrc  Jsrc  Dsrc Qbar [m3/s]    Lsrc        Tsrc\n')
    fw.write('            '+str(len(river))+'\n')
    for il,ll in enumerate(river.keys()):
        fw.write('                %5i %5i   %1i   %6.1d          T T      14.   5.\n' % ( river[ll]['ii'],river[ll]['jj'],river[ll]['dsrc'],river[ll]['qbardir']*np.nanmean(river[ll]['flow'])))


    fw.close()

    return



def locateji_croco(river,grd=None,graph=False,subreg=''):
    '''
    Correction of river indexes
    
    Inputs:
      river    River dictionnary
      grd      Croco grid class
      graph    True, plot the position with mask before and after correction
    
    Outputs
      river    River dictionnary with updated indexes
    '''

    print('Correction of the position of the rivers ...')

    # list_river=[]
    # for ll in river.keys():
    #     list_river.append(ll)

    lon = grd.lon
    lat = grd.lat
    lonu = grd.lonu
    latu = grd.latu
    lonv = grd.lonv
    latv = grd.latv
    mask = grd.maskr
    masku = grd.umask
    maskv = grd.vmask

    river_e=dict()

    if graph:
        print( '-- Plot --')
        fig, ax = py.subplots(1, 1)
        ax.pcolor(masku,edgecolor='w')


        fig2, ax2 = py.subplots(1, 1)
        ax2.pcolor(maskv,edgecolor='w')

        fig3, ax3 = py.subplots(1, 1)
        ax3.pcolor(lonu,latu,masku,edgecolor='w')

        fig4, ax4 = py.subplots(1, 1)
        ax4.pcolor(lonv,latv,maskv,edgecolor='w')

    print('-- Optimisation du positionnement --')
    for ll in river.keys():
        print(ll)
        # ll = ll.split()

        # if len(ll) > 0:
        #
        #      if ll[2] in list_river:

        # print ll[2]

        river_e[ll]=river[ll]

        print( 'i: ',river_e[ll]['ii'])
        print( 'j: ',river_e[ll]['jj'])

        direc=np.array([river[ll]['dsrc'],river[ll]['qbardir']])

        [j2,i2]=locate_runoff(direc,river_e[ll]['jj'],river_e[ll]['ii'],mask,masku,maskv)

        if graph:

            ax.plot(river_e[ll]['ii']-1,river_e[ll]['jj']-1,'or')
            ax.plot(i2-1,j2-1,'+c')

            ax2.plot(river_e[ll]['ii']-1,river_e[ll]['jj']-1,'or')
            ax2.plot(i2-1,j2-1,'+c')
            # print river_e[ll]['jj']-1,river_e[ll]['ii']-1
            ax3.plot(lon[int(river_e[ll]['jj']-1),int(river_e[ll]['ii']-1)],lat[int(river_e[ll]['jj']-1),int(river_e[ll]['ii']-1)],'or')
            ax3.plot(lon[j2-1,i2-1],lat[j2-1,i2-1],'+c')
            ax3.plot(river_e[ll]['longitude'],river_e[ll]['latitude'],'+b')

            ax4.plot(lon[int(river_e[ll]['jj']-1),int(river_e[ll]['ii']-1)],lat[int(river_e[ll]['jj']-1),int(river_e[ll]['ii']-1)],'or')
            ax4.plot(lon[j2-1,i2-1],lat[j2-1,i2-1],'+c')
            ax4.plot(river_e[ll]['longitude'],river_e[ll]['latitude'],'+b')

        print( 'i2: ',i2)
        print( 'j2: ',j2)


        river_e[ll]['ii']=i2
        river_e[ll]['jj']=j2

    if graph:
        ax.set_title('mask_u')
        ax2.set_title('mask_v')
        # mp.legend('before','after')
        py.show()

    return river_e



def qtriver(maskr,masku, maskv,lonr,latr, lonu, lonv, latu, latv, lonriver, latriver, rivername, ii, jj, dsrc, qbardir):
    

    # ------------------
    
    # ------------------
    # Adapte de: 
    # embedding_in_qt5.py --- Simple Qt5 application embedding matplotlib canvases
    #
    # Copyright (C) 2005 Florent Rougon
    #               2006 Darren Dale
    #               2015 Jens H Nielsen
    #
    # This file is an example program for matplotlib. It may be used and
    # modified with no restriction; raw copies as well as modified versions
    # may be distributed without limitation.

    # from __future__ import unicode_literals
    import sys
    import os
    import random
    import matplotlib
    from PyQt5 import QtCore, QtWidgets
    from PyQt5.QtWidgets import  QPushButton, QSizePolicy, QLabel, QLineEdit,QComboBox,QMessageBox
    from PyQt5.QtCore import pyqtSlot, Qt

    from numpy import arange, sin, pi
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.figure import Figure

    # Make sure that we are using ST5
    #matplotlib.use('Qt5Agg')



    progname = os.path.basename(sys.argv[0])
    progversion = "0.1"
  
    # --- Stockage des arguments dans un dictionnaire:
    rr=dict()
    rr['maskr']=maskr
    rr['masku']=masku
    rr['maskv']=maskv
    rr['lonr']=lonr
    rr['lonu']=lonu
    rr['lonv']=lonv
    rr['latr']=latr
    rr['latu']=latu
    rr['latv']=latv
    rr['lonriver']=lonriver
    rr['latriver']=latriver
    rr['rivername']=rivername
    rr['ii']=int(ii)
    rr['jj']=int(jj)
    rr['dsrc']=dsrc
    rr['qbardir']=qbardir
    rr['newname']=np.array([rivername])
    rr['inew']=np.array([int(ii)])
    rr['jnew']=np.array([int(jj)])
    rr['dsrcnew']=np.array([dsrc])
    rr['qbardirnew']=np.array([qbardir])
    rr['stop']=False

    class MaskProcess(FigureCanvas):

        def __init__(self, parent=None, width=5, height=5, dpi=200):
            
            self.fig = Figure(figsize=(width, height), dpi=dpi)
           
            FigureCanvas.__init__(self, self.fig)
            self.setParent(parent)
            FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
            FigureCanvas.updateGeometry(self)
            self.fig.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
            self.fig.canvas.setFocus()
            self.locworkind=0
            self.locnbriv=1
            # self.addToolBar(QtCore.Qt.BottomToolBarArea,NavigationToolbar(dynamic_canvas, self))
            # self.cid=fig.canvas.mpl_connect('key_press_event', self.process_key)
            #self.cid=fig.canvas.mpl_connect('button_release_event', self.zoom)
            self.fig.canvas.mpl_connect('button_press_event', self.on_click)
            self.input_plot()
            self.point_plot()

        def input_plot(self):
            # volume = io.imread("14.mha", plugin='simpleitk')
            self.ax = self.fig.add_subplot(111)
            self.ax.pcolormesh(rr['lonr'],rr['latr'],rr['maskr'],zorder=0)#,edgecolor='black',lw=0.01)
            self.ax.set_title(rr['rivername'])
            self.draw()
 
        def point_plot(self):
            self.ax.plot(rr['lonriver'],rr['latriver'],'+g')
            if rr['dsrcnew']==0:
                self.ax.plot(rr['lonu'][rr['jj']-1,rr['ii']-1],rr['latu'][rr['jj']-1,rr['ii']-1],'or')           
                self.ax.quiver(rr['lonu'][rr['jj']-1,rr['ii']-1],rr['latu'][rr['jj']-1,rr['ii']-1],[rr['qbardir']],[0])
            elif rr['dsrcnew']==1:
                self.ax.plot(rr['lonv'][rr['jj']-1,rr['ii']-1],rr['latv'][rr['jj']-1,rr['ii']-1],'or')
                self.ax.quiver(rr['lonv'][rr['jj']-1,rr['ii']-1],rr['latv'][rr['jj']-1,rr['ii']-1],[0],[rr['qbardir']])
            self.draw()
    
#        def zoom(self,event):
#           #  Pour eviter que l'action lorsqu'on zoom ne fasse aussi comme si on cliquait
#            self.press=None
#            self.draw()

        def on_click(self, event):
            if event.button == 3:
#                print('on_click')
#                print('%s click: button=%d, xdata=%f, ydata=%f' %
#                      ('double' if event.dblclick else 'single', event.button,
#                       event.xdata, event.ydata))
                dist = lonlat_to_m(rr['lonr'],rr['latr'],event.xdata,event.ydata)
                [self.jz,self.iz]=np.where(dist == np.min(dist)) 
            
                # Convert from j,i for array 0:.... (Python) to j,i for array 1:... (ROMS)
                self.jz = self.jz[0]+1
                self.iz = self.iz[0]+1
                if len(self.ax.lines) > self.locnbriv+1:
                    self.ax.lines[-1].remove() #Del previous pont on plot
                    self.ax.collections[-1].remove() # Del associate quiver

                if rr['dsrcnew'][self.locworkind]==0:
                    self.ax.plot(rr['lonu'][self.jz-1,self.iz-1],rr['latu'][self.jz-1,self.iz-1],'oc')
                    self.ax.quiver(rr['lonu'][self.jz-1,self.iz-1],rr['latu'][self.jz-1,self.iz-1],rr['qbardirnew'][self.locworkind],0)
                elif rr['dsrcnew'][self.locworkind]==1:
                    self.ax.plot(rr['lonv'][self.jz-1,self.iz-1],rr['latv'][self.jz-1,self.iz-1],'oc')
                    self.ax.quiver(rr['lonv'][self.jz-1,self.iz-1],rr['latv'][self.jz-1,self.iz-1],0,rr['qbardirnew'][self.locworkind])

                rr['inew'][self.locworkind]=self.iz
                rr['jnew'][self.locworkind]=self.jz
            self.draw()

    class ApplicationWindow(QtWidgets.QMainWindow):
        def __init__(self):
            QtWidgets.QMainWindow.__init__(self)
            # self.left = 10
            #  self.top = 10
            #  self.width = 600
            #  self.height = 400
            self.showmesh=0
            self.countriv=1
            self.workind=0

            self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
            self.setWindowTitle("application main window")
            # self.setGeometry(self.left, self.top, self.width, self.height)

            self.file_menu = QtWidgets.QMenu('&Folder', self)
            self.file_menu.addAction('&Exit', self.fileQuit,
                                     QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
            self.menuBar().addMenu(self.file_menu)

            self.help_menu = QtWidgets.QMenu('&Help', self)
            self.menuBar().addSeparator()
            self.menuBar().addMenu(self.help_menu)

            self.help_menu.addAction("&How to use", self.about)

            self.main_widget = QtWidgets.QWidget(self)
            

            
            # creation du premier bouton
            self.quitsav = QPushButton("Save chosen position and Leave")
            self.quitwosav = QPushButton("Leave without modification")
            self.stop = QPushButton("Stop corrections")

            l = QtWidgets.QVBoxLayout(self.main_widget)
            

            hbox1 = QtWidgets.QHBoxLayout()
            hbox1.addStretch(1)

            self.dd = MaskProcess(self.main_widget, width=10, height=4)
            self.addToolBar(NavigationToolbar(self.dd, self))

            hbox1.addWidget(self.dd)         

            vbox1 = QtWidgets.QVBoxLayout()
            self.meshbox=QPushButton("Show/Hide mesh")
            self.add_riv=QPushButton("Add point for current river")
            self.del_riv=QPushButton("Remove point")
            self.rivpoint=QComboBox()
            self.rivpoint.addItems([rr['rivername']])
            vbox1.addStretch()
            vbox1.addWidget(self.meshbox)
            vbox1.addWidget(self.add_riv)
            vbox1.addWidget(self.rivpoint)
            vbox1.addWidget(self.del_riv)
            vbox1.addStretch()
            hbox1.addLayout(vbox1)
            l.addLayout(hbox1)
           
            self.meshbox.clicked.connect(self.on_mesh_click)
            self.add_riv.clicked.connect(self.on_addriv_click)
            self.del_riv.clicked.connect(self.on_delriv_click)
            self.rivpoint.activated.connect(self.on_selectriv)
            ###
            hbox2 = QtWidgets.QHBoxLayout()
            hbox2.addStretch(1)          
            self.label1 = QLabel("dsrc (0 - East-West; 1 - South-North): ")
            self.textbox = QLineEdit(str(rr['dsrc']))
            self.label2 = QLabel("qbardir ( 1 is positive [S-N or W-E], -1 negative [N-S or E-W] ): ")
            self.textbox2 = QLineEdit(str(rr['qbardir']))
            self.dirbutton = QPushButton("Update (dsrc / qbardir)")
            
            hbox2.addWidget(self.label1)
            hbox2.addWidget(self.textbox)
            hbox2.addWidget(self.label2)
            hbox2.addWidget(self.textbox2)
            hbox2.addWidget(self.dirbutton)           
            l.addLayout(hbox2)
            ###          
            hbox = QtWidgets.QHBoxLayout()
            hbox.addStretch(1)
            hbox.addWidget(self.quitsav)
            hbox.addWidget(self.quitwosav)
            hbox.addWidget(self.stop)

            l.addLayout(hbox)
            ###
            
            self.dirbutton.clicked.connect(self.on_click_dir)
            self.quitsav.clicked.connect(self.on_click)
            self.quitwosav.clicked.connect(self.on_click_wo)
            self.stop.clicked.connect(self.on_click_stop)

            self.main_widget.setFocus()
            self.setCentralWidget(self.main_widget)

            # self.statusBar().showMessage("All hail matplotlib!", 2000)

        def fileQuit(self):
#            self.close()
            qApp.quit()
        def closeEvent(self, ce):
            self.fileQuit()
            
        def on_click(self):
            print('Saving ...')
            rr['dsrcnew'][self.workind]=int(self.textbox.text())
            rr['qbardirnew'][self.workind]=int(self.textbox2.text())
            # print rr['inew'], rr['jnew']
            self.fileQuit()
                      
        def on_click_wo(self):
            rr['inew'][self.workind]=rr['ii']
            rr['jnew'][self.workind]=rr['jj']
            rr['dsrcnew'][self.workind]=rr['dsrc']
            rr['qbardirnew'][self.workind]=rr['qbardir']
            self.fileQuit()
            
        def on_click_stop(self):
            rr['inew'][self.workind]=rr['ii']
            rr['jnew'][self.workind]=rr['jj']
            rr['dsrcnew'][self.workind]=rr['dsrc']
            rr['qbardirnew'][self.workind]=rr['qbardir']
            rr['stop']=True
            self.fileQuit()


        def on_click_dir(self):
            rr['dsrcnew'][self.workind]=int(self.textbox.text())
            rr['qbardirnew'][self.workind]=int(self.textbox2.text())

            # Update plot 
            self.update_plot()

            print( 'dsrc/qbardir updated !')
      
        def on_mesh_click(self):
            # Add a new quadmesh to collections. Then copy it (so they do not have same id) at the first index of collection.
            if self.showmesh==0:
                self.showmesh=1
                self.dd.ax.collections[0]=copy.copy(self.dd.ax.pcolormesh(rr['lonr'],rr['latr'],rr['maskr'],edgecolor='black',lw=0.01,zorder=0))
            else:
                self.showmesh=0
                self.dd.ax.collections[0]=copy.copy(self.dd.ax.pcolormesh(rr['lonr'],rr['latr'],rr['maskr'],zorder=0))
            
            self.dd.ax.collections[-1].remove() # To finish we can delete the last collections quadmesh created (as it is now store in the fisrt index)
            self.dd.draw()

        def on_addriv_click(self):
            count = self.rivpoint.count()
            self.rivpoint.addItems([rr['rivername']+'.'+str(self.countriv)])
            self.rivpoint.setCurrentIndex(count)
            self.countriv+=1
            # update plot stuffs
            self.dd.locnbriv=self.rivpoint.count()
            rr['newname']=np.append(rr['newname'],rr['rivername']+'.'+str(self.countriv-1))
            rr['inew']=np.append(rr['inew'],np.nan)
            rr['jnew']=np.append(rr['jnew'],np.nan)
            rr['dsrcnew']=np.append(rr['dsrcnew'],int(self.textbox.text()))
            rr['qbardirnew']=np.append(rr['qbardirnew'],int(self.textbox2.text()))

            self.on_selectriv()
          

        def on_delriv_click(self):
            cindex = self.rivpoint.currentIndex()
            if self.rivpoint.count()==1:
                # faire un popup puis return
                msg=QMessageBox()            
                msg.setText("Need at least one point for this river!")
                x = msg.exec_()
                return
            self.rivpoint.removeItem(cindex)
            # update plot stuffs
            self.dd.locnbriv=self.rivpoint.count()
            rr['newname']=np.delete(rr['newname'],cindex)
            rr['inew']=np.delete(rr['inew'],cindex)
            rr['jnew']=np.delete(rr['jnew'],cindex)
            rr['dsrcnew']=np.delete(rr['dsrcnew'],cindex)
            rr['qbardirnew']=np.delete(rr['qbardirnew'],cindex)

            self.on_selectriv()


        def update_plot(self):
            
           # Del all points in order to reorganize the line2D object and put the current point at the end
            while len(self.dd.ax.lines)>2:
                self.dd.ax.lines[-1].remove() #remove point
                self.dd.ax.collections[-1].remove() # remove quiver
 
            for ii in range(self.rivpoint.count()+1):
                if ii == self.workind :#or np.isnan(rr['inew'][ii]):
                    continue
                if ii==self.rivpoint.count():
                   locind=self.workind
                   col='oc'
                else:
                   locind=ii
                   col='og' 
                if np.isnan(rr['inew'][locind]): # Point not defined yet
                    continue
   
                jz,iz=int(rr['jnew'][locind]),int(rr['inew'][locind])
                if rr['dsrcnew'][locind]==0:
                    self.dd.ax.plot(rr['lonu'][jz-1,iz-1],rr['latu'][jz-1,iz-1],col)
                    self.dd.ax.quiver(rr['lonu'][jz-1,iz-1],rr['latu'][jz-1,iz-1],rr['qbardirnew'][locind],0)
                elif rr['dsrcnew'][locind]==1:
                    self.dd.ax.plot(rr['lonv'][jz-1,iz-1],rr['latv'][jz-1,iz-1],col)
                    self.dd.ax.quiver(rr['lonv'][jz-1,iz-1],rr['latv'][jz-1,iz-1],0,rr['qbardirnew'][locind])
            self.dd.draw()

            
            
        
        def on_selectriv(self):
            self.workind=self.rivpoint.currentIndex()
            # update plot vars
            self.dd.locworkind=self.workind
            self.update_plot()
        

        def about(self):
            QtWidgets.QMessageBox.about(self, "How to use the GUI ?",
                "Preselected River position "+rr['rivername']+
                                " (green cross). First guess Runoff position (red dot).\n \n"+
                                " Right click to select new runoff position (cyan dot).\n\n"+
                                " Add/Delete a point for the current runoff with buttons on the right. The active point is colored cyan, the others are colored green.\n\n"+
                                " Yellow/mauve areas are corresponding to water/land mask."
                                    )
                                    
    # qApp = QtWidgets.QApplication(sys.argv)

    os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
#    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True) #enable highdpi scalin
#    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True) #use highdpi icons

    # To avoid seg fault when multiple rivers (cf :https://stackoverflow.com/questions/29451285/loading-a-pyqt-application-multiple-times-cause-segmentation-fault)
    global qApp
    qApp = QtWidgets.QApplication.instance()
    if qApp is None:
        qApp = QtWidgets.QApplication([rr])
    aw = ApplicationWindow() 
    aw.setWindowTitle("Riviere: %s" % rivername)
    aw.show()
    qApp.exec_()

    return rr['newname'],rr['inew'],rr['jnew'],rr['dsrcnew'],rr['qbardirnew'],rr['stop']



def correc_qtriver(river,grd=None):#, load_previous_corrections=False):
    # ---- Correction manuelle des j et i des rivieres ----
    # -------------------------------------------------------------------------------------
    # river : river dictionnary
    # 
    # Utilisation de l'interface graphique: qtriver

    print( 'MANUAL correction of runoff positions ...')


    # -- En cas de traitement partie, option pour pouvoir recharger les corrections manuelles precedentes si besoin:

#    if load_previous_corrections:
#        river=pickle.load(open('save_manual_correction.pkl','r'))

    lon_r = grd.lon
    lat_r = grd.lat
    lon_u = grd.lonu
    lat_u = grd.latu
    lon_v = grd.lonv
    lat_v = grd.latv
    mask_r = grd.maskr 
    mask_u = grd.umask # 0 masque / 1 non masque
    mask_v = grd.vmask # 0 masque / 1 non masque
    river_out=dict()
    for ll in river.keys():
        stop=False
        print( ' --- '+ll+' --- ')
        print( 'Before:',river[ll]['ii'],river[ll]['jj'])
        rivnamenew,inew,jnew,dsrcnew,qbardirnew,stop= qtriver(mask_r,mask_u, mask_v,lon_r,lat_r,lon_u, lon_v, lat_u, lat_v,
            river[ll]['longitude'], river[ll]['latitude'], ll, river[ll]['ii'], river[ll]['jj'],river[ll]['dsrc'],river[ll]['qbardir'])
        print( 'After:',inew,jnew,dsrcnew,qbardirnew)
        if len(inew)==1:
            river_out[ll]=dict()
            river_out[ll]['time'] = river[ll]['time']
            river_out[ll]['flow'] = river[ll]['flow']
            river_out[ll]['longitude']  = river[ll]['longitude']
            river_out[ll]['latitude']  = river[ll]['latitude']
            river_out[ll]['ii']=int(np.squeeze(inew))
            river_out[ll]['jj']=int(np.squeeze(jnew))
            river_out[ll]['dsrc']=int(np.squeeze(dsrcnew))
            river_out[ll]['qbardir']=float(np.squeeze(qbardirnew))
        else:
            print('River %s is being divided in %i points (and so is the flow)' %(ll,len(inew)))
            for rivnb,rivname in enumerate(rivnamenew):
                 river_out[rivname]=dict()
                 river_out[rivname]['time']=river[ll]['time']
                 river_out[rivname]['flow']=river[ll]['flow']/len(rivnamenew)
                 river_out[rivname]['longitude']=river[ll]['longitude']
                 river_out[rivname]['latitude']=river[ll]['latitude']
                 river_out[rivname]['ii']=int(np.squeeze(inew[rivnb]))
                 river_out[rivname]['jj']=int(np.squeeze(jnew[rivnb]))
                 river_out[rivname]['dsrc']=np.squeeze(dsrcnew[rivnb])
                 river_out[rivname]['qbardir']=np.squeeze(qbardirnew[rivnb])
        if stop:
            break

#    # --- Save manual modifications
#    pickle.dump(river, open('save_manual_correction.pkl','wb'))

    return river_out


