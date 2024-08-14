__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.le.corre@shom.fr'
__date__   = '2022-09'
__license__='GPL3'

'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO bounday files
In the current state the script can handle:
    - mercator (glorys)
    - soda
    - eccov4  (V4.3)

To add a new dataset you just have to go in Readers/inputs_readers.py and
create a dico with correspondance between croco_var(ssh,u,v,temp,salt)
and dataset variable names.

The periodicity of a dataset (if it is -180/180 or 0/360) is handle by the
script. Just download the dataset and it should be fine

The script works as follow:
    - reads croco_grd
    - reads input data and restricts the area of spatial coverage
    - creates croco_bry.nc
    - Loop on open boundaries with:
          check if there is data before or after for continuity, if not duplicate first or last
          loop on var with:
              * horizontal interpolation
              * vertical interpolation
    - Writes data in netcdf

===========================================================================
'''

#--- Dependencies ---------------------------------------------------------

import netCDF4 as netcdf
import xarray as xr
import pylab as plt
import numpy as np
import glob as glob
from dateutil.relativedelta import relativedelta
import sys
sys.path.append("./Modules/")
sys.path.append("./Readers/")
import interp_tools
import sigmagrid_tools as sig_tools
import Cgrid_transformation_tools as grd_tools
import croco_class as Croco
import ibc_class as Inp

#--- USER CHANGES ---------------------------------------------------------

# Dates
Yorig = 2013                    # year defining the origin of time as: days since Yorig-01-01
Ystart, Mstart = '2013', '01'   # Starting month
Yend, Mend  = '2013','03'       # Ending month 

# Input data information and formating
inputdata = 'mercator_croco'    # Input data dictionnary as defined in the Readers/ibc_reader.py
input_dir = '../../MERCATOR_GLOB_2013/'
input_prefix = 'mercator_*'  # Please use * to include all files
multi_files = False
if multi_files: # Multiple data files. Time is read in ssh file
    input_file = {'ssh':sorted(glob.glob(input_dir+input_prefix+'ETAN.*.nc')),\
                  'temp':sorted(glob.glob(input_dir+input_prefix+'THETA.*.nc')),\
                  'salt':sorted(glob.glob(input_dir+input_prefix+'SALT.*.nc')),\
                  'u':sorted(glob.glob(input_dir+input_prefix+'EVEL.*.nc')),\
                  'v':sorted(glob.glob(input_dir+input_prefix+'NVEL.*.nc'))\
                }
else:  # glob all files
    input_file  = sorted(glob.glob(input_dir + input_prefix))

# default value to consider a z-level fine to be used
Nzgoodmin = 4

# Tracers
tracers = ['temp', 'salt']

# CROCO grid informations
croco_dir = '../../CROCO_FILES/' 
croco_grd = 'croco_grd.nc'
sigma_params = dict(theta_s=7, theta_b=2, N=32, hc=200) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth

# Bry file informations
bry_filename = 'croco_bry.nc' # output will be put in croco_dir by default
obc_dict = dict(south=1, west=1, east=1, north=1) # open boundaries (1=open , [S W E N])
output_file_format = "MONTHLY" # How outputs are spit (MONTHLY,YEARLY,FULL)
cycle_bry = 0.

#--- END USER CHANGES -----------------------------------------------------


if __name__ == '__main__':

    # Put origin date to the right format
    day_zero   = str(Yorig)+'0101'    
    day_zero_num = plt.datetime.datetime(int(day_zero[:4]),
                                         int(day_zero[4:6]),
                                         int(day_zero[6:8]))
    day_zero_num = plt.date2num(day_zero_num)

    # Put start and end date to the right format
    start_date = Ystart+Mstart+'01'+'12'  # defaut start day is 1st
   

    dtstrdt = plt.datetime.datetime(int(start_date[:4]),
                                    int(start_date[4:6]),
                                    int(start_date[6:8]),
                                    int(start_date[8:]))

    dtenddt = plt.datetime.datetime(int(Yend),int(Mend),1,12) \
            + relativedelta(months=1,days=-1) # Last day of the ending month


    dtstr, dtend = plt.date2num(dtstrdt), plt.date2num(dtenddt)

    # --- Load croco_grd --------------------------------------------------

    crocogrd = Croco.CROCO_grd(''.join((croco_dir, croco_grd)), sigma_params)

    # --- Initialize boundary vars ----------------------------------------

    crocogrd.WEST_grid()
    crocogrd.EAST_grid()
    crocogrd.SOUTH_grid()
    crocogrd.NORTH_grid()
    
    # --- Initialize input data class -------------------------------------

    inpdat = Inp.getdata(inputdata,input_file,crocogrd,multi_files,
                         tracers,
                         bdy=[obc_dict,cycle_bry])

    # --- Work on date format for the loop in time ------------------------

    startloc=plt.datetime.datetime(int(start_date[:4]),
                                   int(start_date[4:6]),
                                   1)
    if output_file_format.upper() == "MONTHLY":
        endloc= startloc+relativedelta(months=1,days=-1,hours=12)
    elif output_file_format.upper() == "YEARLY":
        if plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S").year == int(Yend) :
            endloc=plt.num2date(dtend).replace(tzinfo=None)
        else:
            endloc= plt.datetime.datetime(int(start_date[:4]), 12,31,12)

    elif output_file_format.upper() == "FULL":
        endloc=plt.num2date(dtend).replace(tzinfo=None)
    else:
        print("\n Output file format \"%s\" is not setup. Pease change it to MONTHLY, YEARLY or FULL")
        sys.exit()
 
    # --- Start time loop loop in time ------------------------------------

    while plt.date2num(endloc) <= dtend:

        # Load full time dataset
        time = plt.date2num(inpdat.ncglo['time'].values)
        # find index for the time range 
        ind= np.where((time>plt.date2num(startloc)) & (time<=plt.date2num(endloc)))
 
        if len(ind[0])==0 :
            print('\nData is missing for range %s to %s' % (startloc ,endloc))
            sys.exit()
            
        [dtmin,dtmax]=np.min(ind),np.max(ind)
        # create monthly file
        tmp_date = plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S")
            # file name depending on format chosen
        if output_file_format.upper() == "MONTHLY":
            bdy_filename = croco_dir+bry_filename.replace('.nc', '_%s_Y%sM%02i.nc' %(inputdata,tmp_date.year,tmp_date.month))
        elif output_file_format.upper() == "YEARLY":
            bdy_filename = croco_dir+bry_filename.replace('.nc', '_%s_Y%s.nc' %(inputdata,tmp_date.year))
        elif output_file_format.upper() == "FULL":
            bdy_filename = croco_dir+bry_filename.replace('.nc', '_%s.nc' %(inputdata))

        Croco.CROCO.create_bry_nc(None,bdy_filename,crocogrd,obc_dict,cycle_bry,tracers=tracers)
        #
        print('\n-----------------------------------')
        if output_file_format.upper() == "MONTHLY":
            print(' Processing Year %s - Month %02i' %(tmp_date.year,tmp_date.month))
        elif output_file_format.upper() == "YEARLY":
            print(' Processing Year %s' %(tmp_date.year))
        elif output_file_format.upper() == "FULL": 
            tmp_end_date = plt.datetime.datetime.strptime(str(endloc), "%Y-%m-%d %H:%M:%S")
            print(' Processing from Year %s - Month %02i  to Year %s - Month %02i' %(tmp_date.year,tmp_date.month,tmp_end_date.year,tmp_end_date.month))
            del tmp_end_date

        # Check if there is at least 1 point by month when doing Yearly or full
        if output_file_format.upper() != "MONTHLY":
            tmp_str_date=startloc
            tmp_end_date=tmp_str_date+relativedelta(months=1,days=-1,hours=12)
            while(plt.date2num(tmp_end_date) <= plt.date2num(endloc)):
                ind_tmp= np.where((time>plt.date2num(tmp_str_date)) & (time<=plt.date2num(tmp_end_date)))
                tmp_date=plt.datetime.datetime.strptime(str(tmp_str_date), "%Y-%m-%d %H:%M:%S")
                if len(ind_tmp[0])==0:
                    print('\nLacking %s data for  Y%s - M%02i' %(inputdata,tmp_date.year,tmp_date.month))
                    sys.exit()            
                tmp_str_date=tmp_end_date+relativedelta(days=1,hours=-12)
                tmp_end_date=tmp_str_date+relativedelta(months=1,days=-1,hours=12)

            del tmp_str_date,tmp_end_date

        # --- Check if data availabla for the surrounded months -----------
        print('-----------------------------------')
        tmp_str_date=startloc ; tmp_end_date = endloc
        prev_month_str = tmp_str_date+relativedelta(months=-1)
        prev_month_end = tmp_str_date+relativedelta(days=-1,hours=12)
        ind_prev       = np.where((time>plt.date2num(prev_month_str)) & (time<plt.date2num(prev_month_end)))
        if len(ind_prev[0])==0:
            print('   No data for the previous month: using current month')
            prev=1
        else:
            prev=0
            dtmin=dtmin-1 # create overlap before (in this case it takes the previous month)
        del prev_month_str,prev_month_end,ind_prev
        #
        next_month_str = tmp_end_date+relativedelta(days=1)
        next_month_end = next_month_str+relativedelta(months=1,days=-1,hours=12)
        ind_next       = np.where((time>plt.date2num(next_month_str)) & (time<plt.date2num(next_month_end)))
        if len(ind_next[0])==0:
            print('   No data for the next month: using current month')
            nxt=1
        else:
            nxt=0
            dtmax=dtmax+1 # create overlap after (in this case it takes the next month)

        del next_month_str,next_month_end,ind_next,tmp_str_date,tmp_end_date
        
        if np.nanvar(np.gradient(time[dtmin:dtmax+1])) >=5: # Abnormal distribution of days
            Question = input( "Abnormal distribution of days (variance to high) \
                    \nThis may be due to the use of different temproral resolution dataset.\
                    \n Do you want to proceed?: y,[n] ") or 'no'
            if Question.lower() == ("n") or Question.lower() == ("no"):
                print('Aborting')
                sys.exit()
         
        ## --- Handle bry_time --------------------------------------------

        bry_time= time[dtmin:dtmax+1] - day_zero_num 
        if prev == 1 and len(bry_time)==1:
            prev_time = plt.date2num(plt.num2date(bry_time[0]) + relativedelta(days=-30))
            bry_time=np.append(prev_time,bry_time)
            del prev_time
        elif prev == 1 and len(bry_time)>1: 
            date_dt=np.gradient(bry_time)[0]
            prev_time = plt.date2num(plt.num2date(bry_time[0]) + relativedelta(days=-date_dt))
            bry_time=np.append(prev_time,bry_time)
            del prev_time

        if nxt == 1 and len(bry_time)==1:
            nxt_time = plt.date2num(plt.num2date(bry_time[-1]) + relativedelta(days=30))
            bry_time=np.append(bry_time,nxt_time)
            del nxt_time
        elif nxt == 1 and len(bry_time)>1:
            date_dt=np.gradient(bry_time)[-1]
            nxt_time = plt.date2num(plt.num2date(bry_time[-1]) + relativedelta(days=date_dt))
            bry_time=np.append(bry_time,nxt_time)
            del nxt_time
        
        nc=netcdf.Dataset(bdy_filename, 'a')

        nc.Input_data_type=inputdata 
        nc.variables['bry_time'].cycle=cycle_bry
        nc.variables['bry_time'][:]=bry_time
        if cycle_bry==0:
            nc.variables['bry_time'].units='days since %s-01-01 00:00:00' %(Yorig)
        # --- Loop on boundaries ------------------------------------------
  
        if len(tracers) == 0:
            var_loop = ['ssh','velocity']
        else:
            var_loop = ['ssh','tracers','velocity']

        for boundary, is_open in zip(obc_dict.keys(), obc_dict.values()):
            if is_open:
                for vars in var_loop:
                    print('\n     Processing *%s* for %sern boundary' %(vars, boundary))
                    print('     ------------------------------------------')
                    if vars == 'ssh': 
                        (zeta,NzGood) = interp_tools.interp_tracers(inpdat,vars,-1,crocogrd,dtmin,dtmax,prev,nxt,boundary[0].upper())
                        z_rho = crocogrd.scoord2z_r(zeta=zeta,bdy="_"+boundary)
                        z_w   = crocogrd.scoord2z_w(zeta=zeta,bdy="_"+boundary)
    
                    elif vars == 'tracers':
                        trac_dict = dict()                        
                        for trc in tracers:
                            print(f'\nIn tracers processing {trc}')
                            trac_dict[trc] = interp_tools.interp(inpdat,trc,Nzgoodmin,z_rho,crocogrd,dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())
        
                    elif vars == 'velocity':

                        cosa=np.cos(eval(''.join(('crocogrd.angle_',boundary))) )
                        sina=np.sin(eval(''.join(('crocogrd.angle_',boundary))) )

                        [u,v,ubar,vbar]=interp_tools.interp_uv(inpdat,Nzgoodmin,z_rho,cosa,sina,crocogrd,dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())

                        conserv=1  # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport          
                        if conserv == 1:
                            ubar_croco=sig_tools.vintegr4D(u,grd_tools.rho2u(z_w),grd_tools.rho2u(z_rho),np.nan,np.nan)[0]/grd_tools.rho2u(eval(''.join(('crocogrd.h_'+boundary))))
                            vbar_croco=sig_tools.vintegr4D(v,grd_tools.rho2v(z_w),grd_tools.rho2v(z_rho),np.nan,np.nan)[0]/grd_tools.rho2v(eval(''.join(('crocogrd.h_'+boundary))))
                            
                            u = u - np.tile(ubar_croco[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                            u = u + np.tile(ubar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))

                            v = v - np.tile(vbar_croco[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                            v = v + np.tile(vbar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))


    # --- Saving in netcdf ------------------------------------------------

                print('\nSaving %sern boundary in Netcdf' % boundary)
                print('----------------------------------')
  
                 # handle indices (as 2 points where taken next to bdy)
                if str(boundary) == 'west' and is_open:
                    indices3D="[:,:,:,0]" # T,N,J,i=0
                    indices2D="[:,:,0]"   # T,J,i=0
                elif str(boundary) == 'east' and is_open:
                    indices3D="[:,:,:,-1]" # T,N,J,i=last
                    indices2D="[:,:,-1]"   # T,J,i=last
                elif str(boundary) == 'south' and is_open:
                    indices3D="[:,:,0,:]" # T,N,j=0,I
                    indices2D="[:,0,:]"   # T,j=0,I
                elif str(boundary) == 'north' and is_open:
                    indices3D="[:,:,-1,:]" # T,N,j=last,I
                    indices2D="[:,-1,:]"   # T,j=last,I

                mask_zet = np.tile(eval(''.join(('crocogrd.maskr_',boundary))),[zeta.shape[0],1,1])
                if "velocity" in var_loop:
                    mask_u   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],u.shape[1],1,1]) 
                    mask_v   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[u.shape[0],u.shape[1],1,1])
                    mask_ubar   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],1,1])
                    mask_vbar   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[v.shape[0],1,1])
                
                nc.variables['zeta_'+str(boundary)][:]=eval(''.join(('zeta',indices2D)))*eval(''.join(('mask_zet',indices2D)))
                nc.variables['u_'+str(boundary)][:]   =eval(''.join(('u',indices3D)))*eval(''.join(('mask_u',indices3D)))
                nc.variables['v_'+str(boundary)][:]   =eval(''.join(('v',indices3D)))*eval(''.join(('mask_v',indices3D)))
                nc.variables['ubar_'+str(boundary)][:]=eval(''.join(('ubar',indices2D)))*eval(''.join(('mask_ubar',indices2D)))
                nc.variables['vbar_'+str(boundary)][:]=eval(''.join(('vbar',indices2D)))*eval(''.join(('mask_vbar',indices2D)))

                if 'tracers' in var_loop:
                    for varname, value in zip(trac_dict.keys(), trac_dict.values()):
                        mask_tra = np.tile(eval(''.join(('crocogrd.maskr_',boundary))),[value.shape[0],value.shape[1],1,1])
                        nc.variables[f"{varname}_{boundary}"][:] = eval(f'value{indices3D}')*eval(''.join(('mask_tra',indices3D)))

                # handle prev and nxt + save

        nc.close()
    # --- END writting netcdf ---------------------------------------------
    
    # --- Preparing time for next loop ------------------------------------
        startloc=endloc+relativedelta(days=1,hours=-12)
        if output_file_format.upper() == "MONTHLY":
            endloc= startloc+relativedelta(months=1,days=-1,hour=12)
        elif output_file_format.upper() == "YEARLY":
            yearloc=plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S")
            if yearloc.year == int(Yend) :
                endloc=plt.num2date(dtend).replace(tzinfo=None) 
            else:
                endloc= plt.datetime.datetime(int(yearloc.year), 12,31,12)
        elif output_file_format.upper() == "FULL":
            endloc=startloc

