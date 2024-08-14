__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.le.corre@shom.fr'
__date__   = '2022-09'
__license__='GPL3'
'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO initial file
In the current state the script can handle:
    - mercator (glorys)
    - soda
    - eccov4  (V4.3)

To add a new dataset you just have to go in Readers/ibc_readers.py and
create a dico with correspondance between croco_var(ssh,u,v,temp,salt)
and dataset variable names.

The periodicity of a dataset (if it is -180/180 or 0/360) is handle by the
script. Just download the dataset and it should be fine

The script works as follow:
    - reads croco_grd
    - reads input data and restricts the area of spatial coverage
    - creates croco_ini.nc
    - Loop on var with:
        * horizontal interpolation
        * vertical interpolation
    - Writes data in netcdf
===========================================================================
'''

#--- Dependencies ---------------------------------------------------------

import netCDF4 as netcdf
import pylab as plt
import numpy as np
import glob as glob
from datetime import datetime
import sys
sys.path.append("./Modules/")
sys.path.append("./Readers/")
import Cgrid_transformation_tools as grd_tools
import interp_tools
import sigmagrid_tools as sig_tools
import croco_class as Croco
import ibc_class as Inp

#--- USER CHANGES ---------------------------------------------------------

# Dates
# starting date
Yini, Mini, Dini  = '2013','01','01' # Month and days need to be 2-digits format
# reference time (default = ini time)
Yorig, Morig, Dorig = Yini, Mini, Dini # Month and days need to be 2-digits format

# Input data information and formating
inputdata = 'mercator_croco' # Input data dictionnary as defined in the Readers/ibc_reader.py
input_dir = '../../MERCATOR_GLOB_2013/'
input_prefix='mercator_'

input_file  = f'{input_dir}{input_prefix}Y2013M1.cdf'
multi_files=False # If variables are in different netcdf
if multi_files: # Mutiple files
    input_file = { 'ssh'  : input_dir + input_prefix + 'ETAN.Y2013M01.nc',\
                   'temp' : input_dir + input_prefix + 'THETA.Y2013M01.nc',\
                   'salt' : input_dir + input_prefix + 'SALT.Y2013M01.nc',\
                   'u'    : input_dir + input_prefix + 'EVEL.Y2013M01.nc',\
                   'v'    : input_dir + input_prefix + 'NVEL.Y2013M01.nc'\
                }

# time index to use in the file
tndx = 0

# default value to consider a z-level fine to be used
Nzgoodmin = 4 

# tracers
tracers = ['temp','salt']

# CROCO grid informations
croco_dir = '../../CROCO_FILES/'
croco_grd = 'croco_grd.nc'
sigma_params = dict(theta_s=7, theta_b=2, N=32, hc=200) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth

# Ini file informations
ini_filename = 'croco_ini.nc' # output will be put in croco_dir by default

#--- END USER CHANGES -----------------------------------------------------

if __name__ == '__main__':
    # edit ini_filename to add starting date
    ini_filename = ini_filename.replace('.nc', '_%s_Y%sM%s.nc' %(inputdata,Yini, Mini))

    # Load croco_grd
    crocogrd = Croco.CROCO_grd(''.join((croco_dir, croco_grd)), sigma_params)

    # --- Load input (restricted to croco_grd) ----------------------------

    inpdat=Inp.getdata(inputdata,input_file,crocogrd,multi_files,tracers)

    print(' ')
    print(' Making initial file: '+ini_filename)
    print(' ')

    # --- Create the initial file -----------------------------------------

    Croco.CROCO.create_ini_nc(None,''.join((croco_dir + ini_filename)),crocogrd,
                              tracers=tracers)

    # --- Handle initial time ---------------------------------------------
    ini_date_num = datetime(int(Yini), int(Mini), int(Dini))
    ini_date_num = plt.date2num(ini_date_num) + 0.5

    day_zero_num = datetime(int(Yorig), int(Morig), int(Dorig))
    day_zero_num = plt.date2num(day_zero_num)
    
    tstart=0
    
    if ini_date_num != day_zero_num:
        tstart = ini_date_num - day_zero_num # days

    scrumt = tstart*3600*24 # convert in second
    oceant = tstart*3600*24
    tend=0.

   #  --- Compute and save variables on CROCO grid ---------------

    for vars in ['ssh','tracers','velocity']:
        print('\nProcessing *%s*' %vars)
        nc=netcdf.Dataset(croco_dir+ini_filename, 'a')
        if vars == 'ssh' :
            (zeta,NzGood) = interp_tools.interp_tracers(inpdat,vars,-1,crocogrd,tndx,tndx)
            nc.variables['zeta'][0,:,:] = zeta*crocogrd.maskr
            nc.Input_data_type=inputdata
            nc.variables['ocean_time'][:] = oceant
            nc.variables['scrum_time'][:] = scrumt
            nc.variables['scrum_time'].units='seconds since %s-%s-%s 00:00:00' %(Yorig,Morig,Dorig)
            nc.variables['tstart'][:] = tstart
            nc.variables['tend'][:] = tend

            z_rho = crocogrd.scoord2z_r(zeta=zeta)
            z_w   = crocogrd.scoord2z_w(zeta=zeta)
            
        elif vars == 'tracers':
            for tra in tracers:
                print(f'\nIn tracers processing {tra}')
                trac_3d= interp_tools.interp(inpdat,tra,Nzgoodmin,z_rho,crocogrd,tndx,tndx)
                nc.variables[tra][0,:,:,:] = trac_3d*crocogrd.mask3d()

        elif vars == 'velocity':

            cosa=np.cos(crocogrd.angle)
            sina=np.sin(crocogrd.angle)

            [u,v,ubar,vbar]=interp_tools.interp_uv(inpdat,Nzgoodmin,z_rho,cosa,sina,crocogrd,tndx,tndx)
              
            conserv=1  # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport          
            if conserv == 1:
                (ubar_croco,h0)=sig_tools.vintegr(u,grd_tools.rho2u(z_w),grd_tools.rho2u(z_rho),np.nan,np.nan)/grd_tools.rho2u(crocogrd.h)
                (vbar_croco,h0)=sig_tools.vintegr(v,grd_tools.rho2v(z_w),grd_tools.rho2v(z_rho),np.nan,np.nan)/grd_tools.rho2v(crocogrd.h)

                u = u - ubar_croco ; u = u + np.tile(ubar,(z_rho.shape[0],1,1))
                v = v - vbar_croco ; v = v + np.tile(vbar,(z_rho.shape[0],1,1))
           
            nc.variables['u'][0,:,:,:] = u *crocogrd.umask3d()
            nc.variables['v'][0,:,:,:] = v * crocogrd.vmask3d()
            nc.variables['ubar'][0,:,:] = ubar *crocogrd.umask
            nc.variables['vbar'][0,:,:] = vbar * crocogrd.vmask

   
    nc.close()



