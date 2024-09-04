__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.le.corre@shom.fr'
__date__   = '2022-09'
__license__='GPL3'
'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO initial and bry file for zoom

===========================================================================
'''

#--- Dependencies ---------------------------------------------------------
import glob
import numpy as np
import sys
sys.path.append("./Modules/")
import os
import toolsf

#--- USER CHANGES ---------------------------------------------------------

prt_grd = '../../CROCO_FILES/croco_grd.nc'         # Parent grid file
chd_grd = '../../CROCO_FILES/croco_chd_grd.nc'     # Child grid file
ts = 6; tb = 4; hc = 75; n = 50  # Child vertical coordinates parameters

tracers=['temp','salt'] # if no tracers leave empty. string must be 4 or less

### Where to find and put croco files
croco_dir = '../../CROCO_FILES/'

### Ini file ###
make_ini = True  # Do you want to build ini file
prt_his_ini = '../../CROCO_FILES/croco_his_20050201_20050205.nc' # History file to start child simulation
rec = 1          # record index in the ini file

### Bry file ###
make_bry = False  # Do you want to build bry file
obc_cond = 'SWEN' # First letters of the boundaries that are opened (S=South, W=West, E=East, N=North)
# list of all the files you desired to build you bry file
#    can handle all bash sign such as ?,*,[]
#    the only constraint is that each string may be less than 160 characters
prt_his_bry = ['../../CROCO_FILES/croco_his_2005021*','../../CROCO_FILES/croco_his_2005020?*']

#--- END USER CHANGES -----------------------------------------------------

if len(tracers)>0:
    all_tracers = np.zeros((20,len(tracers) ), dtype='c')
    for i in range(len(tracers)):
        len_trac=len(tracers[i])
        if len_trac>20:
            print('Tracer character string must be smaller than 20.\nError with %s' %tracers[i])
        all_tracers[0:len_trac,i]=tracers[i]
        all_tracers[len_trac:,i]=' '
else:
    all_tracers=[]

# --- Make ini ------------------------------------------------------------
if make_ini:
    toolsf.r2r_init(chd_grd,ts,tb,hc,n,prt_grd,prt_his_ini,rec,all_tracers)

    # move files into croco_dir
    os.rename('croco_chd_ini.nc', croco_dir + 'croco_chd_ini.nc')
    os.rename('croco_init_diag.nc', croco_dir + 'croco_init_diag.nc')

# --- Make bry ------------------------------------------------------------
if make_bry:

    # --- Build input file list from prt_his_bry --------------------------

    inputfiles=[]
    for i in range(len(prt_his_bry)):
        if len(prt_his_bry[i])>160:
            print('The size of his path need to be less than 160. It pos %i it is %i' % (i,len(prt_his_bry[i])))
            sys.exit()
        inputfiles+=glob.glob(prt_his_bry[i])
    
    inputfiles=list(set(inputfiles))

    all_files = np.zeros((160,len(inputfiles) ), dtype='c')

    for i in range(len(inputfiles)):
        ll=len(inputfiles[i])
        all_files[0:ll,i]=inputfiles[i]

    # --- Create child bry ------------------------------------------------
    toolsf.r2r_bry(chd_grd,ts,tb,hc,n,obc_cond,prt_grd,all_files,all_tracers)

    # move files into croco_dir
    os.rename('croco_chd_bry.nc', croco_dir + 'croco_chd_bry.nc')
