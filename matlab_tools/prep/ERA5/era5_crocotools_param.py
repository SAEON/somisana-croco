#
# For ERA5 python crocotools parameters list
#
# CAUTION IT MUST BE CONSISTENT with your MATLAB CROCOTOOLS_PARAM.m file in Run directory
# *******************************************************************************
#                         U S E R  *  O P T I O N S
# *******************************************************************************
#
# General path
#
##config_dir = '../croco/Run_TEST/'           # must be the same than crocotools_param
config_dir = 'path_to_my_run_dir/'            # must be the same than crocotools_param
config_name = 'Benguela_LR'
#
# Original ERA5 directory
#
era5_dir_raw = config_dir + 'DATA/ERA5_native_' + config_name
#
# Output ERA5 directory
#
era5_dir_processed = config_dir + 'DATA/ERA5_' + config_name
#
#
# Dates limits
#
year_start = 2005
month_start = 1
year_end = 2005
month_end = 3
#
# Year origin of time
#
Yorig=2000
#
# Overlapping days (at the beginning/end of each month)
#
n_overlap = 0
#
# Request time (daily hours '00/01/.../23')
#
time = '00/01/02/03/04/05/06/07/08/09/10/11/12/13/14/15/16/17/18/19/20/21/22/23'
#
# Request variables (see available at ERA5_variables.json)
variables = ['lsm','tp','strd','ssr','t2m','q','u10','v10'] #note lsm is land_sea_mask
#
# Request area ([north, west, south, east])
#
ownArea = 0 	# 0 if area from a crocotools_param.m file
                # 1 if own area

if ownArea == 0: 
    # To complete if ownArea==0
    paramFile = config_dir + 'crocotools_param.m' # path the crocotools_param file of the simulation
    
else:
    # To complete if ownArea==1
    lonmin=7
    lonmax=23
    latmin=-45
    latmax=-20
#
# Variable names and conversion coefficients  
# TP: convert from accumlated m in a hour into   kg m-2 s-1
#
cff_tp=1000./3600. # m in 1 hour -> kg m-2 s-1
# Heat flux J m-2 in one hour into W m-2
#
cff_heat=1./3600.   # J m-2 in 1 hour -> W m-2
# Names, conversion coefficients and new units
#
variables = ['lsm'  , 'sst' , 'tp'        ,'strd'   ,'ssr'     ,'t2m'  ,'q'      ,'u10'  ,'v10'  ]
conv_cff  = [1.     ,  1.   ,  cff_tp     ,cff_heat ,cff_heat  ,1.     ,1.       ,1.     ,1.     ] 
units     = ['(0-1)',  'K'  , 'kg m-2 s-1','W m-2'  ,'W m-2'   ,'K'    ,'kg kg-1','m s-1','m s-1']



# *******************************************************************************
#                         E N D     U S E R  *  O P T I O N S
# *******************************************************************************
