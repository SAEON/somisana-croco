#
# For ERA5 python crocotools parameters list
#
# CAUTION IT MUST BE CONSISTENT with your MATLAB CROCOTOOLS_PARAM.m file in Run directory
# *******************************************************************************
#                         U S E R  *  O P T I O N S
# *******************************************************************************
#
# Original ERA5 directory
#
era5_dir_raw = './eez' # this is where the output of ERA5_request.py goes
#
# Output ERA5 directory
#
era5_dir_processed = './eez_for_croco' # this is where the output of ERA5_convert.py goes 
#
# should we extract wave variables?
#
wave_extract=False # True to extract wave variables
#
# should we extract sea level pressure?
#
pressure_extract=True
#
# Dates limits
#
year_start = 2007
month_start = 12
year_end = 2014
month_end = 1
#
# Year origin of time
#
Yorig=1993
#
# Overlapping days (at the beginning/end of each month)
#
n_overlap = 0
#
# Request time (daily hours '00/01/.../23')
#
times = ['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00','08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00']
#
# Request area ([north, west, south, east])
#
lonmin=11
lonmax=36
latmin=-39
latmax=-25

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
variables = ['lsm'  , 'sst' , 'tp'        ,'strd'   ,'ssr'     ,'t2m'  ,'q'      ,'u10'  ,'v10']
conv_cff  = [1.     ,  1.   ,  cff_tp     ,cff_heat ,cff_heat  ,1.     ,1.       ,1.     ,1.] 
units     = ['(0-1)',  'K'  , 'kg m-2 s-1','W m-2'  ,'W m-2'   ,'K'    ,'kg kg-1','m s-1','m s-1']

if wave_extract:
    ## append waves variables
    wave_var=['swh', 'mwd', 'pp1d' ,'cdww'];variables.extend(wave_var)
    wave_conv_cff=[1.,  1., 1. , 1.];  conv_cff.extend(wave_conv_cff)
    wave_units=['m','Degrees true','s', 'dimensionless']; units.extend(wave_units)
    
if pressure_extract:
    ## append the pressure variable
    msl_var=['msl'];variables.extend(msl_var)
    msl_conv_cff=[1.];  conv_cff.extend(msl_conv_cff)
    msl_units=['Pa']; units.extend(msl_units)


# *******************************************************************************
#                         E N D     U S E R  *  O P T I O N S
# *******************************************************************************
