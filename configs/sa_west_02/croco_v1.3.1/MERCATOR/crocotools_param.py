# params required by make_ini() and make_bry()

# reference time
Yorig, Morig, Dorig = '2000', '01', '01' # Month and days need to be 2-digits format

# Input data information and formating
inputdata = 'mercator' # Input data dictionnary as defined in the Readers/ibc_reader.py
input_dir = '/home/gfearon/tmp/test_forecast_run/20240417_06/downloaded_data/'

input_file  = f'{input_dir}mercator_20240417.nc'

# default value to consider a z-level fine to be used
Nzgoodmin = 4 

# tracers
tracers = ['temp','salt']

# CROCO grid information
croco_grd = '/home/gfearon/code/somisana-croco/configs/swcape_02/croco_v1.3.1/GRID/croco_grd.nc'
sigma_params = dict(theta_s=5, theta_b=7, N=30, hc=200) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth

# Ini file information
ini_filename = 'croco_ini.nc'
