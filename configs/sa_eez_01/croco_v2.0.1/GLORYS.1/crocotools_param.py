'''
params required by make_ini() and make_bry()
'''

# Input data information and formating
inputdata = 'mercator' # Input data dictionnary as defined in the Readers/ibc_reader.py

# default value to consider a z-level fine to be used
Nzgoodmin = 4 

# multi_files
multi_files = False

# format of monthly files
# (used in converting the datetime for the month to the file name
input_file_fmt = '%Y_%m.nc'

# tracers
tracers = ['temp','salt']

# relative path to the CROCO grid
croco_grd = '../GRID.1/croco_grd.nc.1'
sigma_params = dict(theta_s=7, theta_b=2, N=40, hc=200) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth

# Ini filename prefix
ini_prefix = 'croco_ini_GLORYS'
ini_suffix = '.1'

# Bry filename prefinformations
bry_prefix = 'croco_bry_GLORYS' 
obc_dict = dict(south=1, west=1, east=1, north=0) # open boundaries (1=open , [S W E N])
cycle_bry=0


