'''
params required by make_ini() and make_bry()
'''

# Input data information and formating
inputdata = 'mercator' # Input data dictionnary as defined in the Readers/ibc_reader.py

# default value to consider a z-level fine to be used
Nzgoodmin = 4 

# multi_files
multi_files = False

# tracers
tracers = ['temp','salt']

# relative path to the CROCO grid
croco_grd = '../GRID/croco_grd.nc'
sigma_params = dict(theta_s=5, theta_b=7, N=30, hc=200) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth

# Ini filename prefix
ini_prefix = 'croco_ini_MERCATOR'

# Bry filename prefinformations
bry_prefix = 'croco_bry_MERCATOR' 
obc_dict = dict(south=1, west=1, east=0, north=1) # open boundaries (1=open , [S W E N])
cycle_bry=0


