# params required by make_ini() and make_bry()

# Dates
# starting date
Yorig, Morig, Dorig = '2000','01','01' # Month and days need to be 2-digits format

# Input data information and formating
inputdata = 'mercator' # Input data dictionnary as defined in the Readers/ibc_reader.py

# default value to consider a z-level fine to be used
Nzgoodmin = 4 

# multi_files
multi_files = False

# tracers
tracers = ['temp','salt']

# CROCO grid information
croco_dir = '/home/g.rautenbach/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/GRID/'
croco_grd = croco_dir+'croco_grd.nc'
sigma_params = dict(theta_s=5, theta_b=7, N=30, hc=200) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth

# Ini file information
ini_filename = 'croco_ini.nc'

# Bry file informations
bry_filename = 'croco_bry.nc' # output will be put in croco_dir by default
obc_dict = dict(south=1, west=1, east=0, north=1) # open boundaries (1=open , [S W E N])
cycle_bry=0
