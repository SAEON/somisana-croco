
import copernicusmarine
from dateutil import rrule,relativedelta
from datetime import datetime
import os
ENV_CP = os.path.dirname(__file__)
ENV_MOD = os.path.join(ENV_CP, '../Modules')
import sys
sys.path.append(ENV_MOD)
import croco_class as Croco

'''
This script allows to download dataset from CMEMS website using the new
Copernicus Marine API. Further informations on the API is given here:
https://help.marine.copernicus.eu/en/articles/7949409-copernicus-marine-toolbox-introduction

This script uses "subset" function of copernicusmarine API. This function allows 
to download a subset of the original files.

!!!! WARNING !!!!
When using subset, dataset is converted to Analysis-Ready Cloud-Optimized 
(ARCO) format. Data are the same except for time variable where date is 
centered at midnight instead of noon for original files.
'''



#
Ystart,Mstart = 2013,1 # Starting month
Yend,Mend = 2013,3    # Ending month 

product_id = 'cmems_mod_glo_phy_my_0.083deg_P1M-m'
#cmems_mod_glo_phy_my_0.083deg_P1D-m --> Reana DAILY data
#cmems_mod_glo_phy_my_0.083deg_P1M-m --> Reana MONTHLY data
#cmems_mod_glo_phy-cur_anfc_0.083deg_P1M-m --> Currents Forecast MONTHLY data
#cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m --> Current Forecast DAILY data
# 
variables = ["zos","uo","vo","thetao","so"]
multi_files = False

croco_dir = '../'
croco_grd = 'croco_grd.nc'

###
output_dir = './GLORYS_DATA/'
output_prefix = 'glo12-reana-daily'

################################################################################
###########################  END OF USER DATA  #################################
################################################################################

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

grd = Croco.CROCO_grd(f"{croco_dir}/{croco_grd}")
area = [grd.latmax()+0.1, grd.lonmin()-0.1, grd.latmin()-0.1, grd.lonmax()+0.1]


date_start = datetime(Ystart,Mstart,1)
date_end = datetime(Yend,Mend,1)

print('Please specify you CMEMS login/password')
copernicusmarine.login()

for dt in rrule.rrule(rrule.MONTHLY, dtstart=date_start, until=date_end):
    month_beg = datetime.strftime(dt,'%Y-%m-%d %H:%M:%S')
    month_end = datetime.strftime(dt+
                            relativedelta.relativedelta(days=1),'%Y-%m-%d %H:%M:%S'
                                 )
    if multi_files:
        vars_loop = variables
    else:
        vars_loop = [variables]


    for vv in vars_loop:
        if multi_files:
            out_name = f"{output_prefix}_{vv}_{datetime.strftime(dt,'%Y_%m')}.nc"
            vars_down = [vv]
        else:
            out_name = f"{output_prefix}_{datetime.strftime(dt,'%Y_%m')}.nc"
            vars_down = vv

        copernicusmarine.subset(
                        dataset_id=product_id,
                        variables=vars_down,
                        minimum_longitude = area[1],
                        maximum_longitude = area[3],
                        minimum_latitude = area[2],
                        maximum_latitude = area[0],
                        start_datetime=month_beg,
                        end_datetime=month_end,
                        output_directory=output_dir,
                        output_filename=out_name,
                        force_download=True
                    )