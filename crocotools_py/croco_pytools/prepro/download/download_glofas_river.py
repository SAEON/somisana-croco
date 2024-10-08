import cdsapi
import datetime
import os
ENV_CP = os.path.dirname(__file__)
ENV_MOD = os.path.join(ENV_CP, '../Modules')
import sys
sys.path.append(ENV_MOD)
import croco_class as Croco

'''
download_glofas_river.py
========================

Download monthly grib files of daily river discharges from CMEMS-GLOFAS
https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-glofas-historical?tab=overview
Tranform each grib file into a netcdf file

Require ECMWF CDS API for data request
https://cds.climate.copernicus.eu/api-how-to#
Require CDO for grib to netcdf conversions

Once the files are dowloaded a mean file should be computed with ncra

Pierrick Penven 17/11/2022
'''


#
################################################################################
################################################################################
#
#
#
# Main Program
# 
# 
#
################################################################################
###########################  USER DATA  ########################################
################################################################################
#
Ystart,Mstart = 2013, 1   # Starting month
Yend,Mend  = 2013, 1      # Ending month 

product = 'cems-glofas-historical'

croco_dir = '../'
croco_grd = 'croco_grd.nc'

###
convert2netcdf = True
output_dir = './DATA_RIVER/'
output_name = 'cems_glofas'

#
################################################################################
###########################  END OF USER DATA  #################################
################################################################################
#

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

grd = Croco.CROCO_grd(f"{croco_dir}/{croco_grd}")
area = [grd.latmax()+0.1, grd.lonmin()-0.1, grd.latmin()-0.1, grd.lonmax()+0.1]

months=['january','february','march','april',\
        'may','june','july','august',\
	'september','october','november','december']

#
# Start the CDS API
#

c = cdsapi.Client()

#
# Loop on the years and months
#

for Y in range(Ystart,Yend+1):

  if Y==Ystart: 
    mo_min=Mstart
  else:
    mo_min=1

  if Y==Yend:
    mo_max=Mend
  else:
    mo_max=12

  for M in range(mo_min,mo_max+1):

    outfile_grib= f"{output_dir}/{output_name}_Y{Y}_M{M:02d}.grb"
    outfile_nc  = f"{output_dir}/{output_name}_Y{Y}_M{M:02d}.nc"

#
# Define the request
#
    options = {
        'variable': 'river_discharge_in_the_last_24_hours',
        'format': 'grib',
        'system_version': 'version_4_0',
        'hydrological_model': 'lisflood',
        'product_type': 'consolidated',
        'hyear': str(Y),
        'hmonth': months[(M-1)],
        'hday': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',  ],
        'format': 'grib',
        'area': area,}

#
# Printing message on screen
#

    info_time_clock = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print('                                                           ')
    print('-----------------------------------------------------------')
    print(''+info_time_clock)
    print(' GLOFAS data request, please wait...                       ')
    print('-----------------------------------------------------------')
    print('Request options: ')
    print(options)

#
# Send the request
#

    c.retrieve(product,options,outfile_grib)

#
# Convert grib files into Netcdf files
#
    if convert2netcdf:
        cmd = 'cdo -f nc copy ' + outfile_grib + ' ' + outfile_nc
        os.system(cmd)
#
# Done
#
if convert2netcdf:
    merge_cmd = \
        f"ncrcat {output_dir}/{output_name}*.nc {output_dir}/{output_name}.nc"
    os.system(cmd)
