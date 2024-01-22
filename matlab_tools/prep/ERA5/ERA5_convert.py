#!/usr/bin/env python
#
# ERA5_convert.py
#
# Script to convert ERA5 original files obtained from the Climate Data
# Store (CDS) of Copernicus https://cds.climate.copernicus.eu into
# a format and using units which can be used by the online interpolation of CROCO
#
#
# -------------------------------------------------
# Getting libraries and utilities
# -------------------------------------------------

from datetime import date
import json
from netCDF4 import Dataset as netcdf
import numpy as np
import os

# -------------------------------------------------
# Import my crocotools_param_python file
from era5_crocotools_param import *
# -------------------------------------------------

# -------------------------------------------------
# Setting processed output directory
# -------------------------------------------------
# Get the current directory
os.makedirs(era5_dir_processed,exist_ok=True)

# -------------------------------------------------
# Loading ERA5 variables's information as 
# python Dictionary from JSON file
# -------------------------------------------------

with open('ERA5_variables.json', 'r') as jf:
    era5 = json.load(jf)

#
# -------------------------------------------------
# Loop on Years and Months
# -------------------------------------------------
# 

for iyear in range(year_start,year_end+1):
  for imonth in range(month_start,month_end+1):

#
# -------------------------------------------------
# Loop on variables names
# -------------------------------------------------
# 

    for k in range(len(variables)):

#
# Variable's name, long-name and level-type
#
        vname = variables[k]
        vlong = era5[vname][0]

        print('  Processing variable: '+vname)

#
# Read input filedate.toordinal(date(Yorig,1,1))
#

        fname_in = era5_dir_raw + '/ERA5_ecmwf_' + vname.upper() + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
        nc = netcdf(fname_in,'r+',format='NETCDF4')
        time = nc.variables['time'][:]
        lat = nc.variables['latitude'][:]
        lon = nc.variables['longitude'][:]
        data = nc.variables[vname][:,:,:]
        nc.close()

#
# Flip latitudes (to have increasing latitudes...)
#

        lat = np.flip(lat, axis=0)
        data = np.flip(data, axis=1)

#
# Missing values and multiply by cff to change unit
#

        try:
            mvalue=data.fill_value
        except AttributeError:
            print ('No fill value.. use nan')
            mvalue=np.nan

        data=np.array(data)
        data=conv_cff[k]*data
        data[np.where(data==mvalue)]=9999.

#
# Convert time from hours since 1900-1-1 0:0:0 into days since Yorig-1-1 0:0:0
#

        time = time / 24.
        time = time - date.toordinal(date(Yorig,1,1)) \
	            + date.toordinal(date(1900,1,1))

#
# Changes names
#

        if vname=='u10':
            vname='u10m'

        if vname=='v10':
            vname='v10m'
#
# Create and write output netcdf file
#

        fname_out = era5_dir_processed + '/' + vname.upper() + '_Y' + str(iyear) + 'M' + str(imonth) + '.nc'

        nw = netcdf(fname_out,mode='w',format='NETCDF4')

        dimlon  = nw.createDimension('lon',  len(lon))
        dimlat  = nw.createDimension('lat',  len(lat))
        dimtime = nw.createDimension('time', None)

        varlon = nw.createVariable('lon', 'f4',('lon',))
        varlat = nw.createVariable('lat', 'f4',('lat',))
        vartime = nw.createVariable('time', 'f4',('time',))
        vardata = nw.createVariable(vname.upper(), 'f4',('time','lat','lon'))
        varlon.long_name = 'longitude of RHO-points'
        varlat.long_name = 'latitude of RHO-points'
        vartime.long_name = 'Time'
        varlon.units = 'degree_east'
        varlat.units = 'degree_north'
        vartime.units = 'days since '+str(Yorig)+'-1-1'
        vardata.missing_value = 9999.
        vardata.units = units[k]
        vardata.long_name = vlong
	
        varlon[:]=lon
        varlat[:]=lat
        vartime[:]=time
        vardata[:]=data
	
        nw.close()


# Print last message on screen
print(' ')
print(' ERA5 files conversion done ')
print(' ')




