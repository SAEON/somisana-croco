#!/usr/bin/env python

# Script to download ECMWF ERA5 reanalysis datasets from the Climate Data
#  Store (CDS) of Copernicus https://cds.climate.copernicus.eu
#
#  This script use the CDS Phyton API[*] to connect and download specific ERA5 
#  variables, for a chosen area and monthly date interval, required by CROCO to 
#  perform simulations with atmospheric forcing. Furthermore, this script use 
#  ERA5 parameter names and not parameter IDs as these did not result in stable 
#  downloads. 
#
#  Tested using Python 3.8.6 and Python 3.9.1. This script need the following
#  python libraries pre-installed: "calendar", "datetime", "json" and "os".
#
#  [*] https://cds.climate.copernicus.eu/how-to-api
#
#  Copyright (c) DDONOSO February 2021
#  e-mail:ddonoso@dgeo.udec.cl  
#

#  You may see all available ERA5 variables at the following website
#  https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Parameterlistings

# -------------------------------------------------
# Getting libraries and utilities
# -------------------------------------------------
import cdsapi
from ERA5_utilities import *
import calendar
import datetime
import json
import os

# -------------------------------------------------
# Import my crocotools_param_python file
from era5_crocotools_param import *
print('year_start is '+str(year_start))

# -------------------------------------------------
dl=2

lonmin = str(float(lonmin)-dl)
lonmax = str(float(lonmax)+dl)
latmin = str(float(latmin)-dl)
latmax = str(float(latmax)+dl)
print ('lonmin-dl = ', lonmin)
print ('lonmax+dl =', lonmax)
print ('latmin-dl =', latmin)
print ('latmax+dl =', latmax)
# -------------------------------------------------

area = [latmax, lonmin, latmin, lonmax]

# -------------------------------------------------
# Setting raw output directory
# -------------------------------------------------
# Get the current directory
os.makedirs(era5_dir_raw,exist_ok=True)

# -------------------------------------------------
# Loading ERA5 variables's information as 
# python Dictionary from JSON file
# -------------------------------------------------
with open('ERA5_variables.json', 'r') as jf:
    era5 = json.load(jf)

# -------------------------------------------------
# Downloading ERA5 datasets
# -------------------------------------------------
# Monthly dates limits
monthly_date_start = datetime.datetime(year_start,month_start,1)
monthly_date_end = datetime.datetime(year_end,month_end,1)

# Length of monthly dates loop
len_monthly_dates = (monthly_date_end.year - monthly_date_start.year) * 12 + \
                    (monthly_date_end.month - monthly_date_start.month) + 1

# Initial monthly date
monthly_date = monthly_date_start

# Monthly dates loop
for j in range(len_monthly_dates):

    # Year and month
    year = monthly_date.year;
    month = monthly_date.month;

    # Number of days in month
    days_in_month = calendar.monthrange(year,month)[1] 
    days = [f"{day:02}" for day in range(1, days_in_month + 1)]

    # Variables/Parameters loop
    for k in range(len(variables)):

        # Variable's name, long-name and level-type
        vname = variables[k]
        vlong = era5[vname][0]

        # Request options
        options = {
             'product_type': ['reanalysis'],
             'variable': [vlong],
             'year': [str(year)],
             'month': [str(month)],
             'day': days,
             # 'time': times,
             'data_format': 'netcdf',
             'download_format': 'unarchived',
             'area': area,
             }
		  
        # # Add options to Variable without "diurnal variations"
        if vlong == 'sea_surface_temperature':
            options['time'] = ['00:00']
   
        elif vlong == 'land_sea_mask':
            options['time'] = ['00:00']

        else:
           options['time'] = times

        # Add options to Product "pressure-levels"
        if vlong == 'specific_humidity' or vlong == 'relative_humidity':
           options['pressure_level'] = ['1000']
           product = 'reanalysis-era5-pressure-levels'

        # Product "single-levels"
        else:
           product = 'reanalysis-era5-single-levels'

        # Output filename
        fname = 'ERA5_ecmwf_' + vname.upper() + '_Y' + str(year) + 'M' + str(month).zfill(2) + '.nc'
        output = era5_dir_raw + '/' + fname

        # Information strings
        info_time_clock = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        info_monthly_date = monthly_date.strftime('%Y-%b')
        info_n_overlap = ' with ' + str(n_overlap) + ' overlapping day(s) '

        # Printing message on screen
        print('                                                           ')
        print('-----------------------------------------------------------')
        print('',info_time_clock,'                                        ')
        print(' Performing ERA5 data request, please wait...              ')
        print(' Date [yyyy-mmm] =',info_monthly_date + info_n_overlap,'   ')
        print(' Variable =',vlong,'                                       ')
        print('-----------------------------------------------------------')
        print('Request options: ')
        print(options)
      
        # Server ECMWF-API
        c = cdsapi.Client()

        # Do the request
        c.retrieve(product, options).download(output)
    
    # ---------------------------------------------------------------------
    # Next iteration to monthly date: add one month to current monthly date
    # ---------------------------------------------------------------------
    monthly_date = addmonths4date(monthly_date,1)


# Print output message on screen
print('                                               ')
print(' ERA5 data request has been done successfully! ')
print('                                               ')




