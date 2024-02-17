# -*- coding: utf-8 -*-
"""

"""
from math import floor,ceil
import numpy as np
from datetime import datetime, timedelta, date, time
import calendar
import sys, os
from pathlib import Path
     
def download_glorys_monthly(usrname, passwd, domain, downloadDate, varList, depths, outputDir):
 
    """
    Download a subset of MERCATOR 1/12 deg reanalysis data (GLORYS) for a calendar month from
    http://marine.copernicus.eu/services-portfolio/access-to-products/
    
    """
    
    os.makedirs(outputDir,exist_ok=True)

    # loop on the variables to generate the variable string so that the number of requested variables is flexible
    # so = Salinity in psu, thetao = Temperature in degrees C, zos = SSH in m, uo = Eastward velocity in m/s, vo = Northward velocity in m/s
    var_str=''
    for var in varList:
        var_str=var_str+' -v '+var  
    
    day_end = str(calendar.monthrange(downloadDate.year, downloadDate.month)[1])

    # output filename
    fname = str(downloadDate.strftime('%Y_%m'))+'.nc'
    
    runcommand = 'copernicus-marine subset -i cmems_mod_glo_phy_my_0.083_P1D-m'+ \
            ' -x '+str(domain[0])+' -X '+str(domain[1])+ \
            ' -y '+str(domain[2])+' -Y '+str(domain[3])+ \
            ' -t \''+str(downloadDate.strftime('%Y-%m-01 00:00:00'))+'\''+ \
            ' -T \''+str(downloadDate.strftime('%Y-%m-'))+day_end+' 23:59:59'+'\''+ \
            ' -z '+str(depths[0])+' -Z '+str(depths[1])+ \
                var_str+ \
            ' -o '+outputDir+' -f '+fname+ \
            ' --force-download'+ \
            ' --username '+usrname+' --password '+passwd
    
    if os.path.exists(outputDir+fname)==False:
    	# run the runcommand, i.e. download the data specified above
        # (could do some error handling here in case the download fails)
        print(runcommand)
        os.system(runcommand)
    else:
        print(outputDir+fname+' already exists - not downloading again')

if __name__ == "__main__":
    
    # Mercator glorys product
    #
    # your CMEMS username and password are command line arguments to this script
    usrname = sys.argv[1]
    passwd = sys.argv[2]
    # spatial extent
    domain = [23, 34, -37, -31]
    outputDir='./algoa/'
    # variables to extract
    varList = ['so', 'thetao', 'zos', 'uo', 'vo']
    # min and max depths
    depths = [0.493, 5727.918]
    # download the data
    print('\n'+'downloading mercator data...')
    #start_date = '1993-01'
    #end_date = '2019-12'
    start_date = sys.argv[3]
    end_date = sys.argv[4]

    # Convert the start and end dates to datetime objects
    start_date = datetime.strptime(start_date, '%Y-%m')
    end_date = datetime.strptime(end_date, '%Y-%m')
    
    downloadDate=start_date
    while downloadDate <= end_date:
        print(downloadDate.strftime('%Y-%m'))
        download_glorys(usrname, passwd, domain, downloadDate, varList, depths, outputDir)
        downloadDate=downloadDate+timedelta(days=32) # 32 days ensures we get to the next month
        downloadDate=datetime(downloadDate.year, downloadDate.month, 1) # set the first day of the month

