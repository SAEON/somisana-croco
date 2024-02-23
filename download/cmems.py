"""
Scripts for downloading data from the Copernic Marine Service
http://marine.copernicus.eu/services-portfolio/access-to-products/
"""
from math import floor,ceil
import numpy as np
from datetime import datetime, timedelta, date, time
import calendar
import sys, os
from pathlib import Path
     
def download_glorys(usrname, 
                            passwd, 
                            domain, 
                            start_date,
                            end_date,
                            varList, 
                            depths, 
                            outputDir):
 
    """
    Download month by month of daily MERCATOR 1/12 deg reanalysis data (GLORYS)
    
    """
    
    os.makedirs(outputDir,exist_ok=True)

    downloadDate=start_date
    
    while downloadDate <= end_date:
        
        print(downloadDate.strftime('%Y-%m'))
        
        # loop on the variables to generate the variable string so that the number of requested variables is flexible
        # so = Salinity in psu, thetao = Temperature in degrees C, zos = SSH in m, uo = Eastward velocity in m/s, vo = Northward velocity in m/s
        var_str=''
        for var in varList:
            var_str=var_str+' -v '+var  
        
        day_end = str(calendar.monthrange(downloadDate.year, downloadDate.month)[1])
    
        # output filename
        fname = str(downloadDate.strftime('%Y_%m'))+'.nc'
        
        runcommand = 'copernicusmarine subset -i cmems_mod_glo_phy_my_0.083_P1D-m'+ \
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
            
        downloadDate=downloadDate+timedelta(days=32) # 32 days ensures we get to the next month
        downloadDate=datetime(downloadDate.year, downloadDate.month, 1) # set the first day of the month


if __name__ == "__main__":
    
    # Mercator glorys product
    #
    # your CMEMS username and password
    usrname = 'your_usrname'
    passwd = 'your_passwd'
    # spatial extent
    domain = [23, 34, -37, -31]
    outputDir='./algoa/'
    # variables to extract
    varList = ['so', 'thetao', 'zos', 'uo', 'vo']
    # min and max depths
    depths = [0.493, 5727.918]
    # time to download
    start_date = datetime(1993,1,1)
    end_date = datetime(2019,12,1)
    download_glorys(usrname, passwd, domain, start_date, end_date, varList, depths, outputDir)
    
    
