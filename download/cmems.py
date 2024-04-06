"""
Scripts for downloading data from the Copernic Marine Service
http://marine.copernicus.eu/services-portfolio/access-to-products/
"""
from math import floor,ceil
import numpy as np
from datetime import datetime, timedelta, date
import calendar
import sys, os
from pathlib import Path
import xarray as xr
import asyncio
import subprocess
import time
# the lib imports are from the 'lib' directory
# which is copied in from the original somisana repo Zach set up
from download.lib.open_datasets import open_datasets
from download.lib.log import log
# import copernicusmarine (couldn't use this due to compatibility issues with my environment?)

def download_mercator(usrname, passwd, domain, run_date, hdays, fdays, outputDir):
    """
    Download the operational Mercator ocean output
    Based largely on Zach's download code, but removing the asyncronous download functionality
    """
    hdays = hdays + 1
    fdays = fdays + 1
    start_date = run_date + timedelta(days=-hdays)
    end_date = run_date + timedelta(days=fdays)
    
    def is_valid_netcdf_file(file_path):
        try:
            with xr.open_dataset(file_path) as ds:
                return True
        except:
            return False
    
    # hard coding the variables to extract
    VARIABLES = [
        {
            "name": "so",
            "#": "Salinity in psu",
            "id": "cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m",
            "vars": ["so"],
        },
        {
            "name": "thetao",
            "#": "Temperature in degrees C",
            "id": "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m",
            "vars": ["thetao"],
        },
        {
            "name": "zos",
            "#": "SSH in m",
            "id": "cmems_mod_glo_phy_anfc_0.083deg_P1D-m",
            "vars": ["zos"],
        },
        {
            "name": "uo_vo",
            "#": "uo:Eastward velocity in m/s | vo:Northward velocity in m/s",
            "id": "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m",
            "vars": ["uo", "vo"],
        },
    ]
    
    # all depths
    DEPTH_RANGE = [0.493, 5727.918]

    MAX_RETRIES = 3
    RETRY_WAIT = 10

    def run_cmd(c, run_date, start_date, end_date, domain, outputDir):
        dataset = c["id"]
        name = c["name"]
        vars = c["vars"]
        variables = f"--variable {' --variable '.join(vars)} "
        fname = f"mercator_{name}_{run_date.strftime('%Y%m%d')}.nc"
        c["fname"] = fname

        runcommand = f"""
            copernicusmarine subset -i {dataset} \
                --force-download \
                --username {usrname} \
                --password {passwd} \
                -x {domain[0]} \
                -X {domain[1]} \
                -y {domain[2]} \
                -Y {domain[3]} \
                -t "{start_date.strftime("%Y-%m-%d")}" \
                -T "{end_date.strftime("%Y-%m-%d")}" \
                -z {DEPTH_RANGE[0]} \
                -Z {DEPTH_RANGE[1]} \
                {variables} \
                -o {os.path.normpath(outputDir)} \
                -f {fname}"""

        log(" ".join(runcommand.split()))
        
        f = os.path.normpath(os.path.join(outputDir, fname))
        if os.path.exists(f) == False:
            for i in range(MAX_RETRIES):
                log(
                    f"downloading latest mercator ocean forecast from CMEMS. Attempt {i + 1} of {MAX_RETRIES}"
                )
                try: 
                    os.system(runcommand)
                    if is_valid_netcdf_file(f):
                        log("Completed", name)
                    else:
                        os.unlink(f)
                        raise Exception(f"Mercator download failed (bad NetCDF output) for variable {name}")
                except Exception as e:  # Catch all potential exceptions here
                    log(f"Error: {e}, retrying in {RETRY_WAIT} seconds...")
                    time.sleep(RETRY_WAIT)
                    raise  # Re-raise the exception for potential retries
                    
            log(f"Failed to download after 3 attempts: {name}")
        else:
            log(
                os.path.normpath(os.path.join(outputDir, fname)),
                "already exists - not downloading mercator data",
            )
    
    def get_path(f, outputDir):
        return os.path.abspath(os.path.join(outputDir, f))

    # Download the separate NetCDF files for each variable
    for c in VARIABLES:
        run_cmd(c, run_date, start_date, end_date, domain, outputDir)

    # Concatenate the separate NetCDF files into the expected single file structure
    log("concatenating NetCDF files")
    output_path = os.path.abspath(
        os.path.join(outputDir, f"mercator_{run_date.strftime('%Y%m%d')}.nc")
    )
    with open_datasets(
        *[get_path(el["fname"], outputDir) for el in VARIABLES]
    ) as datasets:
        encoding = {
            "zos": {
                "dtype": "float32",
                "scale_factor": 0.000305185094475746,
                "add_offset": 0.0,
            },
            "so": {
                "dtype": "float32",
                "add_offset": -0.00152592547237873,
                "scale_factor": 0.00152592547237873,
            },
            "uo": {
                "dtype": "float32",
                "scale_factor": 0.000610370188951492,
                "add_offset": 0.0,
            },
            "vo": {
                "dtype": "float32",
                "scale_factor": 0.000610370188951492,
                "add_offset": 0.0,
            },
            "thetao": {
                "dtype": "float32",
                "scale_factor": 0.000732444226741791,
                "add_offset": 21.0,
            },
        }
        with xr.merge(datasets) as ds:
            ds.to_netcdf(
                output_path, mode="w", encoding=encoding, format="NETCDF3_CLASSIC"
            )

    subprocess.call(["chmod", "-R", "775", output_path])
    
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
        
        runcommand = 'copernicusmarine subset -i cmems_mod_glo_phy_my_0.083deg_P1D-m'+ \
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
    
    
