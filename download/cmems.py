"""
Scripts for downloading data from the Copernic Marine Service (CMEMS)
http://marine.copernicus.eu/services-portfolio/access-to-products/
"""
from math import floor,ceil
import numpy as np
from datetime import datetime, timedelta, date
import calendar
import sys, os
from pathlib import Path
import xarray as xr
import subprocess
import time
import threading

def is_valid_netcdf_file(file_path):
    try:
        with xr.open_dataset(file_path) as ds:
            return True
    except:
        return False

def download_cmems(usrname, passwd, dataset, varlist, start_date, end_date, domain, depths, outputDir, fname, ver=''):
    """
    Generic function to download a subset of a CMEMS dataset
    This is called by other functions in this file
    Input variables should be self-explanatory from the runcommand
    """
    
    # skip this file if it already exists
    f = os.path.normpath(os.path.join(outputDir, fname))
    if is_valid_netcdf_file(f):
        print("file already exists - "+ fname)
        return  
    
    variables = f"-v {' -v '.join(varlist)} "

    if ver: # only use if it's a defined input
        version = f"--dataset-version {ver}"
    else:
        version = ""

    runcommand = f"""
        copernicusmarine subset -i {dataset} \
            {version} \
            --username {usrname} \
            --password {passwd} \
            -x {domain[0]} \
            -X {domain[1]} \
            -y {domain[2]} \
            -Y {domain[3]} \
            -t "{start_date.strftime("%Y-%m-%d 00:00:00")}" \
            -T "{end_date.strftime("%Y-%m-%d 23:59:59")}" \
            -z {depths[0]} \
            -Z {depths[1]} \
            {variables} \
            -o {os.path.normpath(outputDir)} \
            -f {fname}"""

    # allow for a few retries if there was a temporary download error
    MAX_RETRIES = 3
    RETRY_WAIT = 10      
    
    i = 0
    while i < MAX_RETRIES:
        print(
            f"Attempt {i+1} of {MAX_RETRIES}"
        )
        try: 
            os.system(runcommand)
            if is_valid_netcdf_file(f):
                print("Completed "+fname)
                break
            else:
                os.unlink(f)
                raise Exception(f"Mercator download failed (bad NetCDF output): {fname}")
        except Exception as e:  # Catch all potential exceptions here
            print(f"Error: {e}, retrying in {RETRY_WAIT} seconds...")
            time.sleep(RETRY_WAIT)
            raise  # Re-raise the exception for potential retries
        i+=1

def download_mercator(usrname, passwd, domain, run_date, hdays, fdays, outputDir):
    """
    Download the operational Mercator ocean output
    """
    # extend the download range by a day either side so we are guarenteed to cover the required model run time
    hdays = hdays + 1
    fdays = fdays + 1
    start_date = run_date + timedelta(days=-hdays)
    end_date = run_date + timedelta(days=fdays)
    
    # hard coding the variable info to extract (variables are in separate files)
    VARIABLES = [
        {
            "name": "so",
            "#": "Salinity in psu",
            "id": "cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m",
            "vars": ["so"],
            "fname": f"mercator_so_{run_date.strftime('%Y%m%d_%H')}.nc"
        },
        {
            "name": "thetao",
            "#": "Temperature in degrees C",
            "id": "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m",
            "vars": ["thetao"],
            "fname": f"mercator_thetao_{run_date.strftime('%Y%m%d_%H')}.nc"
        },
        {
            "name": "zos",
            "#": "SSH in m",
            "id": "cmems_mod_glo_phy_anfc_0.083deg_P1D-m",
            "vars": ["zos"],
            "fname": f"mercator_zos_{run_date.strftime('%Y%m%d_%H')}.nc"
        },
        {
            "name": "uo_vo",
            "#": "uo:Eastward velocity in m/s | vo:Northward velocity in m/s",
            "id": "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m",
            "vars": ["uo", "vo"],
            "fname": f"mercator_uo_vo_{run_date.strftime('%Y%m%d_%H')}.nc"
        },
    ]
    
    # all depths
    depths = [0.493, 5727.918]
    
    # you can loop on the variables to download them in series
    # for var in VARIABLES:
    #     download_cmems(usrname, passwd,var["id"], var["vars"], start_date, end_date, domain, depths, outputDir, var["fname"])
    
    # but I'm rather doing them in parallel to save time (thanks Gemini)
    def download_worker(var):
        download_cmems(usrname, passwd, var["id"], var["vars"], start_date, end_date, domain, depths, outputDir, var["fname"])
    threads = []
    for var in VARIABLES:
        t = threading.Thread(target=download_worker, args=(var,))
        threads.append(t)
        t.start()
    # Wait for all threads to finish
    for t in threads:
        t.join()
    
    # Concatenate the separate NetCDF files
    print("concatenating NetCDF files")
    output_path = os.path.abspath(os.path.join(outputDir, f"MERCATOR_{run_date.strftime('%Y%m%d_%H')}.nc"))
    with xr.open_mfdataset([os.path.abspath(os.path.join(outputDir, var["fname"])) for var in VARIABLES]) as ds:
        ds.to_netcdf(output_path, mode="w")

    subprocess.call(["chmod", "-R", "775", output_path])
    
def download_glorys(usrname, 
                            passwd, 
                            domain, 
                            start_date,
                            end_date,
                            varlist, 
                            depths, 
                            outputDir):
 
    """
    Download month by month of daily MERCATOR 1/12 deg reanalysis data (GLORYS)
    
    This is actually superseded by download_monthly()
    
    """
   
    if start_date >= datetime(2021,7,1):
        dataset = 'cmems_mod_glo_phy_myint_0.083deg_P1D-m'
    else:
        dataset = 'cmems_mod_glo_phy_my_0.083deg_P1D-m'

    os.makedirs(outputDir,exist_ok=True)

    downloadDate=start_date
    
    while downloadDate <= end_date:
        
        print(downloadDate.strftime('%Y-%m'))
        
        # start and end days of this month
        start_date_download=datetime(downloadDate.year,downloadDate.month,1)
        day_end = calendar.monthrange(downloadDate.year, downloadDate.month)[1]
        end_date_download=datetime(downloadDate.year,downloadDate.month,day_end)
    
        # output filename
        fname = str(downloadDate.strftime('%Y_%m'))+'.nc'
        
        download_cmems(usrname, passwd, dataset, varlist, start_date_download, end_date_download, domain, depths, outputDir, fname)
            
        downloadDate=downloadDate+timedelta(days=32) # 32 days ensures we get to the next month
        downloadDate=datetime(downloadDate.year, downloadDate.month, 1) # set the first day of the month

def download_cmems_monthly(usrname, 
                            passwd,
                            dataset,
                            domain,
                            start_date,
                            end_date,
                            varlist, 
                            depths, 
                            outputDir):
 
    """
    Download month by month for any dataset on CMEMS
    
    """
   
    os.makedirs(outputDir,exist_ok=True)

    downloadDate=start_date
    
    while downloadDate <= end_date:
        
        print(downloadDate.strftime('%Y-%m'))
        
        # start and end days of this month
        start_date_download=datetime(downloadDate.year,downloadDate.month,1)
        day_end = calendar.monthrange(downloadDate.year, downloadDate.month)[1]
        end_date_download=datetime(downloadDate.year,downloadDate.month,day_end)
    
        # output filename
        fname = str(downloadDate.strftime('%Y_%m'))+'.nc'
        
        download_cmems(usrname, passwd, dataset, varlist, start_date_download, end_date_download, domain, depths, outputDir, fname)
            
        downloadDate=downloadDate+timedelta(days=32) # 32 days ensures we get to the next month
        downloadDate=datetime(downloadDate.year, downloadDate.month, 1) # set the first day of the month

#TODO I think we should create a download_cmems_ops function which works in a similar
# way to download_mercator (with the hdays and fdays option) but is generic for any cmems dataset
# A lot of the code in download_mercator can be moved to the cli function, which would call the 
# generic download_cmems_ops function

if __name__ == "__main__":
    
    # Mercator glorys product
    #:wq
    # your CMEMS username and password
    usrname = 'gfearon'
    passwd = 'Isabella_2015'
    # spatial extent
    domain = [47, 59, 22, 32]
    outputDir='~/CMEMS_GLO_BGC_L4/'
    # variables to extract
    varList = ['CHL', 'CHL_uncertainty', 'flags']
    # min and max depths
    depths = [0, 0]
    # time to download
    start_date = datetime(2014,1,1)
    end_date = datetime(2014,1,1)
    dataset='cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D'
    download_cmems_monthly(usrname, passwd, dataset, domain, start_date, end_date, varList, depths, outputDir)
    
    
