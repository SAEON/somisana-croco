from datetime import datetime, timedelta
import os
from pathlib import Path
import urllib

"""
Download GFS forecast data
The GFS model is initialised every 6 hours, and provides hourly forecasts
For the historical data we download the forecast for hours 1 through 6 from each initialisation 
The forecast data gets downloaded from the latest available initialisation

The latest initialisation is used at the time this is run, so that if there is an error on GFS
side then a slightly older initialization will be found
"""

def time_param(dt):
    return dt.strftime("%Y%m%d") + "/" + dt.strftime("%H") + "/atmos"

def create_fname(dt, i):
    return dt.strftime("%Y%m%d") + dt.strftime("%H") + "_f" + str(i).zfill(3) + ".grb"

def validate_download_or_remove(fileout):
    if Path(fileout).stat().st_size < 1000:
        print(
            "WARNING:", fileout, "< 1kB (flagged as invalid)", open(fileout, "r").read()
        )
        os.remove(fileout)
        return False
    else:
        return True

def set_params(_params, dt, i):
    params = dict(_params)
    params["file"] = "gfs.t{h}{z}{f}".format(
        h=dt.strftime("%H"), z="z.pgrb2.0p25.f", f=str(i).zfill(3)
    )
    params["dir"] = "/gfs.{t}".format(t=time_param(dt))
    # return params
    return urllib.parse.urlencode(params)  # Encode the parameters

def download_file(fname, outputDir, encoded_params):
    url = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?" + encoded_params  # Construct URL
    fileout = os.path.join(outputDir, fname)
    if not os.path.isfile(fileout):
        max_retries = 3
        delay = 60
        download_success = False
        for attempt in range(1,max_retries+1):
            try:
                now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[{now}] Downloading {fileout}")
                response = urllib.request.urlopen(url)  # Fetch data
                with open(fileout, 'wb') as f:
                    f.write(response.read())
                download_success = validate_download_or_remove(fileout)

            except urllib.error.URLError as e:
                print(f"Download of {fileout} failed: {e}")

            if download_success:
                return
            else:
                print(f"Retrying {fileout} download in {delay} seconds...")
    else:
        print("File already exists", fileout)

def get_latest_available_dt(dt):
    latest_available_date = datetime(dt.year, dt.month, dt.day, 18, 0, 0)
    gfs_exists = False
    iters = 0
    
    def check_file_exists(url):
        try:
            response = urllib.request.urlopen(url)  # Attempt to open the URL
            return response.status == 200         # Check for HTTP 200 OK 
        except Exception:  # Catch potential network-related errors
            return False

    while not (gfs_exists):
        if iters > 4:
            print("GFS data is not presently available")
            exit(1)

        dataset_url = (
            "https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs"
            + latest_available_date.strftime("%Y%m%d")
            + "/gfs_0p25_"
            + latest_available_date.strftime("%H")
            + "z"
        )

        print("Testing GFS availability", dataset_url)
        response = urllib.request.urlopen(dataset_url)
        data = response.read().decode('utf-8')  # Read the response data
        if "is not an available service" in data: 
            latest_available_date = latest_available_date + timedelta(hours=-6)
            iters += 1
        else:
            # Assume valid if error message is not found
            print(
                "Latest available GFS initialisation found at",
                dataset_url,
                "\n",
                "\n",
            )
            gfs_exists = True

    return latest_available_date

def download_hindcast(start, end, outputDir, params):
    while start < end:
        for i in range(1, 7):  # hours 1 to 6
            download_file(create_fname(start, i), outputDir, set_params(params, start, i))
        start = start + timedelta(hours=6)

def download_forecast(total_forecast_hours, latest_available_date, outputDir, params):
    for i in range(1, total_forecast_hours + 1):
        download_file(create_fname(latest_available_date, i), outputDir, set_params(params, latest_available_date, i))

def download_gfs_atm(domain, run_date, hdays, fdays, outputDir):
    _now = datetime.now()
    hdays = hdays + 0.25
    fdays = fdays + 0.25
    start_date = run_date + timedelta(days=-hdays)

    latest_available_date = get_latest_available_dt(run_date)
    #latest_available_date = run_date ## FOR DEBUGGING - DELETE
    delta_days = (latest_available_date - run_date).total_seconds() / 86400

    params = {
        "lev_10_m_above_ground": "on",
        "lev_2_m_above_ground": "on",
        "lev_surface": "on",
        "lev_mean_sea_level": "on",
        "var_DLWRF": "on",
        "var_DSWRF": "on",
        "var_LAND": "on",
        "var_PRATE": "on",
        "var_RH": "on",
        "var_SPFH": "on",
        "var_TMP": "on",
        "var_UFLX": "on",
        "var_UGRD": "on",
        "var_ULWRF": "on",
        "var_USWRF": "on",
        "var_VFLX": "on",
        "var_VGRD": "on",
        "var_PRMSL": "on",
        "subregion": "",
        "leftlon": str(domain[0]),
        "rightlon": str(domain[1]),
        "toplat": str(domain[3]),
        "bottomlat": str(domain[2]),
    }

    # Download forcing files up to latest available date
    print("\nDOWNLOADING HINDCAST files")
    download_hindcast(start_date, latest_available_date, outputDir, params)

    # Download forecast forcing files
    print("\nDOWNLOADING FORECAST files")
    total_forecast_hours = int((fdays - delta_days) * 24)
    download_forecast(total_forecast_hours, latest_available_date, outputDir, params)    

    print("GFS download completed (in " + str(datetime.now() - _now) + " h:m:s)")
    
    # write a text .env file containing the delta_days variable
    # as we need it in subsequent steps for preparing model input files
    gfs_env = outputDir + '/gfs.env'
    with open(gfs_env, "w+") as env:
        env.writelines(
            [
                "RUN_DATE=" + str(datetime.strftime(run_date, "%Y-%m-%d %H")) + "\n",
                "DELTA_DAYS_GFS=" + str(delta_days) + "\n",
                "HDAYS=" + str(hdays) + "\n",
                "FDAYS=" + str(fdays) + "\n",
            ]
        )
    try:
        os.chmod(gfs_env, 0o777)
    except:
        print(
            "Unable to run os.chmod on the .env file - you may have to do this manually"
        )
    
