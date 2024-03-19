import asyncio
from datetime import datetime, timedelta, time
import os
from pathlib import Path
from datetime import datetime, timedelta
import requests as r
import aiofiles
import aiohttp

"""
Download GFS forecast data
The GFS model is initialised every 6 hours, and provides hourly forecasts
For the historical data we download the forecast for hours 1 through 6 from each initialisation 
The forecast data gets downloaded from the latest available initialisation

The latest initialisation is used at the time this is run, so that if there is an error on GFS
side then a slightly older initialization will be found
"""

def yyyymmdd(dt):
    return dt.strftime("%Y") + dt.strftime("%m") + dt.strftime("%d")


def time_param(dt):
    return yyyymmdd(dt) + "/" + dt.strftime("%H") + "/atmos"


def create_fname(dt, i):
    return yyyymmdd(dt) + dt.strftime("%H") + "_f" + str(i).zfill(3) + ".grb"


def validate_download_or_remove(fileout):
    if Path(fileout).stat().st_size < 1000:
        print(
            "WARNING:", fileout, "< 1kB (flagged as invalid)", open(fileout, "r").read()
        )
        os.remove(fileout)


def set_params(_params, dt, i):
    params = dict(_params)
    params["file"] = "gfs.t{h}{z}{f}".format(
        h=dt.strftime("%H"), z="z.pgrb2.0p25.f", f=str(i).zfill(3)
    )
    params["dir"] = "/gfs.{t}".format(t=time_param(dt))
    return params


async def download_file(semaphore, fname, outputDir, params):
    url = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl"
    fileout = os.path.join(outputDir, fname)
    if not os.path.isfile(fileout):
        async with semaphore:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params) as response:
                    if response.status == 200:
                        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        print(f"[{now}] Downloading {fileout}")
                        async with aiofiles.open(fileout, mode="wb") as f:
                            async for chunk in response.content.iter_chunked(1024):
                                if chunk:
                                    await f.write(chunk)
                        validate_download_or_remove(fileout)
                    else:
                        print(f"Request failed with status code {response.status}")
    else:
        print("File already exists", fileout)


def get_latest_available_dt(dt):
    latest_available_date = datetime(dt.year, dt.month, dt.day, 18, 0, 0)
    gfs_exists = False
    iters = 0

    while not (gfs_exists):
        if iters > 4:
            print("GFS data is not presently available")
            exit(1)

        dataset_url = (
            "https://nomads.ncep.noaa.gov/dods/gfs_0p25_1hr/gfs"
            + yyyymmdd(latest_available_date)
            + "/gfs_0p25_1hr_"
            + latest_available_date.strftime("%H")
            + "z"
        )

        print("Testing GFS availability", dataset_url)
        result = r.head(dataset_url)
        xdap = result.headers.get("XDAP")

        if xdap:
            print(
                "Latest available GFS initialisation found at",
                dataset_url,
                "\n",
                "X-DAP HTTP Header",
                xdap,
                "\n",
                "\n",
            )
            gfs_exists = True
        else:
            latest_available_date = latest_available_date + timedelta(hours=-6)
            iters += 1

    return latest_available_date


# nomads.ncep.noaa.gov has strict rate limits,
# so keep this number low
MAX_CONCURRENT_NET_IO = 2


async def download_hindcast(semaphore, start, end, outputDir, params):
    tasks = []
    while start < end:
        for i in range(1, 7):  # hours 1 to 6
            tasks.append(
                asyncio.create_task(
                    download_file(
                        semaphore,
                        create_fname(start, i),
                        outputDir,
                        set_params(params, start, i),
                    )
                )
            )
        start = start + timedelta(hours=6)
    await asyncio.gather(*tasks)


async def download_forecast(
    semaphore, total_forecast_hours, latest_available_date, outputDir, params
):
    tasks = []
    for i in range(1, total_forecast_hours + 1):
        tasks.append(
            asyncio.create_task(
                download_file(
                    semaphore,
                    create_fname(latest_available_date, i),
                    outputDir,
                    set_params(params, latest_available_date, i),
                )
            )
        )
    await asyncio.gather(*tasks)

def download_gfs_atm(domain, run_date, hdays, fdays, outputDir):
    _now = datetime.now()
    hdays = hdays + 0.25
    fdays = fdays + 0.25
    start_date = run_date + timedelta(days=-hdays)

    latest_available_date = get_latest_available_dt(run_date)
    delta_days = (latest_available_date - run_date).total_seconds() / 86400

    params = {
        "lev_10_m_above_ground": "on",
        "lev_2_m_above_ground": "on",
        "lev_surface": "on",
        "var_DLWRF": "on",
        "var_DSWRF": "on",
        "var_LAND": "on",
        "var_PRATE": "on",
        "var_RH": "on",
        "var_TMP": "on",
        "var_UFLX": "on",
        "var_UGRD": "on",
        "var_ULWRF": "on",
        "var_USWRF": "on",
        "var_VFLX": "on",
        "var_VGRD": "on",
        "subregion": "",
        "leftlon": str(domain[0]),
        "rightlon": str(domain[1]),
        "toplat": str(domain[3]),
        "bottomlat": str(domain[2]),
    }

    # Download forcing files up to latest available date
    print("\nDOWNLOADING HINDCAST files")
    asyncio.run(
        download_hindcast(
            asyncio.BoundedSemaphore(MAX_CONCURRENT_NET_IO),
            start_date,
            latest_available_date,
            outputDir,
            params,
        )
    )

    # Download forecast forcing files
    print("\nDOWNLOADING FORECAST files")
    total_forecast_hours = int((fdays - delta_days) * 24)
    asyncio.run(
        download_forecast(
            asyncio.BoundedSemaphore(MAX_CONCURRENT_NET_IO),
            total_forecast_hours,
            latest_available_date,
            outputDir,
            params,
        )
    )

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
    
