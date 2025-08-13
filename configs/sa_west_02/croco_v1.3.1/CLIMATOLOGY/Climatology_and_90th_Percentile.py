#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 12:08:35 2025

@author: nc.memela
"""

import os
import glob
import xarray as xr
import numpy as np
from datetime import datetime

# === CONFIG ===
INPUT_DIR = "/media/nc.memela/Leviathan/MODELS/somisana-croco/configs/SA_West/croco_v1.3.1/C01_I01_GLORYS_ERA5/output"
OUT_DIR = "/home/nc.memela/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/CLIMATOLOGY"
os.makedirs(OUT_DIR, exist_ok=True)

LOGFILE = os.path.join(OUT_DIR, "climatology_generation_python.log")

def log(msg):
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    line = f"{timestamp} {msg}"
    print(line)
    with open(LOGFILE, "a") as f:
        f.write(line + "\n")

log("Climatology generation started.")

mean_outputs = []
p90_outputs = []

for month in range(1, 13):
    month_str = f"{month:02d}"
    log(f"Processing month {month_str}...")

    files = sorted(glob.glob(os.path.join(INPUT_DIR, f"croco_avg_Y*{month_str}.nc")))
    if not files:
        log(f"[WARN] No files found for month {month_str}. Skipping.")
        continue

    # Validate files
    valid_files = []
    for f in files:
        try:
            xr.open_dataset(f, engine="netcdf4").close()
            valid_files.append(f)
        except Exception as e:
            log(f"[SKIP] Failed to open {f}: {e}")

    if not valid_files:
        log(f"[ERROR] No valid files to process for month {month_str}. Skipping.")
        continue

    # Output paths
    mean_out = os.path.join(OUT_DIR, f"monthly_mean_M{month_str}.nc")
    p90_out = os.path.join(OUT_DIR, f"monthly_p90_M{month_str}.nc")

    # Skip if already done
    if os.path.exists(mean_out) and os.path.exists(p90_out):
        log(f"[SKIP] Outputs already exist for month {month_str}.")
        mean_outputs.append(mean_out)
        p90_outputs.append(p90_out)
        continue

    try:
        log(f"[INFO] Opening {len(valid_files)} files with xarray...")
        ds = xr.open_mfdataset(
            valid_files,
            combine='by_coords',
            parallel=False,
            chunks={},  # lazy loading
            engine='netcdf4'
        )

        # Compute and save mean
        if not os.path.exists(mean_out):
            log(f"[INFO] Computing monthly mean for {month_str}...")
            mean_ds = ds.mean(dim="time", keep_attrs=True)
            mean_ds.to_netcdf(mean_out)
            mean_outputs.append(mean_out)

        # Compute and save 90th percentile
        if not os.path.exists(p90_out):
            log(f"[INFO] Computing 90th percentile (30-day rolling)...")
            window = 30
            ds_rolling = ds.rolling(time=window, center=True).construct("window_dim")
            p90_ds = ds_rolling.reduce(np.nanpercentile, q=90, dim="window_dim", keep_attrs=True)
            p90_ds.to_netcdf(p90_out)
            p90_outputs.append(p90_out)

        ds.close()
    except Exception as e:
        log(f"[ERROR] Failed to process month {month_str}: {e}")
        continue

# Final merge of all means and percentiles
try:
    log("[INFO] Merging all monthly mean climatologies...")
    mean_combined = xr.concat([xr.open_dataset(f) for f in mean_outputs], dim="month")
    mean_combined.to_netcdf(os.path.join(OUT_DIR, "monthly_mean_clim.nc"))
    log("[DONE] monthly_mean_clim.nc created.")

    log("[INFO] Merging all monthly 90th percentile climatologies...")
    p90_combined = xr.concat([xr.open_dataset(f) for f in p90_outputs], dim="month")
    p90_combined.to_netcdf(os.path.join(OUT_DIR, "monthly_p90_clim.nc"))
    log("[DONE] monthly_p90_clim.nc created.")
except Exception as e:
    log(f"[ERROR] Failed during final merge: {e}")

log("✅ Climatology generation completed.")
