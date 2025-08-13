#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 01:17:38 2025

@author: nc.memela
"""


import os
from datetime import datetime, timedelta
import numpy as np
import xarray as xr
import crocotools_py.postprocess as post
from pathlib import Path


def compute_climatology_and_p90(
    years,
    data_dir,
    varname,
    ref_date,
    halfwin=15,
    savefile=None,
    verbose=True,
):
    """
    Compute daily climatology and 90th percentile for a CROCO field using a ±halfwin-day window.

    Baseline preserved:
      - Outer loop = DOY → years
      - post.get_var(ds_m, time=dt) on a single month Dataset
    Efficiency:
      - Tiny cache of monthly files across DOYs (no repeated open/close)
      - Correct Dec/Jan year wrap for file selection

    Metadata:
      - Dynamically capture global attrs + coords + coord-like vars from first successful sample
      - Carry them into the output (excluding 'time' and the data var itself)
      - Copy source var attrs (e.g., units), add clear long_name
    """

    # Build path for a given year, month
    def _make_path(y, m):
        return os.path.join(data_dir, f"croco_avg_Y{y}M{m:02d}.nc")

    # Iterate all days between two dates (inclusive)
    def _iter_days(start_dt, end_dt):
        for i in range((end_dt - start_dt).days + 1):
            yield start_dt + timedelta(days=i)

    # Given a center month in the actual year, return prev/current/next with correct YEARS
    def _months_for_year(center_month, year):
        if center_month == 1:
            return [(year - 1, 12), (year, 1), (year, 2)]
        elif center_month == 12:
            return [(year, 11), (year, 12), (year + 1, 1)]
        else:
            return [(year, center_month - 1), (year, center_month), (year, center_month + 1)]

    # Resolve which YEAR a (month m) belongs to, relative to center_month in year yy
    def _resolve_year_for_md(yy, center_month, m):
        if center_month == 1 and m == 12:
            return yy - 1
        if center_month == 12 and m == 1:
            return yy + 1
        return yy

    # Leap-year template for day-of-year
    cal0 = datetime(2000, 1, 1)

    grid_shape = None
    means_per_day = []
    p90s_per_day = []

    # ---- persistent cache across DOYs ----
    _open_cache = {}  # (year, month) -> xr.Dataset

    def _ensure_open(months_needed):
        """Open any (year, month) in months_needed and keep cached (no closing here)."""
        for (yy, mm) in months_needed:
            key = (yy, mm)
            if key not in _open_cache:
                p = _make_path(yy, mm)
                if os.path.exists(p):
                    _open_cache[key] = xr.open_dataset(p, decode_times=False)

    # ---- dynamic template capture (global attrs, coords, coord-like vars) ----
    template_global_attrs = None
    template_var_attrs = None
    template_coords = None  # dict name->DataArray to attach as coordinates

    def _capture_template(ds_sel, ds_src, data_var):
        """Capture global attrs + coords/coord-like vars once, dynamically."""
        nonlocal template_global_attrs, template_var_attrs, template_coords
        if template_global_attrs is not None:
            return

        # start from source monthly file attrs; overlay any selection attrs
        template_global_attrs = {**getattr(ds_src, "attrs", {}), **getattr(ds_sel, "attrs", {})}
        # carry base attrs (units etc.) from data var
        template_var_attrs = dict(getattr(ds_sel[data_var], "attrs", {}))

        coords = {}

        def collect_from(ds):
            # true coordinates
            for name, da in ds.coords.items():
                if name == "time":
                    continue
                coords[name] = da
            # coord-like data_vars: dims subset of our standard grid dims, 1D or 2D
            for name, da in ds.data_vars.items():
                if name == data_var:
                    continue
                if name in coords:
                    continue
                if set(da.dims).issubset({"s_rho", "eta_rho", "xi_rho"}) and da.ndim in (1, 2):
                    coords[name] = da

        collect_from(ds_src)
        collect_from(ds_sel)
        template_coords = coords

    # ---- DOY loop (set to full year when ready) ----
    for doy in range(1, 367):  # change to range(1, 367) for full 366 days
        center = cal0 + timedelta(days=doy - 1)
        month_c, day_c = center.month, center.day

        start_win = center - timedelta(days=halfwin)
        end_win   = center + timedelta(days=halfwin)
        win_md = [(d.month, d.day) for d in _iter_days(start_win, end_win)]

        samples = []

        if verbose:
            print(f"[{doy:03d}/366] calendar day = {month_c:02d}-{day_c:02d}")

        # loop over years
        for yy in years:
            months_needed = _months_for_year(month_c, yy)
            _ensure_open(months_needed)

            for (m, d) in win_md:
                y_for_md = _resolve_year_for_md(yy, month_c, m)
                try:
                    dt = datetime(y_for_md, m, d)
                except ValueError:
                    continue  # e.g., Feb 29 on non-leap years

                ds_m = _open_cache.get((dt.year, dt.month))
                if ds_m is None:
                    continue

                try:
                    ds_sel = post.get_var(ds_m, varname, time=dt, ref_date=ref_date)
                    if varname not in ds_sel:
                        continue

                    # first successful sample → set shape & capture template
                    if grid_shape is None:
                        arr0 = ds_sel[varname].values
                        grid_shape = arr0.shape  # e.g., (s_rho, eta_rho, xi_rho)
                        _capture_template(ds_sel, ds_m, varname)

                    arr = ds_sel[varname].values
                    samples.append(arr.astype(np.float32, copy=False))

                except Exception as e:
                    if verbose:
                        print(f"  -> skip {dt:%Y-%m-%d}: {e}")

        # aggregate for this DOY
        if len(samples) == 0:
            if grid_shape is None:
                means_per_day.append(None)
                p90s_per_day.append(None)
            else:
                means_per_day.append(np.full(grid_shape, np.float32(np.nan)))
                p90s_per_day.append(np.full(grid_shape, np.float32(np.nan)))
        else:
            stack = np.asarray(samples, dtype=np.float32)  # (nsamples, *grid)
            mean_d = np.nanmean(stack, axis=0, dtype=np.float32)
            p90_d  = np.nanpercentile(stack, 90, axis=0).astype(np.float32)
            means_per_day.append(mean_d)
            p90s_per_day.append(p90_d)

    # ---- finalize & save ----
    # close cache before raising
    if grid_shape is None:
        for ds in _open_cache.values():
            try:
                ds.close()
            except Exception:
                pass
        raise RuntimeError("No samples were collected.")

    # ensure every day has float32 array
    means_per_day = [
        (m if m is not None else np.full(grid_shape, np.float32(np.nan))).astype(np.float32)
        for m in means_per_day
    ]
    p90s_per_day = [
        (p if p is not None else np.full(grid_shape, np.float32(np.nan))).astype(np.float32)
        for p in p90s_per_day
    ]

    # dim names
    if len(grid_shape) == 3:
        dim_names = ("s_rho", "eta_rho", "xi_rho")
    else:
        dim_names = tuple(f"dim{i}" for i in range(len(grid_shape)))

    ndoy = len(means_per_day)
    doy_coord = np.arange(1, ndoy + 1, dtype=np.int16)

    # coords: start with dayofyear, then attach captured coords/coord-like vars
    coords = {"dayofyear": ("dayofyear", doy_coord)}
    if template_coords:
        for name, da in template_coords.items():
            # include only coords whose dims are compatible with our output grid
            if all(dim in ("s_rho", "eta_rho", "xi_rho") for dim in da.dims):
                coords[name] = (da.dims, da.values, dict(getattr(da, "attrs", {})))

    # dataarrays
    clim_mean = xr.DataArray(
        np.stack(means_per_day, axis=0),
        dims=("dayofyear",) + dim_names,
        coords=coords,
        name="climatology",
    ).astype(np.float32)

    clim_p90 = xr.DataArray(
        np.stack(p90s_per_day, axis=0),
        dims=("dayofyear",) + dim_names,
        coords=coords,
        name="threshold_90",
    ).astype(np.float32)

    # variable attrs from source var (+ long_name)
    base_attrs = template_var_attrs or {}
    clim_attrs = dict(base_attrs)
    clim_attrs["long_name"] = f"{varname} daily climatology (±{halfwin}d window)"
    clim_mean.attrs = clim_attrs

    p90_attrs = dict(base_attrs)
    p90_attrs["long_name"] = f"{varname} daily 90th percentile (±{halfwin}d window)"
    clim_p90.attrs = p90_attrs

    # assemble dataset
    out = xr.Dataset({"climatology": clim_mean, "threshold_90": clim_p90}, coords=coords)

    # global attrs: template + light provenance (no hard-coding)
    gattrs = dict(template_global_attrs or {})
    hist = gattrs.get("history", "")
    extra = f"daily climatology & p90 (±{halfwin}d) from years={','.join(map(str, years))}"
    gattrs["history"] = f"{hist} | {extra}".strip(" |")
    gattrs["ref_date"] = str(ref_date)
    gattrs["halfwin_days"] = np.int64(halfwin)
    gattrs["varname"] = varname
    out.attrs = gattrs

    # save
    if savefile:
        savefile = Path(savefile).expanduser().resolve()
        savefile.parent.mkdir(parents=True, exist_ok=True)

        chunksizes = (min(64, ndoy),) + tuple(min(128, s) for s in grid_shape)
        encoding = {
            "climatology":  {"zlib": True, "complevel": 4, "dtype": "float32",
                             "chunksizes": chunksizes, "_FillValue": np.float32(np.nan)},
            "threshold_90": {"zlib": True, "complevel": 4, "dtype": "float32",
                             "chunksizes": chunksizes, "_FillValue": np.float32(np.nan)},
            "dayofyear":    {"dtype": "int16"},
        }

        out.to_netcdf(savefile, engine="h5netcdf", encoding=encoding)
        if verbose:
            print(f"Saved daily climatology and 90th percentile to {savefile}")

    # close cached datasets once
    for ds in _open_cache.values():
        try:
            ds.close()
        except Exception:
            pass

    return out


# ----------------- CALL -----------------
if __name__ == "__main__":
    years = list(range(1993, 1995))  # expand as needed
    data_dir = "/media/disk01/ERA5_run_SA_West_Hindcast"
    varname = "temp"
    ref_date = datetime(2000, 1, 1)

    OUT_DIR = Path("/home/nc.memela/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/CLIMATOLOGY/Climatology_with_90th/monthly_climatology/p90")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SAVEFILE = OUT_DIR / "climatology_with_p90_93-94_late4.nc"

    result = compute_climatology_and_p90(
        years=years,
        data_dir=data_dir,
        varname=varname,
        ref_date=ref_date,
        halfwin=15,
        savefile=str(SAVEFILE),
        verbose=True,
    )
    
#%%

import os
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr
import crocotools_py.postprocess as post


def compute_climatology_and_p90(
    years,
    data_dir,
    varname,
    ref_date,
    halfwin=15,
    savefile=None,
    verbose=True,
):
    """
    Daily climatology and 90th percentile using a ±halfwin-day window (Hobday-style)
    from CROCO monthly files that contain *daily* data.

    - Loop order preserved: day-of-year (DOY) -> years
    - For each month/year needed, load full month once via post.get_var(..., time=None)
      and slice with xarray: .sel(time=slice(start, end))
    - Mean and p90 computed over concatenated daily samples across all years in window
    """

    # ---------- helpers ----------
    def _make_path(y, m):
        return os.path.join(data_dir, f"croco_avg_Y{y}M{m:02d}.nc")

    def _months_for_year(center_month, year):
        if center_month == 1:
            return [(year - 1, 12), (year, 1), (year, 2)]
        elif center_month == 12:
            return [(year, 11), (year, 12), (year + 1, 1)]
        else:
            return [(year, center_month - 1), (year, center_month), (year, center_month + 1)]

    def _map_year(yy, center_month, other_month):
        if center_month == 1 and other_month == 12:
            return yy - 1
        if center_month == 12 and other_month == 1:
            return yy + 1
        return yy

    # capture coords and attrs once
    template_global_attrs = None
    template_var_attrs = None
    template_coords = None

    def _capture_template(ds_sel, ds_src, data_var):
        nonlocal template_global_attrs, template_var_attrs, template_coords
        if template_global_attrs is not None:
            return
        template_global_attrs = {**getattr(ds_src, "attrs", {}), **getattr(ds_sel, "attrs", {})}
        template_var_attrs = dict(getattr(ds_sel[data_var], "attrs", {}))
        coords = {}
        # true coords (except time)
        for name, da in ds_src.coords.items():
            if name != "time":
                coords[name] = da
        for name, da in ds_sel.coords.items():
            if name != "time":
                coords[name] = da
        # coord-like data_vars (grid stuff)
        for src in (ds_src, ds_sel):
            for name, da in src.data_vars.items():
                if name == data_var:
                    continue
                if name in coords:
                    continue
                if set(da.dims).issubset({"s_rho", "eta_rho", "xi_rho"}) and da.ndim in (1, 2):
                    coords[name] = da
        template_coords = coords

    # ---------- main ----------
    cal0 = datetime(2000, 1, 1)  # leap template
    grid_shape = None
    means_per_day = []
    p90s_per_day = []

    # cache of opened monthly datasets
    _open_cache = {}  # (year, month) -> xr.Dataset

    def _ensure_open(months_needed):
        for (yy, mm) in months_needed:
            key = (yy, mm)
            if key not in _open_cache:
                p = _make_path(yy, mm)
                if os.path.exists(p):
                    _open_cache[key] = xr.open_dataset(p, decode_times=False)

    for doy in range(1, 7):
        center = cal0 + timedelta(days=doy - 1)
        month_c, day_c = center.month, center.day
        start_win = center - timedelta(days=halfwin)
        end_win = center + timedelta(days=halfwin)

        if verbose:
            print(f"[{doy:03d}/366] calendar day = {month_c:02d}-{day_c:02d}")

        samples = []  # list of DataArrays with dim 'time'

        for yy in years:
            # map window bounds to the correct year (handle Dec/Jan wrap)
            yL = _map_year(yy, month_c, start_win.month)
            yR = _map_year(yy, month_c, end_win.month)
            win_start = datetime(yL, start_win.month, start_win.day)
            win_end = datetime(yR, end_win.month, end_win.day)

            months_needed = _months_for_year(month_c, yy)
            _ensure_open(months_needed)

            for (y_m, m_m) in months_needed:
                ds_m = _open_cache.get((y_m, m_m))
                if ds_m is None:
                    continue

                # month bounds
                m_first = datetime(y_m, m_m, 1)
                m_last = datetime(y_m + (m_m // 12), (m_m % 12) + 1, 1) - timedelta(days=1)

                # intersect with window
                sub_start = max(win_start, m_first)
                sub_end = min(win_end, m_last)
                if sub_start > sub_end:
                    continue

                try:
                    # IMPORTANT: get full month via post.get_var(..., time=None), then slice with xarray
                    ds_month = post.get_var(ds_m, varname, ref_date=ref_date)  # full month selection
                    da = ds_month[varname]
                    if "time" in da.dims:
                        da = da.sel(time=slice(np.datetime64(sub_start), np.datetime64(sub_end)))
                        # if empty (no days fell in slice), skip
                        if da.sizes.get("time", 0) == 0:
                            continue
                    else:
                        # rare single-record case: treat as one step at sub_start
                        da = da.expand_dims(time=[np.datetime64(sub_start)])

                    if grid_shape is None:
                        grid_shape = da.isel(time=0).shape
                        _capture_template(ds_month, ds_m, varname)

                    samples.append(da.astype(np.float32))

                except Exception as e:
                    if verbose:
                        print(f"    skip {y_m}-{m_m:02d} slice: {e}")

        # aggregate this DOY
        if not samples:
            if grid_shape is None:
                means_per_day.append(None)
                p90s_per_day.append(None)
            else:
                means_per_day.append(np.full(grid_shape, np.float32(np.nan)))
                p90s_per_day.append(np.full(grid_shape, np.float32(np.nan)))
        else:
            window_da = xr.concat(samples, dim="time").sortby("time")
            if verbose:
                print(f"  samples in window (time steps): {int(window_da.sizes['time'])}")
            mean_d = window_da.mean(dim="time", skipna=True).values.astype(np.float32)
            p90_d = window_da.quantile(0.9, dim="time", skipna=True).values.astype(np.float32)
            means_per_day.append(mean_d)
            p90s_per_day.append(p90_d)

    # finalize
    if grid_shape is None:
        for ds in _open_cache.values():
            try:
                ds.close()
            except Exception:
                pass
        raise RuntimeError("No samples were collected.")

    means_per_day = [
        (m if m is not None else np.full(grid_shape, np.float32(np.nan))).astype(np.float32)
        for m in means_per_day
    ]
    p90s_per_day = [
        (p if p is not None else np.full(grid_shape, np.float32(np.nan))).astype(np.float32)
        for p in p90s_per_day
    ]

    if len(grid_shape) == 3:
        dim_names = ("s_rho", "eta_rho", "xi_rho")
    else:
        dim_names = tuple(f"dim{i}" for i in range(len(grid_shape)))

    ndoy = len(means_per_day)
    doy_coord = np.arange(1, ndoy + 1, dtype=np.int16)

    # coords to attach
    coords = {"dayofyear": ("dayofyear", doy_coord)}
    if template_coords:
        for name, da in template_coords.items():
            if all(dim in ("s_rho", "eta_rho", "xi_rho") for dim in da.dims):
                coords[name] = (da.dims, da.values, dict(getattr(da, "attrs", {})))

    clim_mean = xr.DataArray(
        np.stack(means_per_day, axis=0),
        dims=("dayofyear",) + dim_names,
        coords=coords,
        name="climatology",
    ).astype(np.float32)

    clim_p90 = xr.DataArray(
        np.stack(p90s_per_day, axis=0),
        dims=("dayofyear",) + dim_names,
        coords=coords,
        name="threshold_90",
    ).astype(np.float32)

    base_attrs = template_var_attrs or {}
    clim_mean.attrs = {**base_attrs, "long_name": f"{varname} daily climatology (±{halfwin}d window)"}
    clim_p90.attrs = {**base_attrs, "long_name": f"{varname} daily 90th percentile (±{halfwin}d window)"}

    out = xr.Dataset({"climatology": clim_mean, "threshold_90": clim_p90}, coords=coords)

    gattrs = dict(template_global_attrs or {})
    hist = gattrs.get("history", "")
    extra = f"daily climatology & p90 (±{halfwin}d) from years={','.join(map(str, years))}"
    gattrs["history"] = f"{hist} | {extra}".strip(" |")
    gattrs["ref_date"] = str(ref_date)
    gattrs["halfwin_days"] = np.int64(halfwin)
    gattrs["varname"] = varname
    out.attrs = gattrs

    if savefile:
        savefile = Path(savefile).expanduser().resolve()
        savefile.parent.mkdir(parents=True, exist_ok=True)
        chunksizes = (min(64, ndoy),) + tuple(min(128, s) for s in grid_shape)
        encoding = {
            "climatology":  {"zlib": True, "complevel": 4, "dtype": "float32",
                             "chunksizes": chunksizes, "_FillValue": np.float32(np.nan)},
            "threshold_90": {"zlib": True, "complevel": 4, "dtype": "float32",
                             "chunksizes": chunksizes, "_FillValue": np.float32(np.nan)},
            "dayofyear":    {"dtype": "int16"},
        }
        out.to_netcdf(savefile, engine="h5netcdf", encoding=encoding)
        if verbose:
            print(f"Saved daily climatology and 90th percentile to {savefile}")

    for ds in _open_cache.values():
        try:
            ds.close()
        except Exception:
            pass

    return out


# ----------------- CALL -----------------
if __name__ == "__main__":
    years = list(range(1993, 1995))  # e.g., 1993 & 1994
    data_dir = "/media/disk01/ERA5_run_SA_West_Hindcast"
    varname = "temp"
    ref_date = datetime(2000, 1, 1)

    OUT_DIR = Path("/home/nc.memela/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/CLIMATOLOGY/Climatology_with_90th/monthly_climatology/p90")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SAVEFILE = OUT_DIR / "climatology_with_p90_window-slice.nc"

    result = compute_climatology_and_p90(
        years=years,
        data_dir=data_dir,
        varname=varname,
        ref_date=ref_date,
        halfwin=15,
        savefile=str(SAVEFILE),
        verbose=True,
    )


#%%

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---- settings ----
NC = Path("/home/nc.memela/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/CLIMATOLOGY/Climatology_with_90th/monthly_climatology/p90/climatology_with_p90_93-94_late3.nc")

# choose ONE of the two:
DAY_VALUE  = 6      # actual dayofyear value in the file (e.g., 1..366); set to None to use index
DAY_INDEX  = None   # 0-based index into dayofyear (e.g., 0 for first entry); set to None to use value

S_RHO_INDEX = -5     # 0-based depth index (0 is top if positive-up)

SAVE_PNG = None     # e.g., Path("clim_vs_p90_day001_sigma0.png") or None to just show

# ---- load ----
ds = xr.open_dataset(NC)

# select day
if DAY_VALUE is not None:
    # pick the nearest matching dayofyear value
    day = int(DAY_VALUE)
    day_idx = int((np.abs(ds["dayofyear"].astype(int) - day)).argmin())
else:
    day_idx = int(DAY_INDEX if DAY_INDEX is not None else 0)

# extract fields
clim = ds["climatology"].isel(dayofyear=day_idx, s_rho=S_RHO_INDEX)
thr  = ds["threshold_90"].isel(dayofyear=day_idx, s_rho=S_RHO_INDEX)

# ensure last two dims are (eta_rho, xi_rho)
def as_eta_xi(da):
    if list(da.dims)[-2:] != ["eta_rho", "xi_rho"]:
        return da.transpose(..., "eta_rho", "xi_rho")
    return da

clim = as_eta_xi(clim)
thr  = as_eta_xi(thr)

# lon/lat (present as variables in your file)
lon = as_eta_xi(ds["lon_rho"])
lat = as_eta_xi(ds["lat_rho"])

# apply mask if available (prefer mask_rho; fall back to mask). Treat >0 as water.
if "mask_rho" in ds:
    mask2d = as_eta_xi(ds["mask_rho"]).astype(bool)
elif "mask" in ds:
    mask2d = as_eta_xi(ds["mask"]).astype(bool)
else:
    mask2d = None

if mask2d is not None:
    clim = xr.where(mask2d, clim, np.nan)
    thr  = xr.where(mask2d,  thr, np.nan)

# shared color scale
vmin = np.nanmin([clim.min().item(), thr.min().item()])
vmax = np.nanmax([clim.max().item(), thr.max().item()])

# titles
day_val = int(ds["dayofyear"].isel(dayofyear=day_idx).item())
t_left  = f"Climatology (day {day_val}, s_rho idx {S_RHO_INDEX})"
t_right = f"90th percentile (day {day_val}, s_rho idx {S_RHO_INDEX})"

# ---- plot ----
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

m0 = ax0.pcolormesh(lon, lat, clim, shading="auto", vmin=vmin, vmax=vmax)
ax0.set_title(t_left)
ax0.set_xlabel("Longitude"); ax0.set_ylabel("Latitude")

m1 = ax1.pcolormesh(lon, lat, thr,  shading="auto", vmin=vmin, vmax=vmax)
ax1.set_title(t_right)
ax1.set_xlabel("Longitude"); ax1.set_ylabel("Latitude")

cbar = fig.colorbar(m0, ax=[ax0, ax1], orientation="vertical", fraction=0.046, pad=0.04)
units = clim.attrs.get("units") or thr.attrs.get("units") or ""
cbar.set_label(f"{ds.attrs.get('varname','temp')}" + (f" ({units})" if units else ""))

if SAVE_PNG:
    plt.savefig(SAVE_PNG, dpi=180)
else:
    plt.show()


#%%

import xarray as xr
import numpy as np

fpath = "/home/nc.memela/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/CLIMATOLOGY/Climatology_with_90th/monthly_climatology/p90/climatology_with_p90_93-94_late2.nc"
out = xr.open_dataset(fpath)

cl = out["climatology"].isel(s_rho=0)     # surface for a quick check
p9 = out["threshold_90"].isel(s_rho=0)

absdiff = np.abs(cl - p9)
print("max |clim - p90|:", float(absdiff.max()))
print("mean |clim - p90|:", float(absdiff.mean()))



# #!/usr/bin/env python3
# """
# Compute Hobday et al. 90th percentile threshold for a CROCO variable.

# Key features:
# - Loops over 366 calendar days (uses a leap base year so Feb 29 exists).
# - For each calendar day, and for each year:
#   - Open exactly THREE monthly files: previous, current, next month.
#     * January => (Dec Y-1, Jan Y, Feb Y)
#     * December => (Nov Y, Dec Y, Jan Y+1)
#   - Sample a ±window day range (default: ±15 days) around the calendar day.
#   - Collect samples across ALL years for that day.
# - Compute climatological mean and 90th percentile for each calendar day.
# - Save results to a NetCDF file (no changes to your original CROCO files).

# No grid file required. Uses crocotools_py.postprocess.get_var on already-opened datasets.
# """

# import os
# from pathlib import Path
# from datetime import datetime, timedelta
# import numpy as np
# import xarray as xr
# import crocotools_py.postprocess as post


# # ----------------- small utility helpers -----------------

# def make_file_path(year: int, month: int, data_dir: str) -> str:
#     """
#     Build the exact CROCO monthly file path matching your naming convention.

#     Example: Y=1993, M=02 -> '.../croco_avg_Y1993M02.nc'
#     """
#     return str(Path(data_dir) / f"croco_avg_Y{year}M{month:02d}.nc")


# def three_month_keys(base_month: int) -> tuple[int, int, int]:
#     """
#     Given a calendar month (1..12), return the (prev, this, next) month numbers (1..12).
#     This is purely month arithmetic; year boundaries are handled elsewhere.
#     """
#     m = base_month
#     prev_m = 12 if m == 1 else m - 1
#     next_m = 1  if m == 12 else m + 1
#     return prev_m, m, next_m


# def open_three_months_for_year(year: int, base_month: int, data_dir: str, engine: str | None = "netcdf4"):
#     """
#     Open the three monthly files needed for a given (year, base_month):
#       - prev month file (might live in year-1 when base_month is January)
#       - this month file (always in 'year')
#       - next month file (might live in year+1 when base_month is December)

#     Returns:
#       dict keyed by (yy, mm) -> xr.Dataset for files that exist.
#       (Missing files are simply skipped.)

#     Notes:
#     - We open with decode_times=False; crocotools handles time via 'ref_date'.
#     - We keep datasets open just long enough to extract the day-window,
#       then we close them explicitly.
#     """
#     prev_m, this_m, next_m = three_month_keys(base_month)

#     # Build the three (year, month) keys, adjusting year at boundaries
#     keys = []
#     keys.append((year - 1, 12) if base_month == 1  else (year, prev_m))  # prev
#     keys.append((year, this_m))                                          # this
#     keys.append((year + 1, 1)  if base_month == 12 else (year, next_m))  # next

#     opened = {}
#     for (yy, mm) in keys:
#         path = make_file_path(yy, mm, data_dir)
#         if os.path.exists(path):
#             # Open plainly; no chunking or Dask. Engine explicit for clarity.
#             ds = xr.open_dataset(path, decode_times=False, engine=engine)
#             opened[(yy, mm)] = ds
#     return opened


# def iter_days(start_d: datetime, end_d: datetime):
#     """Yield each date between start_d and end_d inclusive (day by day)."""
#     n = (end_d - start_d).days + 1
#     for i in range(n):
#         yield start_d + timedelta(days=i)


# # ----------------- main computation -----------------

# def hobday_90th_percentile(
#     years: list[int],
#     data_dir: str,
#     varname: str,
#     ref_date: datetime,
#     window: int = 15,
#     savefile: str = "hobday_p90.nc",
#     verbose: bool = False,
#     dtype=np.float32,
# ):
#     """
#     Compute the climatological mean and 90th percentile for 'varname' using a
#     ±window day sampling around each calendar day, across all specified years.

#     Parameters
#     ----------
#     years : list[int]
#         Years to include in the climatology (e.g., list(range(1993, 2020))).
#     data_dir : str
#         Directory containing CROCO monthly files named like 'croco_avg_YYYYYMmm.nc'.
#     varname : str
#         CROCO variable to extract (e.g., 'temp').
#     ref_date : datetime
#         CROCO reference datetime used by post.get_var to decode time indices.
#     window : int, default 15
#         Half-window in days; ±window yields a 2*window+1 day sampling window.
#     savefile : str, default "hobday_p90.nc"
#         Output NetCDF file path.
#     verbose : bool, default False
#         If True, prints progress and simple counters.
#     dtype : numpy dtype, default np.float32
#         Output data type for saved arrays (float32 keeps file sizes manageable).

#     Returns
#     -------
#     xarray.Dataset
#         Dataset with variables:
#           - 'climatology'  : (dayofyear, [s_rho,] eta_rho, xi_rho)
#           - 'threshold_90' : (dayofyear, [s_rho,] eta_rho, xi_rho)
#         and coordinate 'dayofyear' = 1..366.

#     Notes
#     -----
#     - This function never modifies your original monthly CROCO files.
#     - Only three files are opened per (year, calendar day), then closed promptly.
#     - Early January and late December windows automatically wrap across New Year:
#         * Jan-01 window includes Dec(Y-1) days
#         * Dec-31 window includes Jan(Y+1) days
#     - Feb-29 is skipped for non-leap years.
#     """

#     # Use a leap base year so all 366 calendar days exist for iteration
#     base_year = 2000
#     doy_start = datetime(base_year, 1, 1)

#     # We collect per-day results into Python lists and stack at the end.
#     # We'll allocate NaNs lazily using the shape from the first successful sample.
#     day_means = []
#     day_p90s = []
#     template_shape = None   # shape of one sample (so we can fill NaNs when needed)
#     is_3d = None            # True if var has vertical dimension (s_rho), else False

#     # Iterate over all calendar days 1..366 (covers leap day)
#     for doy in range(1, 367):
#         base_day = doy_start + timedelta(days=doy - 1)  # e.g., 2000-01-01, 2000-01-02, ...
#         base_month = base_day.month
#         base_dom   = base_day.day

#         if verbose:
#             print(f"[{doy:03d}/366] calendar day = {base_month:02d}-{base_dom:02d}")

#         # Container for all samples across all years for this one calendar day
#         samples = []

#         # Loop through each year once per calendar day
#         for Y in years:
#             # Some years won't have Feb-29; skip those gracefully
#             try:
#                 center = datetime(Y, base_month, base_dom)
#             except ValueError:
#                 # e.g., trying Feb-29 on a non-leap year
#                 continue

#             # Open exactly the three monthly files that can cover ±window around 'center'
#             opened = open_three_months_for_year(Y, base_month, data_dir)

#             # If none of the three exists, skip this year
#             if not opened:
#                 if verbose:
#                     print(f"  (Y={Y}) no monthly files found for months around {base_month:02d}")
#                 continue

#             # Define the ±window range around 'center'
#             start_win = center - timedelta(days=window)
#             end_win   = center + timedelta(days=window)

#             # For each day in this window, pick the already-opened month dataset and extract
#             for d in iter_days(start_win, end_win):
#                 ds_month = opened.get((d.year, d.month))  # exact (year, month)
#                 if ds_month is None:
#                     # window day falls in a month file we didn't/couldn't open (missing file)
#                     continue

#                 try:
#                     # Call crocotools get_var on the open Dataset.
#                     # It will slice to the nearest time for 'd' (using ref_date) and return an xarray.Dataset.
#                     ds_sel = post.get_var(ds_month, varname, time=d, ref_date=ref_date)
#                     arr = ds_sel[varname].values  # force load as a NumPy array now (no lazy graph)
#                     samples.append(arr)

#                     # On the very first successful sample, record the shape and whether it’s 3D
#                     if template_shape is None:
#                         template_shape = arr.shape
#                         is_3d = (arr.ndim == 3)
#                 except Exception:
#                     # If a particular timestamp is missing, just skip it; keep this lean and robust
#                     continue

#             # IMPORTANT: close the three open month files before moving to the next year
#             for ds in opened.values():
#                 try:
#                     ds.close()
#                 except Exception:
#                     pass

#         # After collecting all samples for this calendar day across all years:
#         if samples:
#             # Stack into a (N, ...) array and compute stats along axis 0
#             stack = np.stack(samples, axis=0).astype(dtype, copy=False)
#             mean_d = np.nanmean(stack, axis=0).astype(dtype, copy=False)
#             p90_d  = np.nanpercentile(stack, 90, axis=0).astype(dtype, copy=False)
#         else:
#             # No samples: fill with NaNs (use discovered template shape if we have one)
#             if template_shape is None:
#                 # Still no successful sample at all in the whole run so far
#                 mean_d = np.array(np.nan, dtype=dtype)
#                 p90_d  = np.array(np.nan, dtype=dtype)
#             else:
#                 mean_d = np.full(template_shape, np.nan, dtype=dtype)
#                 p90_d  = np.full(template_shape, np.nan, dtype=dtype)

#         # Store results for this calendar day
#         day_means.append(mean_d)
#         day_p90s.append(p90_d)

#         if verbose:
#             print(f"  collected samples: {len(samples):4d}")

#     # If we never got any data, stop here with a helpful error
#     if template_shape is None:
#         raise RuntimeError("No data collected—check paths, years, and filenames.")

#     # Build xarray outputs with correct dimensions:
#     # 3D variable -> (dayofyear, s_rho, eta_rho, xi_rho)
#     # 2D variable -> (dayofyear,         eta_rho, xi_rho)
#     if is_3d:
#         dims = ("dayofyear", "s_rho", "eta_rho", "xi_rho")
#     else:
#         dims = ("dayofyear", "eta_rho", "xi_rho")

#     # Convert lists to arrays of target dtype
#     arr_mean = np.array(day_means, dtype=dtype)
#     arr_p90  = np.array(day_p90s,  dtype=dtype)

#     da_mean = xr.DataArray(arr_mean, dims=dims, name=f"{varname}_climatology")
#     da_p90  = xr.DataArray(arr_p90,  dims=dims, name=f"{varname}_p90")

#     # Attach a 1..366 day-of-year coordinate
#     out = xr.Dataset(
#         {"climatology": da_mean, "threshold_90": da_p90},
#         coords={"dayofyear": np.arange(1, 367, dtype=int)},
#     )

#     # Save final product; this does NOT touch your original CROCO files
#     out.to_netcdf(savefile)
#     if verbose:
#         print(f"Saved climatology & 90th percentile to: {savefile}")

#     return out


# # ----------------- example usage -----------------
# if __name__ == "__main__":
#     # Full 1993–2019 span
#     years    = list(range(1993))
#     data_dir = "/media/nc.memela/Leviathan/MODELS/somisana-croco/configs/SA_West/croco_v1.3.1/C01_I01_GLORYS_ERA5/output"
#     varname  = "temp"
#     ref_date = datetime(2000, 1, 1)   # CROCO ref datetime (adjust if your setup differs)
#     window   = 15                     # ±15-day window (total 31 days)

#     # Run with verbose progress printing
#     ds_out = hobday_90th_percentile(
#         years, data_dir, varname, ref_date,
#         window=window, savefile="hobday_p90.nc",
#         verbose=True, dtype=np.float32,
#     )
