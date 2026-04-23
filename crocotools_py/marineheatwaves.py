import numpy as np
import pandas as pd
import xarray as xr
import scipy.ndimage as ndimage
from datetime import date
from netCDF4 import Dataset
import concurrent.futures
import gc


def detect_events_with_climatology(temp_data, clim_seas, clim_thresh, is_cold, t_dates=None):
    """
    Detect MHW/MCS events using pre-computed climatology.
 
    Parameters
    ----------
    temp_data   : array (time,)  — daily temperature time series
    clim_seas   : array (time,)  — seasonal climatology aligned to temp_data
    clim_thresh : array (time,)  — threshold (90th pct for MHW, 10th pct for MCS),
                                   also aligned to temp_data
    is_cold     : bool           — True for cold spells, False for heat waves
    t_dates     : array, optional — time vector in ordinal format
 
    Returns
    -------
    mhw        : dict   — event properties (same keys as Hobday detect output)
    categories : int8 array (time,)
                 0=none, 1=Moderate, 2=Strong, 3=Severe, 4=Extreme
                 (always positive; sign convention applied by the caller)
    """
    n_time     = len(temp_data)
    categories = np.zeros(n_time, dtype='int8')
 
    mhw = {
        'time_start': [], 'time_end': [], 'time_peak': [],
        'date_start': [], 'date_end': [], 'date_peak': [],
        'index_start': [], 'index_end': [], 'index_peak': [],
        'duration': [], 'duration_moderate': [], 'duration_strong': [],
        'duration_severe': [], 'duration_extreme': [],
        'intensity_max': [], 'intensity_mean': [], 'intensity_var': [],
        'intensity_cumulative': [],
        'intensity_max_relThresh': [], 'intensity_mean_relThresh': [],
        'intensity_var_relThresh': [], 'intensity_cumulative_relThresh': [],
        'intensity_max_abs': [], 'intensity_mean_abs': [],
        'intensity_var_abs': [], 'intensity_cumulative_abs': [],
        'category': [],
        'rate_onset': [], 'rate_decline': [],
        'n_events': 0,
    }
 
    if np.all(np.isnan(temp_data)) or np.all(temp_data == 0):
        return mhw, categories
 
    # Replace NaN with climatology so exceedance detection is clean
    temp_clean = temp_data.copy()
    temp_clean[np.isnan(temp_clean)] = clim_seas[np.isnan(temp_clean)]
 
    # Flip sign for cold spells so all subsequent logic is identical
    if is_cold:
        temp_clean     = -1.0 * temp_clean
        clim_seas_use  = -1.0 * clim_seas
        clim_thresh_use = -1.0 * clim_thresh
    else:
        clim_seas_use  = clim_seas
        clim_thresh_use = clim_thresh
 
    # Boolean exceedance array
    exceed_bool = (temp_clean - clim_thresh_use).copy()
    exceed_bool[exceed_bool <= 0]          = False
    exceed_bool[exceed_bool > 0]           = True
    exceed_bool[np.isnan(exceed_bool)]     = False

    # --- NEW: Bridge gaps of 1 or 2 days (Hobday criteria) ---
    true_indices = np.where(exceed_bool)[0]
    if len(true_indices) > 0:
        for i in range(len(true_indices) - 1):
            idx_current = true_indices[i]
            idx_next = true_indices[i+1]
            gap = idx_next - idx_current - 1
            if 1 <= gap <= 2:
                exceed_bool[idx_current+1:idx_next] = True
    # ---------------------------------------------------------
 
    events, n_events = ndimage.label(exceed_bool)
 
    cat_names = ['Moderate', 'Strong', 'Severe', 'Extreme']
    min_duration = 5
 
    for ev in range(1, n_events + 1):
        event_idx  = np.where(events == ev)[0]
        event_dur  = len(event_idx)
        if event_dur < min_duration:
            continue
 
        tt_start = event_idx[0]
        tt_end   = event_idx[-1]
 
        temp_mhw   = temp_clean[tt_start:tt_end + 1]
        thresh_mhw = clim_thresh_use[tt_start:tt_end + 1]
        seas_mhw   = clim_seas_use[tt_start:tt_end + 1]
 
        mhw_relSeas       = temp_mhw - seas_mhw
        mhw_relThresh     = temp_mhw - thresh_mhw
        mhw_relThreshNorm = (temp_mhw - thresh_mhw) / (thresh_mhw - seas_mhw)
        mhw_abs           = temp_mhw
 
        tt_peak = int(np.argmax(mhw_relSeas))
 
        mhw['index_start'].append(int(tt_start))
        mhw['index_end'].append(int(tt_end))
        mhw['index_peak'].append(int(tt_start + tt_peak))
 
        if t_dates is not None:
            mhw['time_start'].append(int(t_dates[tt_start]))
            mhw['time_end'].append(int(t_dates[tt_end]))
            mhw['time_peak'].append(int(t_dates[tt_start + tt_peak]))
            mhw['date_start'].append(date.fromordinal(int(t_dates[tt_start])))
            mhw['date_end'].append(date.fromordinal(int(t_dates[tt_end])))
            mhw['date_peak'].append(date.fromordinal(int(t_dates[tt_start + tt_peak])))
        else:
            for key in ('time_start', 'time_end', 'time_peak'):
                mhw[key].append(None)
            for key in ('date_start', 'date_end', 'date_peak'):
                mhw[key].append(None)
 
        mhw['duration'].append(event_dur)
 
        mhw['intensity_max'].append(float(mhw_relSeas[tt_peak]))
        mhw['intensity_mean'].append(float(mhw_relSeas.mean()))
        mhw['intensity_var'].append(float(np.sqrt(mhw_relSeas.var())))
        mhw['intensity_cumulative'].append(float(mhw_relSeas.sum()))
 
        mhw['intensity_max_relThresh'].append(float(mhw_relThresh[tt_peak]))
        mhw['intensity_mean_relThresh'].append(float(mhw_relThresh.mean()))
        mhw['intensity_var_relThresh'].append(float(np.sqrt(mhw_relThresh.var())))
        mhw['intensity_cumulative_relThresh'].append(float(mhw_relThresh.sum()))
 
        mhw['intensity_max_abs'].append(float(mhw_abs[tt_peak]))
        mhw['intensity_mean_abs'].append(float(mhw_abs.mean()))
        mhw['intensity_var_abs'].append(float(np.sqrt(mhw_abs.var())))
        mhw['intensity_cumulative_abs'].append(float(mhw_abs.sum()))
 
        tt_peakCat    = int(np.argmax(mhw_relThreshNorm))
        cats          = np.clip(np.floor(1.0 + mhw_relThreshNorm), 1, 5)
        peak_cat_val  = int(np.clip(cats[tt_peakCat], 1, 4))
        mhw['category'].append(cat_names[peak_cat_val - 1])
 
        mhw['duration_moderate'].append(int(np.sum(cats == 1)))
        mhw['duration_strong'].append(int(np.sum(cats == 2)))
        mhw['duration_severe'].append(int(np.sum(cats == 3)))
        mhw['duration_extreme'].append(int(np.sum(cats >= 4)))
 
        if event_dur > 1:
            mhw['rate_onset'].append(float(mhw_relSeas[tt_peak] / (tt_peak + 1)))
            mhw['rate_decline'].append(float(mhw_relSeas[tt_peak] / (event_dur - tt_peak)))
        else:
            mhw['rate_onset'].append(0.0)
            mhw['rate_decline'].append(0.0)
 
        categories[tt_start:tt_end + 1] = cats.astype('int8')
 
    mhw['n_events'] = len(mhw['duration'])
    return mhw, categories
 
 
def align_climatology_to_temp(temp_time, doy_values, clim_data):
    """
    Map a 366-element day-of-year climatology onto a temperature time series.
 
    Parameters
    ----------
    temp_time  : array of datetime64  — time axis of the temperature data
    doy_values : array (366,)         — day-of-year labels (1–366)
    clim_data  : array (366,)         — climatology values for each DOY
 
    Returns
    -------
    aligned_clim : float32 array (len(temp_time),)
    """
    clim_dict   = {int(d): clim_data[i] for i, d in enumerate(doy_values)}
    n_time      = len(temp_time)
    aligned     = np.empty(n_time, dtype='float32')
 
    for i, t in enumerate(temp_time):
        doy = int(pd.Timestamp(t).dayofyear)
        if doy in clim_dict:
            aligned[i] = clim_dict[doy]
        elif doy == 60 and 60 not in clim_dict:
            # Feb-29 in a non-leap year climatology: use Feb-28
            aligned[i] = clim_dict.get(59, np.nan)
        else:
            aligned[i] = np.nan
 
    return aligned
 
 
def process_level_batch(temp_slice, clim_seas_slice, clim_thresh_slice, is_cold, t_dates):
    """
    Detect MHW/MCS events for a (time, eta, xi) slab.
 
    Parameters
    ----------
    temp_slice        : float array (time, eta, xi)
    clim_seas_slice   : float array (time, eta, xi) — seas clim aligned to temp
    clim_thresh_slice : float array (time, eta, xi) — thresh aligned to temp
    is_cold           : bool
    t_dates           : ordinal int array (time,)
 
    Returns
    -------
    categories : int8 array (time, eta, xi)
    mhw_dicts  : list of dicts, one per pixel (None for all-NaN pixels)
    """
    n_time, n_eta, n_xi = temp_slice.shape
    categories = np.zeros((n_time, n_eta, n_xi), dtype='int8')
    mhw_dicts  = []
 
    for i in range(n_eta):
        for j in range(n_xi):
            temp_ts        = temp_slice[:, i, j]
            clim_seas_ts   = clim_seas_slice[:, i, j]
            clim_thresh_ts = clim_thresh_slice[:, i, j]
 
            if np.all(np.isnan(temp_ts)):
                mhw_dicts.append(None)
                continue
 
            mhw_ev, cats           = detect_events_with_climatology(
                temp_ts, clim_seas_ts, clim_thresh_ts, is_cold, t_dates
            )
            categories[:, i, j]   = cats
            mhw_dicts.append(mhw_ev)
 
    return categories, mhw_dicts
 
 
def create_mhw_output_netcdf(output_file, n_time_daily, n_levels, n_eta, n_xi,
                              ds_temp_daily, ds_temp, ds_clim, mode_name):
    """
    Create and initialise the output NetCDF file for MHW/MCS categories.
 
    Parameters
    ----------
    output_file   : str or Path
    n_time_daily  : int          — number of daily time steps
    n_levels      : int          — number of s_rho levels
    n_eta, n_xi   : int          — grid dimensions
    ds_temp_daily : xr.Dataset   — daily-averaged temperature (for time coord)
    ds_temp       : xr.Dataset   — original temperature file (for lat/lon)
    ds_clim       : xr.Dataset   — climatology file (fallback lat/lon source)
    mode_name     : str          — "MHW" or "MCS"
 
    Returns
    -------
    nc_out : netCDF4.Dataset  — open, ready for writing
    """
    nc_out = Dataset(output_file, mode='w', format='NETCDF4')
 
    nc_out.createDimension('time',    n_time_daily)
    nc_out.createDimension('s_rho',   n_levels)
    nc_out.createDimension('eta_rho', n_eta)
    nc_out.createDimension('xi_rho',  n_xi)
 
    time_var       = nc_out.createVariable('time', 'f8', ('time',))
    time_var.units    = 'days since 1993-01-01'
    time_var.calendar = 'standard'
    time_var[:]    = (
        (pd.to_datetime(ds_temp_daily.time.values) - pd.Timestamp('1993-01-01'))
        .total_seconds() / 86400
    )
 
    s_var    = nc_out.createVariable('s_rho', 'i4', ('s_rho',))
    s_var[:] = np.arange(n_levels, dtype='int32')
 
    for src in (ds_temp, ds_clim):
        if 'lat_rho' in src and 'lon_rho' in src:
            lv       = nc_out.createVariable('lat_rho', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
            lv[:]    = src['lat_rho'].values
            lov      = nc_out.createVariable('lon_rho', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
            lov[:]   = src['lon_rho'].values
            break
 
    cat_var = nc_out.createVariable(
        'category', 'i1',
        ('time', 's_rho', 'eta_rho', 'xi_rho'),
        zlib=True, complevel=4, fill_value=np.int8(0),
        chunksizes=(n_time_daily, 1, n_eta, n_xi),
    )
    cat_var.long_name     = f'{mode_name} event category'
    cat_var.flag_values   = np.array([-4, -3, -2, -1, 0, 1, 2, 3, 4], dtype='int8')
    cat_var.flag_meanings = (
        'extreme_MCS severe_MCS strong_MCS moderate_MCS '
        'no_event '
        'moderate_MHW strong_MHW severe_MHW extreme_MHW'
    )
    cat_var.description   = (
        'Positive = MHW (above 90th-pct threshold); '
        'Negative = MCS (below 10th-pct threshold). '
        'Magnitude: 1=Moderate, 2=Strong, 3=Severe, 4=Extreme.'
    )
 
    nc_out.title       = f'{mode_name} event detection'
    nc_out.description = 'Detected using pre-computed climatology'
    nc_out.created     = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
 
    return nc_out
 
 
def process_single_level(level, n_levels, ds_temp_daily, ds_clim, temp_var_name,
                          doy_values, is_cold, t_dates, batch_size, nc_out):
    """
    Detect MHW/MCS events for one s_rho level and write categories to nc_out.
 
    Parameters
    ----------
    level          : int         — s_rho index to process
    n_levels       : int         — total number of levels (for progress display)
    ds_temp_daily  : xr.Dataset  — daily temperature data
    ds_clim        : xr.Dataset  — climatology (day_of_year, s_rho, eta_rho, xi_rho)
                                   must contain 'climatology' and 'threshold_90'
                                   (or 'threshold_10' for MCS — selected via is_cold)
    temp_var_name  : str
    doy_values     : int array   — day-of-year axis of ds_clim (1–366)
    is_cold        : bool        — True → use threshold_10 and flip sign to negative
    t_dates        : ordinal int array (time,)
    batch_size     : int         — eta_rho rows per chunk
    nc_out         : netCDF4.Dataset — output file (must already have 'category' var)
 
    Notes
    -----
    - Reads climatology and threshold from ds_clim['climatology'] and either
      ds_clim['threshold_90'] (MHW) or ds_clim['threshold_10'] (MCS).
    - Aligns 366-element DOY arrays to the full temperature time axis via
      align_climatology_to_temp() for every pixel.
    - Writes int8 categories: +1…+4 for MHW, −1…−4 for MCS, 0 for no event.
    - Calls nc_out.sync() every 5 levels.
    """
    n_time_daily = ds_temp_daily.sizes['time']
    n_eta        = ds_temp_daily.sizes['eta_rho']
    n_xi         = ds_temp_daily.sizes['xi_rho']
    cat_var      = nc_out.variables['category']
 
    thresh_key = 'threshold_10' if is_cold else 'threshold_90'
 
    clim_seas_level   = ds_clim['climatology'].isel(s_rho=level).values   # (366, eta, xi)
    clim_thresh_level = ds_clim[thresh_key].isel(s_rho=level).values      # (366, eta, xi)
 
    print(f"\n   Level {level}/{n_levels - 1}:")
 
    for i in range(0, n_eta, batch_size):
        end_i = min(i + batch_size, n_eta)
 
        temp_slice        = ds_temp_daily[temp_var_name].isel(
            s_rho=level, eta_rho=slice(i, end_i)
        ).values                                          # (time, slab_eta, xi)
 
        clim_seas_slice   = np.empty((n_time_daily, end_i - i, n_xi), dtype='float32')
        clim_thresh_slice = np.empty((n_time_daily, end_i - i, n_xi), dtype='float32')
 
        temp_times = ds_temp_daily.time.values
        for li in range(end_i - i):
            for lj in range(n_xi):
                clim_seas_slice[:, li, lj] = align_climatology_to_temp(
                    temp_times, doy_values, clim_seas_level[:, i + li, lj]
                )
                clim_thresh_slice[:, li, lj] = align_climatology_to_temp(
                    temp_times, doy_values, clim_thresh_level[:, i + li, lj]
                )
 
        categories, _ = process_level_batch(
            temp_slice, clim_seas_slice, clim_thresh_slice, is_cold, t_dates
        )
 
        # Apply sign convention: MCS categories are negative
        if is_cold:
            categories = -np.abs(categories).astype('int8')
 
        cat_var[:, level, i:end_i, :] = categories
        print(f"      Rows {i:3d}-{end_i:3d} complete", end='\r')
 
    print(f"      Level {level} complete" + " " * 20)
 
    if level % 5 == 0:
        nc_out.sync()
 
    gc.collect()
 
 
def resample_to_daily(ds_temp, temp_var_name):
    """
    Resample temperature data to daily means and gap-fill short gaps (≤7 days).
    Assumes time coordinate is already datetime64 (e.g. as returned by get_var()).
 
    Parameters
    ----------
    ds_temp       : xr.Dataset
    temp_var_name : str
 
    Returns
    -------
    ds_temp_daily : xr.Dataset
    """
    ds_temp_daily = ds_temp.resample(time='1D').mean()
    ds_temp_daily = ds_temp_daily.interpolate_na(dim='time', limit=7)
    return ds_temp_daily
 
 