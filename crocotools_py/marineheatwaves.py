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

def load_and_harmonize_baselines(clim_file, thresh_file):
    """
    Centralized, memory-safe reader that safely extracts baseline datasets,
    drops heavy velocity arrays, and fixes integer vs float coordinate mismatch.
    """
    import xarray as xr
    print(f'Opening climatology (dropping heavy variables): {clim_file}')
    vars_to_drop = ['u', 'v', 'salt', 'ubar', 'vbar']
    ds_clim_raw = xr.open_dataset(clim_file, drop_variables=vars_to_drop)
    
    print(f'Opening thresholds: {thresh_file}')
    ds_thresh_raw = xr.open_dataset(thresh_file, drop_variables=vars_to_drop)

    # Unify dimensions
    if 'dayofyear' in ds_clim_raw.dims:
        ds_clim_raw = ds_clim_raw.rename_dims({'dayofyear': 'day_of_year'}).rename({'dayofyear': 'day_of_year'})
    if 'dayofyear' in ds_thresh_raw.dims:
        ds_thresh_raw = ds_thresh_raw.rename_dims({'dayofyear': 'day_of_year'}).rename({'dayofyear': 'day_of_year'})
    if 'temp' in ds_clim_raw.data_vars and 'climatology' not in ds_clim_raw.data_vars:
        ds_clim_raw = ds_clim_raw.rename({'temp': 'climatology'})

    # Enforce type harmonization to prevent alignment outer-join memory inflation
    for c in ['s_rho', 'eta_rho', 'xi_rho']:
        if c in ds_thresh_raw.coords and c in ds_clim_raw.coords:
            ds_thresh_raw = ds_thresh_raw.assign_coords({c: ds_clim_raw[c].values})

    # Assemble lightweight reference container
    ds_clim = xr.Dataset(coords=ds_clim_raw.coords)
    ds_clim['climatology'] = ds_clim_raw['climatology']
    if 'zeta' in ds_clim_raw.data_vars:
        ds_clim['zeta'] = ds_clim_raw['zeta']
    
    ds_clim['threshold_90'] = ds_thresh_raw['threshold_90']
    ds_clim['threshold_10'] = ds_thresh_raw['threshold_10']
    
    for v in ['lon_rho', 'lat_rho', 'day_of_year']:
        if v in ds_clim_raw.coords and v not in ds_clim.coords:
            ds_clim = ds_clim.assign_coords({v: ds_clim_raw[v]})
            
    return ds_clim


def detect_mhw_forecast(temp_file, clim_file, thresh_file, fname_out, temp_var='temp', Yorig=2000, batch_size=5):
    """
    Centralized operational tracking pipeline. Computes daily averages,
    MHW/MCS category bounds, temperature anomalies, and frontal gradients.
    """
    import os
    import gc
    import numpy as np
    import pandas as pd
    import xarray as xr
    from pathlib import Path
    from netCDF4 import Dataset
    from crocotools_py.postprocess import get_ds, handle_time

    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    out_file = Path(fname_out)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    # Stream clean baseline arrays
    ds_clim = load_and_harmonize_baselines(clim_file, thresh_file)

    num_levels = ds_clim.sizes['s_rho']
    n_eta      = ds_clim.sizes['eta_rho']
    n_xi       = ds_clim.sizes['xi_rho']
    doy_values = ds_clim['day_of_year'].values

    print(f'Loading temperature: {temp_file}')
    ds_temp = get_ds(temp_file, temp_var)
    ds_temp = handle_time(ds_temp, Yorig=Yorig)

    raw_times    = ds_temp.time.values
    first_day    = pd.Timestamp(raw_times[0]).normalize()
    last_day     = pd.Timestamp(raw_times[-1]).normalize()
    target_dates = pd.date_range(start=first_day, end=last_day, freq='1D')
    T_daily      = len(target_dates)
    t_dates      = np.array([d.toordinal() for d in target_dates], dtype=int)
    daily_doy_map = np.clip(target_dates.dayofyear.values - 1, 0, 365)

    if out_file.exists():
        out_file.unlink()

    nc_out = Dataset(str(out_file), mode='w', format='NETCDF4')
    nc_out.createDimension('xi_rho', n_xi)
    nc_out.createDimension('eta_rho', n_eta)
    nc_out.createDimension('s_rho', num_levels)
    nc_out.createDimension('time', T_daily)
    
    lon_var = nc_out.createVariable('lon_rho', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
    lat_var = nc_out.createVariable('lat_rho', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
    h_var   = nc_out.createVariable('h', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
    mask_v  = nc_out.createVariable('mask_rho', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
    
    lon_var[:] = ds_temp['lon_rho'].values if 'lon_rho' in ds_temp else ds_clim['lon_rho'].values
    lat_var[:] = ds_temp['lat_rho'].values if 'lat_rho' in ds_temp else ds_clim['lat_rho'].values
    h_var[:]   = ds_temp['h'].values
    mask_v[:]  = ds_temp['mask_rho'].values

    time_var = nc_out.createVariable('time', 'f8', ('time',))
    time_var.units = f'days since {first_day.strftime("%Y-%m-%d")}'
    time_var[:] = (target_dates - first_day).total_seconds() / 86400.0

    nc_out.createVariable('s_rho', 'f4', ('s_rho',))[:] = np.arange(num_levels)
    for attr in ['hc', 'Vtransform', 'theta_s', 'theta_b']:
        if attr in ds_temp:
            nc_out.setncattr(attr, float(ds_temp[attr].values))
        elif attr in ds_temp.attrs:
            nc_out.setncattr(attr, float(ds_temp.attrs[attr]))
    
    cat_var = nc_out.createVariable('category', 'i1', ('time', 's_rho', 'eta_rho', 'xi_rho'), 
                                     zlib=True, fill_value=-127, chunksizes=(T_daily, 1, n_eta, n_xi))
    cat_var.long_name = 'MHW_MCS Combined Event Categories'
    
    temp_anom_var = nc_out.createVariable('temp_anom', 'f4', ('time', 's_rho', 'eta_rho', 'xi_rho'), 
                                           zlib=True, fill_value=np.nan, chunksizes=(T_daily, 1, n_eta, n_xi))
    temp_anom_var.long_name = "Sea Water Temperature Daily Anomaly"
    temp_anom_var.units = "degC"
    
    zeta_var = nc_out.createVariable('zeta', 'f4', ('time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=np.nan)
    zeta_anom_var = nc_out.createVariable('zeta_anom', 'f4', ('time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=np.nan)

    sst_front_var = nc_out.createVariable('sst_front', 'f4', ('time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=np.nan)
    sst_front_var.long_name = "Sea Surface Temperature Horizontal Front Magnitude"
    sst_front_var.units = "degC / km"

    print("Computing Daily Zeta and Zeta Anomalies...")
    if 'zeta' in ds_temp:
        ds_zeta_daily = (ds_temp['zeta'].resample(time='1D').mean()
                         .reindex(time=target_dates, method='nearest')
                         .interpolate_na(dim='time', limit=None).compute())
        zeta_var[:] = ds_zeta_daily.values
        if 'zeta' in ds_clim:
            zeta_anom_var[:] = ds_zeta_daily.values - ds_clim['zeta'].values[daily_doy_map, :, :]

    print("Processing vertical planes...")
    try:
        for k in range(num_levels - 1, -1, -1):
            process_single_level(k, num_levels, ds_temp, ds_clim, temp_var, doy_values, False, t_dates, batch_size, nc_out)
            mhw_cats = cat_var[:, k, :, :][:]

            cat_var[:, k, :, :] = np.zeros_like(mhw_cats)
            process_single_level(k, num_levels, ds_temp, ds_clim, temp_var, doy_values, True, t_dates, batch_size, nc_out)
            mcs_cats = cat_var[:, k, :, :][:]

            combined = mhw_cats.copy()
            combined[combined == 0] = mcs_cats[combined == 0]
            
            mask_rho_2d = ds_temp['mask_rho'].values
            if mask_rho_2d.ndim > 2: mask_rho_2d = mask_rho_2d[0]
            combined[:, mask_rho_2d == 0] = -127
            cat_var[:, k, :, :] = combined

            print(f"      -> Calculating daily temperature anomalies for level {k}...")
            ds_level_daily = (ds_temp[temp_var].isel(s_rho=k).resample(time='1D').mean()
                              .reindex(time=target_dates, method='nearest')
                              .interpolate_na(dim='time', limit=None).compute())
            
            temp_anom_var[:, k, :, :] = ds_level_daily.values - ds_clim['climatology'].isel(s_rho=k).values[daily_doy_map, :, :]
            
            if k == num_levels - 1:
                print("      -> Calculating daily surface thermal fronts (SST fronts)...")
                pm = ds_temp['pm'].values if 'pm' in ds_temp else np.ones((n_eta, n_xi))
                pn = ds_temp['pn'].values if 'pn' in ds_temp else np.ones((n_eta, n_xi))
                if pm.ndim > 2: pm, pn = pm[0], pm[0]
                
                sst_front_data = np.zeros((T_daily, n_eta, n_xi), dtype='float32')
                sst_vals = ds_level_daily.values
                for t_idx in range(T_daily):
                    d_eta, d_xi = np.gradient(sst_vals[t_idx])
                    sst_front_data[t_idx] = np.hypot(d_xi * pm, d_eta * pn) * 1000.0
                    
                sst_front_var[:] = np.where(mask_rho_2d[np.newaxis, :, :] == 1, sst_front_data, np.nan)

            nc_out.sync()
            gc.collect()
    finally:
        nc_out.sync()
        nc_out.close()

    ds_clim.close()
    print(f'Done: {fname_out} ({out_file.stat().st_size / (1024 ** 2):.1f} MB)')
    
# =====================================================================
# OPERATIONAL PLOTTING & ANIMATION MODULE
# =====================================================================

import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.patheffects as pe
from matplotlib.animation import FuncAnimation
from matplotlib.patches import FancyBboxPatch, Wedge
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import crocotools_py.postprocess as post

TARGETS = {
    "Kleinsee":       (17.030382, -29.680623), "Hondeklipbaai":  (17.252461, -30.315292),
    "Doringbaai":     (18.213554, -31.814509), "Elandsbaai":     (18.30165,  -32.312317),
    "Laaiplek":       (18.125354, -32.742041), "Paternoster":    (17.870305, -32.777566),
    "Saldanha":       (17.929861, -33.074807), "Yzerfontein":    (18.13382,  -33.361876),
    "Bloubergstrand": (18.443896, -33.803906), "Oudekraal":      (18.342541, -33.980098),
    "Cape Point":     (18.46024,  -34.358313), "Simonstown":     (18.442294, -34.176514),
    "Strand":         (18.810174, -34.120553), "Hangklip":       (18.803882, -34.374716),
    "Kleinmond":      (19.026591, -34.355882), "Hermanus":       (19.256989, -34.425957),
    "Gansbaai":       (19.323381, -34.576985),
}

WINDOW_DAYS = 10
FILL_MOD   = "#ffc73e";  FILL_STR   = "#f77819"; FILL_SEV   = "#bf460c";  FILL_EXT   = "#4e1909"
FILL_C_MOD = "#a6d3e8";  FILL_C_STR = "#5da6c9"; FILL_C_SEV = "#2074a3";  FILL_C_EXT = "#103c68"

MHW_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_MOD,   2: FILL_STR,   3: FILL_SEV,   4: FILL_EXT}
MCS_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_C_MOD, 2: FILL_C_STR, 3: FILL_C_SEV, 4: FILL_C_EXT}

CMAP_9 = mplc.ListedColormap([FILL_C_EXT, FILL_C_SEV, FILL_C_STR, FILL_C_MOD, "#ffffff", FILL_MOD, FILL_STR, FILL_SEV, FILL_EXT])
CMAP_9.set_bad("white")
BNORM_9 = mplc.BoundaryNorm([-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5], CMAP_9.N)

def nearest(lon2d, lat2d, lon0, lat0):
    d2 = (lon2d - lon0)**2 + (lat2d - lat0)**2
    return np.unravel_index(np.argmin(d2), d2.shape)

def doy_index(times):
    return np.clip(pd.to_datetime(times).dayofyear.values - 1, 0, 365)

def compute_site_flag_data(sites, cat_ds, lev):
    site_data = {}
    for site_name, data in sites.items():
        pj, pi = data["pj"], data["pi"]
        cat = (cat_ds["category"].isel(s_rho=lev, eta_rho=pj, xi_rho=pi).load().values.astype(float))
        cat[cat == -127] = 0
        mhw_days = cat[cat > 0]
        mcs_days = cat[cat < 0]
        max_mhw = float(np.max(mhw_days)) if len(mhw_days) > 0 else 0.0
        max_mcs = float(np.abs(np.min(mcs_days))) if len(mcs_days) > 0 else 0.0

        if max_mhw == 0.0 and max_mcs == 0.0:
            site_data[site_name] = {"mode": "None", "max_cat": 0}
        elif max_mhw >= max_mcs:
            site_data[site_name] = {"mode": "MHW", "max_cat": max_mhw}
        else:
            site_data[site_name] = {"mode": "MCS", "max_cat": max_mcs}
    return site_data

def plot_timeseries_multisite(sites, today, output_dir, depth_name):
    out_dir = Path(output_dir); out_dir.mkdir(parents=True, exist_ok=True)
    today = pd.Timestamp(today)

    for site_name, data in sites.items():
        fct_dates = pd.to_datetime(data["fct_dates"])
        fct_temp  = np.atleast_1d(data["fct_temp"])
        fct_seas  = np.atleast_1d(data["fct_seas"])
        fct_h_thr = np.atleast_1d(data["fct_h_thr"])
        fct_c_thr = np.atleast_1d(data["fct_c_thr"])

        obs_dates = pd.to_datetime(data["obs_dates"])
        obs_temp  = np.atleast_1d(data["obs_temp"])
        obs_seas  = np.atleast_1d(data["obs_seas"])
        obs_h_thr = np.atleast_1d(data["obs_h_thr"])
        obs_c_thr = np.atleast_1d(data["obs_c_thr"])

        all_dates = np.concatenate([obs_dates, fct_dates])
        all_temp  = np.concatenate([obs_temp, fct_temp])
        all_seas  = np.concatenate([obs_seas, fct_seas])
        all_h_thr = np.concatenate([obs_h_thr, fct_h_thr])
        all_c_thr = np.concatenate([obs_c_thr, fct_c_thr])

        fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
        ax.yaxis.grid(True, color="#cccccc", linewidth=0.7, zorder=0)
        ax.set_facecolor("white"); fig.patch.set_facecolor("white")

        ax.plot(all_dates, all_seas, ":", color="gray", label="Climatology", lw=1.5, zorder=2)
        ax.plot(all_dates, all_h_thr, "--", color="#d9534f", label="MHW threshold", lw=1.2, zorder=2)
        ax.plot(all_dates, all_c_thr, "--", color="#337ab7", label="MCS threshold", lw=1.2, zorder=2)

        if len(obs_dates) > 0:
            ax.plot(obs_dates, obs_temp, color="#777777", lw=2.5, label="SST observed", zorder=5)
        if len(fct_dates) > 0:
            ax.plot(fct_dates, fct_temp, color="black", lw=2.5, label="SST forecast", zorder=5)

        ax.axvline(today, color="black", lw=1.2, zorder=6)
        ax.text(today + pd.Timedelta(hours=4), ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1]-ax.get_ylim()[0]), "Today", va="top", ha="left", fontsize=9, fontweight="bold")

        ax.set_title(f"{site_name}  ({abs(data['lat']):.3f}°S, {data['lon']:.3f}°E)", fontsize=14, fontweight="bold", pad=10, color="#1a3a5c")
        ax.set_ylabel(f"Temperature [°C]", fontsize=11, fontweight="bold", color="#1a3a5c")
        ax.set_xlim(all_dates[0], all_dates[-1])
        
        ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d'))
        for spine in ("top", "right"): ax.spines[spine].set_visible(False)
        ax.spines['left'].set_color('#1a3a5c'); ax.spines['bottom'].set_color('#1a3a5c')

        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=9, frameon=False)
        plt.savefig(out_dir / f"{site_name.replace(' ', '_')}_{depth_name}_{today.strftime('%Y%m%d')}.png", dpi=150, bbox_inches="tight")
        plt.close()

def plot_flag_map(site_data, today, start_date, end_date, out_path, lat, lon, depth_name="Surface"):
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    
    def _flag_col(mode, cat):
        if mode == "None" or cat == 0: return "#4CAF7D"
        return (MHW_FLAG_COLOURS if mode == "MHW" else MCS_FLAG_COLOURS)[max(0, min(4, int(round(cat))))]

    def _draw_gauge(ax_g):
        cat_labels = {1: "Mod", 2: "Str", 3: "Sev", 4: "Ext"}
        ax_g.set_xlim(-1.55, 1.55); ax_g.set_ylim(-1.55, 1.55); ax_g.set_aspect("equal"); ax_g.axis("off")
        n, r_out, r_in = 4, 1.30, 0.48

        for k, (th1, th2) in enumerate(zip(np.degrees(np.linspace(0, np.pi, n + 1)), np.degrees(np.linspace(0, np.pi, n + 1))[1:])):
            ax_g.add_patch(Wedge((0, 0), r_out, th1, th2, width=r_out - r_in, fc=MHW_FLAG_COLOURS[k+1], ec="white", lw=0.8, zorder=1))
            mid = np.radians((th1 + th2) / 2); rl = (r_out + r_in) / 2
            ax_g.text(rl * np.cos(mid), rl * np.sin(mid), cat_labels[k+1], ha="center", va="center", fontsize=7.0, fontweight="bold", color="white", rotation=np.degrees(mid) - 90, zorder=3)

        for k, (th1, th2) in enumerate(zip(np.degrees(np.linspace(np.pi, 2 * np.pi, n + 1)), np.degrees(np.linspace(np.pi, 2 * np.pi, n + 1))[1:])):
            ax_g.add_patch(Wedge((0, 0), r_out, th1, th2, width=r_out - r_in, fc=MCS_FLAG_COLOURS[n-k], ec="white", lw=0.8, zorder=1))
            mid = np.radians((th1 + th2) / 2); rl = (r_out + r_in) / 2
            ax_g.text(rl * np.cos(mid), rl * np.sin(mid), cat_labels[n-k], ha="center", va="center", fontsize=7.0, fontweight="bold", color="white", rotation=np.degrees(mid) - 90, zorder=3)

        ax_g.add_patch(plt.Circle((0, 0), r_in, fc=MHW_FLAG_COLOURS[0], ec="white", lw=1.0, zorder=2))
        ax_g.text(0, 0, "None", ha="center", va="center", fontsize=8, fontweight="bold", color="white", zorder=4)
        ax_g.text(0,  r_out + 0.10, "MHW", ha="center", va="bottom", fontsize=7.5, fontweight="bold", color=MHW_FLAG_COLOURS[2])
        ax_g.text(0, -(r_out + 0.10), "MCS", ha="center", va="top", fontsize=7.5, fontweight="bold", color=MCS_FLAG_COLOURS[3])
        ax_g.set_title("Max Intensity\n(Discrete Flags)", fontsize=7, fontweight="bold", pad=3, color="#1a3a5c")

    coast_order = ["Kleinsee", "Hondeklipbaai", "Doringbaai", "Elandsbaai", "Laaiplek", "Paternoster", "Saldanha", "Yzerfontein", "Bloubergstrand", "Oudekraal", "Cape Point", "Simonstown", "Strand", "Hangklip", "Kleinmond", "Hermanus", "Gansbaai"]
    BOX_SIZE, OFFSHORE, BOX_STEP_DIST = 0.45, -0.10, 0.50
    all_boxes = []
    
    dense_lons, dense_lats = [], []
    for k in range(len(coast_order) - 1):
        lon0, lat0 = TARGETS[coast_order[k]]; lon1, lat1 = TARGETS[coast_order[k + 1]]
        ts = np.linspace(0, 1, 100)
        dense_lons.extend(lon0 + ts * (lon1 - lon0)); dense_lats.extend(lat0 + ts * (lat1 - lat0))
        
    dense_lons, dense_lats = np.array(dense_lons), np.array(dense_lats)
    dists = np.zeros(len(dense_lons)); dists[1:] = np.cumsum(np.hypot(np.diff(dense_lons), np.diff(dense_lats)))
    
    for bd in np.arange(0, dists[-1], BOX_STEP_DIST):
        cx, cy = np.interp(bd, dists, dense_lons), np.interp(bd, dists, dense_lats)
        nearest_site = min(coast_order, key=lambda s: np.hypot(cx - TARGETS[s][0], cy - TARGETS[s][1]))
        info = site_data.get(nearest_site, {"mode": "None", "max_cat": 0})
        idx = max(1, min(np.searchsorted(dists, bd), len(dists) - 1))
        slen = np.hypot(dense_lons[idx] - dense_lons[idx-1], dense_lats[idx] - dense_lats[idx-1]) or 1.0
        px, py = -(dense_lats[idx] - dense_lats[idx-1]) / slen, (dense_lons[idx] - dense_lons[idx-1]) / slen
        all_boxes.append((cx + px * OFFSHORE, cy + py * OFFSHORE, _flag_col(info["mode"], info["max_cat"])))

    fig = plt.figure(figsize=(10, 13), dpi=150); ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.set_extent([15.8, 20.5, -36.0, -28.0], crs=ccrs.PlateCarree())

    for cx, cy, col in all_boxes:
        ax.add_patch(FancyBboxPatch((cx - BOX_SIZE/2, cy - BOX_SIZE/2), BOX_SIZE, BOX_SIZE, boxstyle="round,pad=0.04", facecolor=col, edgecolor="white", linewidth=0.6, zorder=3, transform=ccrs.PlateCarree()))

    for site_name, (site_lon, site_lat) in TARGETS.items():
        info = site_data.get(site_name, {"mode": "None", "max_cat": 0})
        cat_int = max(0, min(4, int(round(info["max_cat"]))))
        ax.plot(site_lon, site_lat, "o", ms=4, color="white", zorder=8, mec="black", mew=0.8, transform=ccrs.PlateCarree())
        lbl_txt = f"{site_name}\nNone" if info["mode"] == "None" or cat_int == 0 else f"{site_name}\n{info['mode']} – {['None', 'Moderate', 'Strong', 'Severe', 'Extreme'][cat_int]}"
        ax.text(site_lon + 0.08, site_lat, lbl_txt, ha="left", va="center", fontsize=6.5, fontweight="bold", color="#1a3a5c", zorder=9, transform=ccrs.PlateCarree(), path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=4); ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor="#555544", zorder=5)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="#aaaaaa", alpha=0.8, linestyle="--", zorder=2); gl.top_labels = gl.right_labels = False
    _draw_gauge(fig.add_axes([0.53, 0.61, 0.28, 0.28]))
    ax.set_title(f"SA West Coast  ·  MHW / MCS Flag Map  ·  {depth_name}\nForecast: {pd.to_datetime(start_date).strftime('%d %b')} – {pd.to_datetime(end_date).strftime('%d %b %Y')}", fontsize=12, color="#1a3a5c", pad=8)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()

def _update_spatial_frame(frame, cat_data, time_data, mesh_obj, title_obj, d_name):
    mesh_obj.set_array(cat_data[frame].ravel())
    title_obj.set_text(f"MHW & MCS Categories ({d_name})\nDate: {str(time_data[frame])[:10]}")
    return mesh_obj, title_obj

def animate_spatial_categories(cat_ds, ds_fcst, lat, lon, depth_name, lev, is_varying, idx_2d, out_path):
    out_path = Path(out_path); times = pd.to_datetime(cat_ds.time.values)
    cat = cat_ds["category"].values[:, idx_2d, np.meshgrid(np.arange(idx_2d.shape[0]), np.arange(idx_2d.shape[1]), indexing="ij")[1], np.meshgrid(np.arange(idx_2d.shape[0]), np.arange(idx_2d.shape[1]), indexing="ij")[0]].astype(float) if is_varying else cat_ds["category"].isel(s_rho=lev).values.astype(float)
    cat[cat == -127] = np.nan
    mask = ds_fcst["mask_rho"].values if "mask_rho" in ds_fcst else np.ones_like(lat)
    if mask.ndim > 2: mask = mask[0]
    cat = np.where(mask[np.newaxis, :, :] == 1, cat, np.nan)

    fig = plt.figure(figsize=(9, 8)); ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(lon, lat, cat[0], transform=ccrs.PlateCarree(), cmap=CMAP_9, norm=BNORM_9, shading="auto")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6); ax.add_feature(cfeature.LAND, facecolor="white", edgecolor='black', zorder=2); ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=":")
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="#aaaaaa", alpha=0.8, linestyle="--", zorder=2); gl.top_labels = gl.right_labels = False

    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.04, ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4])
    cbar.ax.set_yticklabels(["Ext", "Sev", "Str", "Mod", "Neut", "Mod", "Str", "Sev", "Ext"])
    cbar.set_label("MCS (Cold)  ←  Intensity  →  MHW (Heat)")

    title = ax.set_title(f"MHW & MCS Categories ({depth_name})\nDate: {str(times[0])[:10]}")
    ani = FuncAnimation(fig, _update_spatial_frame, frames=len(times), fargs=(cat, times, mesh, title, depth_name), blit=False)
    ani.save(out_path, writer='ffmpeg', fps=5, dpi=120); plt.close(fig)
    
def _update_anomaly_frame(frame, anom_data, time_data, mesh_obj, title_obj):
    mesh_obj.set_array(anom_data[frame].ravel())
    title_obj.set_text(f"Sea Water Temperature Daily Anomaly (Surface)\nDate: {str(time_data[frame])[:10]}")
    return mesh_obj, title_obj

def animate_surface_anomalies(cat_ds, lat, lon, out_path):
    out_path = Path(out_path); times = pd.to_datetime(cat_ds.time.values)
    anom = cat_ds["temp_anom"].isel(s_rho=-1).values.astype(float)

    fig = plt.figure(figsize=(9, 8)); ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(lon, lat, anom[0], transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=-3.0, vmax=3.0, shading="auto")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6); ax.add_feature(cfeature.LAND, facecolor="lightgray", edgecolor='black', zorder=2)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="#aaaaaa", alpha=0.8, linestyle="--", zorder=2); gl.top_labels = gl.right_labels = False
    
    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.04); cbar.set_label("Temperature Anomaly (°C)")
    title = ax.set_title(f"Sea Water Temperature Daily Anomaly (Surface)\nDate: {str(times[0])[:10]}")
    ani = FuncAnimation(fig, _update_anomaly_frame, frames=len(times), fargs=(anom, times, mesh, title), blit=False)
    ani.save(out_path, writer='ffmpeg', fps=5, dpi=120); plt.close(fig)

def _update_frontal_frame(frame, front_data, time_data, mesh_obj, title_obj):
    mesh_obj.set_array(front_data[frame].ravel())
    title_obj.set_text(f"Horizontal Thermal Front Magnitude (Surface)\nDate: {str(time_data[frame])[:10]}")
    return mesh_obj, title_obj

def animate_surface_fronts(cat_ds, lat, lon, out_path):
    out_path = Path(out_path); times = pd.to_datetime(cat_ds.time.values)
    front = cat_ds["sst_front"].values.astype(float)

    fig = plt.figure(figsize=(9, 8)); ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(lon, lat, front[0], transform=ccrs.PlateCarree(), cmap='inferno', vmin=0.05, vmax=0.50, shading="auto")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6); ax.add_feature(cfeature.LAND, facecolor="lightgray", edgecolor='black', zorder=2)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="#aaaaaa", alpha=0.8, linestyle="--", zorder=2); gl.top_labels = gl.right_labels = False
    
    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.04); cbar.set_label("SST Gradient Magnitude (°C / km)")
    title = ax.set_title(f"Horizontal Thermal Front Magnitude (Surface)\nDate: {str(times[0])[:10]}")
    ani = FuncAnimation(fig, _update_frontal_frame, frames=len(times), fargs=(front, times, mesh, title), blit=False)
    ani.save(out_path, writer='ffmpeg', fps=5, dpi=120); plt.close(fig)


def plot_operational_mhw_mcs(forecast_file, cat_file, clim_file, thresh_file, out_dir, start_date, end_date, Yorig=2000):
    print("Rendering Operational MHW/MCS Visuals")
    out_dir = Path(out_dir)
    
    ds_clim = load_and_harmonize_baselines(clim_file, thresh_file)
    ds_cat  = xr.open_dataset(cat_file)
    ds_fcst = post.handle_time(post.get_ds(forecast_file, "temp"), Yorig=Yorig)
    
    lat = ds_fcst.lat_rho.values if "lat_rho" in ds_fcst else ds_fcst.lat.values
    if lat.ndim > 2: lat, lon = lat[0], ds_fcst.lon_rho.values[0] if "lon_rho" in ds_fcst else ds_fcst.lon.values[0]
    else: lon = ds_fcst.lon_rho.values if "lon_rho" in ds_fcst else ds_fcst.lon.values
    
    h = ds_fcst.h.values if "h" in ds_fcst else np.zeros_like(lat)
    if h.ndim > 2: h = h[0]
    nlev = len(ds_fcst.s_rho) if "s_rho" in ds_fcst else ds_fcst.dims.get("s_rho", 32)
    today = pd.Timestamp(ds_fcst.time.values[4]).normalize()

    depth_levels = {
        "Surface": {"type": "fixed", "lev": nlev - 1},
        "Bottom":  {"type": "fixed", "lev": 0},
    }

    ds_fcst_single = post.handle_time(xr.open_dataset(forecast_file, decode_times=False), Yorig=Yorig)

    for depth_name, depth_info in depth_levels.items():
        print(f"\nProcessing Depth: {depth_name}...")
        sites = {}
        for site_name, (site_lon, site_lat) in TARGETS.items():
            pj, pi = nearest(lon, lat, site_lon, site_lat)
            lev_site = depth_info["lev"]
            ts = ds_fcst_single["temp"].isel(s_rho=lev_site).resample(time="1D").mean().load()
            all_dates, all_temps = pd.to_datetime(ts.time.values), ts.isel(eta_rho=pj, xi_rho=pi).values
            doy_all = doy_index(all_dates)
            obs_m, fct_m = all_dates < today, all_dates >= today
            
            sites[site_name] = dict(
                pj=int(pj), pi=int(pi), lon=float(lon[pj, pi]), lat=float(lat[pj, pi]),
                obs_dates=pd.DatetimeIndex(all_dates[obs_m]), obs_temp=all_temps[obs_m],
                obs_seas=ds_clim["climatology"].isel(s_rho=lev_site, eta_rho=pj, xi_rho=pi).values[doy_all][obs_m],
                obs_h_thr=ds_clim["threshold_90"].isel(s_rho=lev_site, eta_rho=pj, xi_rho=pi).values[doy_all][obs_m],
                obs_c_thr=ds_clim["threshold_10"].isel(s_rho=lev_site, eta_rho=pj, xi_rho=pi).values[doy_all][obs_m],
                fct_dates=pd.DatetimeIndex(all_dates[fct_m]), fct_temp=all_temps[fct_m],
                fct_seas=ds_clim["climatology"].isel(s_rho=lev_site, eta_rho=pj, xi_rho=pi).values[doy_all][fct_m],
                fct_h_thr=ds_clim["threshold_90"].isel(s_rho=lev_site, eta_rho=pj, xi_rho=pi).values[doy_all][fct_m],
                fct_c_thr=ds_clim["threshold_10"].isel(s_rho=lev_site, eta_rho=pj, xi_rho=pi).values[doy_all][fct_m],
            )

        print("  -> Time Series...")
        plot_timeseries_multisite(sites, today, out_dir / depth_name, depth_name)

        print("  -> Flag Maps...")
        plot_flag_map(compute_site_flag_data(sites, ds_cat, depth_info["lev"]), today, start_date, end_date, out_dir / f"FlagMap_{depth_name}_{today.strftime('%Y%m%d')}.png", lat, lon, depth_name)

        print("  -> Spatial Category Animation...")
        animate_spatial_categories(ds_cat, ds_fcst, lat, lon, depth_name, depth_info["lev"], False, None, out_dir / f"Categories_Animation_{depth_name}.mp4")

        if depth_name == "Surface":
            print("  -> Temperature Anomaly Animation...")
            animate_surface_anomalies(ds_cat, lat, lon, out_dir / "Temperature_Anomaly_Animation_Surface.mp4")
            print("  -> Thermal Front Animation...")
            animate_surface_fronts(ds_cat, lat, lon, out_dir / "Thermal_Front_Animation_Surface.mp4")

    ds_fcst_single.close(); ds_fcst.close(); ds_clim.close(); ds_cat.close()
    print(f"\nAll operational visuals saved cleanly to: {out_dir}")
 
 