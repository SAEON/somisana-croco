"""
Marine Heatwave (MHW) and Marine Cold Spell (MCS) detection.

Event detection against a pre-built day-of-year climatology and percentile
thresholds, following Hobday et al. (2016). This module holds the detection
algorithm only, with no operational plumbing or plotting, so that it can be
imported on its own for hindcast analysis.

The operational forecast products built on top of it (the daily-averaged
products file, anomalies, SST fronts, stratification, and the figures and
animations) live in crocotools_py.products.
"""

from datetime import date
import numpy as np
import pandas as pd
import xarray as xr
import scipy.ndimage as ndimage

# Core Heatwave Algorithms
def empty_event_dict():
    """
    An empty per-event statistics dict, as returned by
    detect_events_with_climatology(). Kept as its own function so the same
    definition is used whether or not statistics are being collected.
    """
    return {'time_start': [], 'time_end': [], 'time_peak': [],
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
        'n_events': 0,}


def detect_events_with_climatology(temp_data, clim_seas, clim_thresh, is_cold, t_dates=None,
                                    stats=True):
    """
    Detect MHW/MCS events using pre-computed climatology.

    Returns (mhw, categories) where 'categories' is the per-timestep signed
    category time-series and 'mhw' is a dict of per-event statistics
    (durations, intensities, onset/decline rates etc - see
    empty_event_dict()).

    Categories follow Hobday et al. (2016): with dT = threshold - climatology,
    a day is category floor((T - climatology)/dT), i.e. 1 = Moderate (T is
    between 1x and 2x dT above the climatology), 2 = Strong, 3 = Severe and
    4 = Extreme. Category 4 is an open-ended top bucket (>= 4x dT) - there is
    no category 5 in the scheme.

    Caveat when applying this through the water column: the category is a
    multiple of dT, so it assumes dT is physically meaningful. That is safe for
    SST, which always has a substantial seasonal cycle, but in deep water with
    almost no seasonal signal dT collapses towards zero (values below 0.001
    degC occur on the sa-west grid) and a trivial temperature wobble normalises
    up into a high category. On a sa-west forecast this affects of order 0.1%
    of ocean cell-levels, essentially all of them in the lower water column, so
    it is left in rather than special-cased - but categories at depth in
    low-variability cells should be read with that in mind.

    stats : if True (the default) the per-event statistics are computed, as
            they always have been - this is what you want for hindcast
            analysis. If False, only 'categories' is computed and the first
            return value is None. Building those statistics is the bulk of
            the cost of this function, and the operational forecast pipeline
            writes only the categories, so it passes stats=False. The
            'categories' output is identical either way - the statistics are
            computed in an 'if stats:' block *after* the category assignment,
            so both settings run the same category code path.
    """
    n_time     = len(temp_data)
    categories = np.zeros(n_time, dtype='int8')

    mhw = empty_event_dict() if stats else None

    if np.all(np.isnan(temp_data)) or np.all(temp_data == 0):
        return mhw, categories

    temp_clean = temp_data.copy()
    temp_clean[np.isnan(temp_clean)] = clim_seas[np.isnan(temp_clean)]

    if is_cold:
        exceed_bool = (temp_clean < clim_thresh).astype(float)
    else:
        exceed_bool = (temp_clean > clim_thresh).astype(float)
    exceed_bool[np.isnan(temp_clean) | np.isnan(clim_thresh)] = 0.0

    exceed_mask = exceed_bool.astype(bool)
    true_indices = np.where(exceed_mask)[0]
    if len(true_indices) > 0:
        for i in range(len(true_indices) - 1):
            gap = true_indices[i+1] - true_indices[i] - 1
            if 1 <= gap <= 2:
                exceed_mask[true_indices[i]+1 : true_indices[i+1]] = True

    events, n_events = ndimage.label(exceed_mask) 

    if is_cold:
        clim_seas_use   = -1.0 * clim_seas
        clim_thresh_use = -1.0 * clim_thresh
        temp_work       = -1.0 * temp_clean
    else:
        clim_seas_use   = clim_seas
        clim_thresh_use = clim_thresh
        temp_work       = temp_clean

    cat_names = ['Moderate', 'Strong', 'Severe', 'Extreme']
    min_duration = 5

    for ev in range(1, n_events + 1):
        event_idx  = np.where(events == ev)[0]
        event_dur  = len(event_idx)
        if event_dur < min_duration:
            continue

        tt_start = event_idx[0]
        tt_end   = event_idx[-1]

        temp_mhw = temp_work[tt_start:tt_end + 1]
        thresh_mhw = clim_thresh_use[tt_start:tt_end + 1]
        seas_mhw   = clim_seas_use[tt_start:tt_end + 1]

        # The category time-series is the only thing the operational pipeline
        # needs, so it is computed (and assigned) before the statistics block
        # below - that way stats=True and stats=False run identical code to
        # get here and are guaranteed to produce the same categories.
        #
        # floor(1 + (T - thresh)/(thresh - clim)) == floor((T - clim)/dT), the
        # Hobday category. It is capped at 4 because Extreme is an open-ended
        # top bucket (>= 4x dT) - which is how duration_extreme counts it
        # below (cats >= 4) and how peak_cat_val names it. Capping at 5 instead
        # let a category with no name, no colour in the plotting colormap and
        # no place in the file's own metadata reach the output.
        mhw_relThreshNorm = (temp_mhw - thresh_mhw) / (thresh_mhw - seas_mhw)
        cats              = np.clip(np.floor(1.0 + mhw_relThreshNorm), 1, 4)
        categories[tt_start:tt_end + 1] = cats.astype('int8')

        if not stats:
            continue

        mhw_relSeas       = temp_mhw - seas_mhw
        mhw_relThresh     = temp_mhw - thresh_mhw
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
            for key in ('time_start', 'time_end', 'time_peak'): mhw[key].append(None)
            for key in ('date_start', 'date_end', 'date_peak'): mhw[key].append(None)

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

    if stats:
        mhw['n_events'] = len(mhw['duration'])
    return mhw, categories

def build_doy_alignment_index(temp_time, doy_values):
    """
    Precompute, once, the mapping from each timestep in temp_time to its
    position along the day-of-year axis of a climatology array (length
    len(doy_values)). Used by load_daily_baselines() to pull out only the
    days we need in a single gather - the mapping is the same for every
    level and every grid cell, so it only ever needs computing once.

    Returns
    -------
    idx   : int array, shape (n_time,) - position into the climatology
            array's day-of-year axis for each timestep (0 where invalid -
            see `valid`)
    valid : bool array, shape (n_time,) - False where no climatology value
            exists for that day (caller should fill those rows with NaN)
    """
    doy_pos = {int(d): i for i, d in enumerate(doy_values)}
    n_time = len(temp_time)
    idx = np.zeros(n_time, dtype='int64')
    valid = np.zeros(n_time, dtype=bool)

    for i, t in enumerate(temp_time):
        doy = int(pd.Timestamp(t).dayofyear)
        if doy in doy_pos:
            idx[i] = doy_pos[doy]
            valid[i] = True
        elif doy == 60 and 60 not in doy_pos and 59 in doy_pos:
            idx[i] = doy_pos[59]
            valid[i] = True
        # else: idx stays 0, valid stays False -> caller fills with NaN

    return idx, valid


def load_daily_baselines(ds_clim, temp_time):
    """
    Read the climatology and the two percentile thresholds ONCE, for only the
    days we actually need, across all vertical levels.

    The day-of-year climatology files are chunked with s_rho spanning a full
    chunk (e.g. 30 dayofyear x 30 s_rho x 5 eta_rho x 130 xi_rho, deflate 4),
    so pulling out a single level with .isel(s_rho=k) has to inflate the whole
    ~10 GB file - about 20-30 s per read, and that is decompression cost, not
    disk, so the OS page cache doesn't help. Doing that once per level per
    MHW/MCS pass was costing over an hour per run. A forecast only spans a
    handful of days, so we gather those days up front instead and keep the
    (n_time, s_rho, eta_rho, xi_rho) arrays in memory (~190 MB for a 9 day
    forecast on the sa-west grid).

    Selection is positional (isel on a precomputed index) rather than
    label-based (sel) so we keep the Feb-29 fallback in
    build_doy_alignment_index for climatologies that don't carry day 366.

    Returns (climatology, threshold_90, threshold_10), each with shape
    (len(temp_time), s_rho, eta_rho, xi_rho), NaN on days with no climatology.
    """
    idx, valid = build_doy_alignment_index(temp_time, ds_clim['dayofyear'].values)

    out = []
    for var in ('climatology', 'threshold_90', 'threshold_10'):
        arr = ds_clim[var].isel(dayofyear=idx).values.astype('float32')
        if not valid.all():
            arr[~valid, ...] = np.nan
        out.append(arr)
    return tuple(out)


def process_level_batch(temp_slice, clim_seas_slice, clim_thresh_slice, is_cold, t_dates,
                         stats=True):
    """
    Detect MHW/MCS events for a (time, eta, xi) slab.

    stats : passed through to detect_events_with_climatology(). With the
            default of True the per-grid-cell statistics dicts are collected
            and returned as before. With stats=False they are not computed at
            all and the second return value is None - the 'categories' output
            is unaffected either way.
    """
    n_time, n_eta, n_xi = temp_slice.shape
    categories = np.zeros((n_time, n_eta, n_xi), dtype='int8')
    mhw_dicts  = [] if stats else None

    for i in range(n_eta):
        for j in range(n_xi):
            temp_ts        = temp_slice[:, i, j]
            clim_seas_ts   = clim_seas_slice[:, i, j]
            clim_thresh_ts = clim_thresh_slice[:, i, j]

            if np.all(np.isnan(temp_ts)):
                if stats:
                    mhw_dicts.append(None)
                continue

            mhw_ev, cats           = detect_events_with_climatology(
                temp_ts, clim_seas_ts, clim_thresh_ts, is_cold, t_dates, stats=stats
            )
            categories[:, i, j]   = cats
            if stats:
                mhw_dicts.append(mhw_ev)
    return categories, mhw_dicts

def process_single_level(level, n_levels, temp_level, clim_seas_level, clim_thresh_level,
                          is_cold, t_dates, batch_size, cat_slice, stats=False):
    """
    Detect MHW/MCS events for one vertical plane and store inside the in-memory array tracker.

    temp_level/clim_seas_level/clim_thresh_level are all (time, eta_rho, xi_rho)
    arrays for this level, already resampled to daily means and already aligned
    onto the forecast days by load_daily_baselines - so nothing is read from
    disk in here.

    stats defaults to False here (unlike the functions this calls, which
    default to True for backwards compatibility): this function only ever
    writes categories into cat_slice, so collecting the per-grid-cell event
    statistics would be pure overhead - and it is the dominant cost of the
    operational run.
    """
    n_eta = temp_level.shape[1]

    print(f"\n   Level {level}/{n_levels - 1}:")

    for i in range(0, n_eta, batch_size):
        end_i = min(i + batch_size, n_eta)

        categories, _ = process_level_batch(temp_level[:, i:end_i, :],
                                            clim_seas_level[:, i:end_i, :],
                                            clim_thresh_level[:, i:end_i, :],
                                            is_cold, t_dates, stats=stats)

        if is_cold:
            categories = -np.abs(categories).astype('int8')

        cat_slice[:, i:end_i, :] = categories
        print(f"      Rows {i:3d}-{end_i:3d} complete", end='\r')
    print(f"      Level {level} complete" + " " * 20)

def load_and_harmonize_baselines(clim_file, thresh_file):
    """
    Centralized memory-safe reader for baseline datasets using native dayofyear variables.
    """
    print(f'Opening climatology (dropping heavy variables): {clim_file}')
    vars_to_drop = ['u', 'v', 'salt', 'ubar', 'vbar']
    ds_clim_raw = xr.open_dataset(clim_file, drop_variables=vars_to_drop)
    
    print(f'Opening thresholds: {thresh_file}')
    ds_thresh_raw = xr.open_dataset(thresh_file)

    # Use native coordinate names directly
    ds_clim = xr.Dataset(coords=ds_clim_raw.coords)
    ds_clim['climatology'] = ds_clim_raw['temp'] if 'temp' in ds_clim_raw.data_vars else ds_clim_raw['climatology']
    
    if 'zeta' in ds_clim_raw.data_vars: 
        ds_clim['zeta'] = ds_clim_raw['zeta']
        
    ds_clim['threshold_90'] = ds_thresh_raw['threshold_90']
    ds_clim['threshold_10'] = ds_thresh_raw['threshold_10']
        
    for v in ['lon_rho', 'lat_rho', 'dayofyear']:
        if v in ds_clim_raw.coords and v not in ds_clim.coords:
            ds_clim = ds_clim.assign_coords({v: ds_clim_raw[v]})
            
    return ds_clim
