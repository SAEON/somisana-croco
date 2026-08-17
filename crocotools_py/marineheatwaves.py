"""
Marine Heatwave (MHW) and Marine Cold Spell (MCS) detection and statistics.

Two separate things, in that order:

  detection  - per day, per grid cell, against a pre-built day-of-year
               climatology and percentile thresholds, following Hobday et al.
               (2016). This is what products.add_mhw_mcs calls to write the
               'category' variable into a products file.
  statistics - per year, per grid cell, computed from those categories once
               they exist. See annual_event_stats() and write_event_stats().

They are deliberately not one operation. Detection is only ever handed one
file's worth of days, and event statistics are only meaningful over a whole
record, so computing them together would count an event once per file it
touches. See the section header above annual_event_stats().

The operational forecast products built on top of the detection (the
daily-averaged products file, anomalies, SST fronts, stratification, and the
figures and animations) live in crocotools_py.products.
"""

import os
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import scipy.ndimage as ndimage

import crocotools_py.postprocess as post
import crocotools_py.define_attrs as cf_attrs

# Two separate exceedances no more than this many days apart are joined into a
# single event, as per Hobday et al. (2016). Kept as a constant so that the
# detection code and the metadata written by category_attrs() cannot drift.
MAX_GAP_DAYS = 2

# Hobday category names, indexed by category number - 1
CATEGORY_NAMES = ['Moderate', 'Strong', 'Severe', 'Extreme']

def category_attrs(min_duration=5, max_gap=MAX_GAP_DAYS):
    """
    The netcdf attributes describing the signed MHW/MCS 'category' variable.

    This is the single definition of that metadata, so that the operational
    products file and any hindcast product file written from this module
    describe the categories the same way. Pass the min_duration actually used
    for the run - it is recorded in the attributes, so a file can be read back
    without guessing which persistence criterion produced it.

    Returns a dict suitable for assigning straight onto DataArray.attrs.
    """
    attrs = {
        'long_name': 'MHW/MCS combined event category',
        'standard_name': 'status_flag',
        'units': '1',
        'flag_values': np.arange(-4, 5, dtype='int8'),
        'flag_meanings': ('marine_cold_spell_extreme marine_cold_spell_severe '
                          'marine_cold_spell_strong marine_cold_spell_moderate none '
                          'marine_heatwave_moderate marine_heatwave_strong '
                          'marine_heatwave_severe marine_heatwave_extreme'),
        'valid_range': np.array([-4, 4], dtype='int8'),
        'description': (
            'Sign: positive = marine heatwave (MHW, above the 90th percentile), '
            'negative = marine cold spell (MCS, below the 10th percentile), '
            '0 = neither.\n '
            'Magnitude: the Hobday et al. (2016) category, derived from the difference between the day-of-year '
            'climatological mean (clim), and the day-of-year 90th (MHW) or 10th (MCS) percentile thresholds (thresh) \n'
            ' dT = |threshold - clim|\n'
            '  1 Moderate: |T - clim| is 1x to 2x dT\n'
            '  2 Strong:   |T - clim| is 2x to 3x dT\n'
            '  3 Severe:   |T - clim| is 3x to 4x dT\n'
            '  4 Extreme:  |T - clim| is 4x dT or more \n'),
        'event_definition': (
            f'A run of days beyond the threshold counts as an event only if it lasts at '
            f'least {min_duration} day(s); exceedances separated by a gap of no more than '
            f'{max_gap} day(s) are joined into a single event, and the bridging days take '
            f'the minimum category of 1. Days outside an event are 0.'),
        'min_duration_days': np.int32(min_duration),
        'max_gap_days': np.int32(max_gap),
        'reference': ('Hobday et al. (2016), Progress in Oceanography 141, 227-238, '
                      'doi:10.1016/j.pocean.2015.12.014'),
    }
    if min_duration != 5:
        attrs['comment'] = (
            f'Here we are using a time persistence of {min_duration} day(s) to count '
            'as an event, rather than the 5 days of Hobday et al. (2016). This is useful '
            'in the context of operational output, where an event can be flagged in '
            'real time rather than waiting 5 days to flag retrospectively.')
    return attrs


# Core Heatwave Algorithms
def detect_events_with_climatology(temp_data, clim_seas, clim_thresh, is_cold,
                                   min_duration=5):
    """
    Detect MHW/MCS events using pre-computed climatology.

    Returns the per-timestep category time-series for one grid cell. The
    magnitude is the Hobday et al. (2016) category and the sign is applied by
    the caller, since a cell cannot be in both states at once.

    Categories follow Hobday et al. (2016): with dT = threshold - climatology,
    a day is category floor((T - climatology)/dT), i.e. 1 = Moderate (T is
    between 1x and 2x dT above the climatology), 2 = Strong, 3 = Severe and
    4 = Extreme. Category 4 is an open-ended top bucket (>= 4x dT) - there is
    no category 5 in the scheme.

    This function deliberately does NOT compute per-event statistics
    (durations, intensities, onset rates). Those are only meaningful over a
    whole record, and this is only ever handed one file's worth of days - so
    computing them here would silently count an event once per file it touches
    and split its duration across them. They are a separate step, run on the
    categories once they exist: see annual_event_stats().

    min_duration : number of days an exceedance must last before it counts as an
            event, default 5 as per Hobday et al. (2016).
    """
    n_time     = len(temp_data)
    categories = np.zeros(n_time, dtype='int8')

    if np.all(np.isnan(temp_data)) or np.all(temp_data == 0):
        return categories

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
            if 1 <= gap <= MAX_GAP_DAYS:
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

        # floor(1 + (T - thresh)/(thresh - clim)) == floor((T - clim)/dT), the
        # Hobday category. It is capped at 4 because Extreme is an open-ended
        # top bucket (>= 4x dT): capping at 5 instead let a category with no
        # name, no colour in the plotting colormap and no place in the file's
        # own metadata reach the output.
        mhw_relThreshNorm = (temp_mhw - thresh_mhw) / (thresh_mhw - seas_mhw)
        cats              = np.clip(np.floor(1.0 + mhw_relThreshNorm), 1, 4)
        categories[tt_start:tt_end + 1] = cats.astype('int8')

    return categories


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


def load_daily_baselines(ds_clim, temp_time,
                         variables=('climatology', 'threshold_90', 'threshold_10')):
    """
    Read the requested baseline fields ONCE, for only the days we actually
    need, across all vertical levels.

    Selection is positional (isel on a precomputed index) rather than
    label-based (sel) so we keep the Feb-29 fallback in
    build_doy_alignment_index for climatologies that don't carry day 366.

    Returns one array per name in `variables` (by default the climatology and
    both thresholds), each with shape (len(temp_time), s_rho, eta_rho, xi_rho)
    and NaN on days with no climatology. Computing anomalies only needs the
    climatology, so it asks for that alone rather than reading thresholds it
    will not use.
    """
    idx, valid = build_doy_alignment_index(temp_time, ds_clim['dayofyear'].values)

    out = []
    for var in variables:
        arr = ds_clim[var].isel(dayofyear=idx).values.astype('float32')
        if not valid.all():
            arr[~valid, ...] = np.nan
        out.append(arr)
    return tuple(out)


def process_level_batch(temp_slice, clim_seas_slice, clim_thresh_slice, is_cold,
                        min_duration=5):
    """
    Detect MHW/MCS events for a (time, eta, xi) slab, returning the categories.

    min_duration is passed through to detect_events_with_climatology().
    """
    n_time, n_eta, n_xi = temp_slice.shape
    categories = np.zeros((n_time, n_eta, n_xi), dtype='int8')

    for i in range(n_eta):
        for j in range(n_xi):
            temp_ts = temp_slice[:, i, j]

            if np.all(np.isnan(temp_ts)):
                continue

            categories[:, i, j] = detect_events_with_climatology(
                temp_ts, clim_seas_slice[:, i, j], clim_thresh_slice[:, i, j],
                is_cold, min_duration=min_duration)
    return categories

def process_single_level(level, n_levels, temp_level, clim_seas_level, clim_thresh_level,
                         is_cold, batch_size, cat_slice, min_duration=5):
    """
    Detect MHW/MCS events for one vertical plane and store inside the in-memory array tracker.

    temp_level/clim_seas_level/clim_thresh_level are all (time, eta_rho, xi_rho)
    arrays for this level, already resampled to daily means and already aligned
    onto the forecast days by load_daily_baselines - so nothing is read from
    disk in here.

    min_duration is passed through to detect_events_with_climatology() - see
    there for what it means and why the operational runs do not use the
    Hobday default.
    """
    n_eta = temp_level.shape[1]

    print(f"\n   Level {level}/{n_levels - 1}:")

    for i in range(0, n_eta, batch_size):
        end_i = min(i + batch_size, n_eta)

        categories = process_level_batch(temp_level[:, i:end_i, :],
                                         clim_seas_level[:, i:end_i, :],
                                         clim_thresh_level[:, i:end_i, :],
                                         is_cold, min_duration=min_duration)

        if is_cold:
            categories = -np.abs(categories).astype('int8')

        cat_slice[:, i:end_i, :] = categories
        print(f"      Rows {i:3d}-{end_i:3d} complete", end='\r')
    print(f"      Level {level} complete" + " " * 20)

def load_and_harmonize_baselines(clim_file, thresh_file=None):
    """
    Centralized memory-safe reader for baseline datasets using native dayofyear
    variables.

    thresh_file is optional: computing anomalies needs only the climatology,
    whereas MHW/MCS detection needs the percentile thresholds as well.
    """
    print(f'Opening climatology (dropping heavy variables): {clim_file}')
    vars_to_drop = ['u', 'v', 'salt', 'ubar', 'vbar']
    ds_clim_raw = xr.open_dataset(clim_file, drop_variables=vars_to_drop)

    # Use native coordinate names directly
    ds_clim = xr.Dataset(coords=ds_clim_raw.coords)
    ds_clim['climatology'] = ds_clim_raw['temp'] if 'temp' in ds_clim_raw.data_vars else ds_clim_raw['climatology']

    if 'zeta' in ds_clim_raw.data_vars:
        ds_clim['zeta'] = ds_clim_raw['zeta']

    if thresh_file is not None:
        print(f'Opening thresholds: {thresh_file}')
        ds_thresh_raw = xr.open_dataset(thresh_file)
        ds_clim['threshold_90'] = ds_thresh_raw['threshold_90']
        ds_clim['threshold_10'] = ds_thresh_raw['threshold_10']


    for v in ['lon_rho', 'lat_rho', 'dayofyear']:
        if v in ds_clim_raw.coords and v not in ds_clim.coords:
            ds_clim = ds_clim.assign_coords({v: ds_clim_raw[v]})
            
    return ds_clim


# Annual event statistics
#
# Once the categories exist, the event statistics are computed from them rather
# than by detecting all over again. A run of days with a non-zero category IS an
# event: min_duration and the gap bridging were already applied when the
# categories were written, so nothing here needs the temperature, the
# climatology or the thresholds - only 'category' and 'temp_anom', which is
# T - climatology and therefore the event intensity by definition.
#
# This is also why it must be a separate step rather than a flag on the
# detection. Detection is only ever handed one file's worth of days, so an event
# crossing a file boundary would be counted once per file it touches and its
# duration split between them. The statistics are computed on the concatenated
# record instead, where a run is the whole event - which is sound precisely
# because products.detection_window() makes the monthly categories identical to
# what a single pass over the whole record would give.
#
# Two attribution conventions, both recorded on the output variables:
#
#   day-based metrics   (days, days_<category>, intensity_*) count each day in
#                       the calendar year that day falls in. Unambiguous.
#   event-based metrics (count, duration_mean, duration_max) attribute the whole
#                       event to the year it STARTED in, so an event is counted
#                       exactly once no matter how many years it spans.
#
# A consequence worth knowing when reading the file: duration_mean is not
# days/count, because the two are binned differently at year boundaries.

# Statistics computed for each sign, as {suffix: (long_name, units, dtype)}.
#
# Note the units on the day counts and durations. 'days' (plural) is NOT usable
# here: xarray reads a variable whose units attribute is a plural time unit as a
# timedelta and decodes it to nanoseconds, so a count of 41724 days comes back
# as 3.6e18 to anyone who opens the file normally. A count is dimensionless, so
# it takes '1' with the meaning in the long_name; a duration really is in days,
# so it takes the UDUNITS canonical singular 'day', which is not decoded.
EVENT_STAT_VARS = {
    'days':                ('days in an event', '1', 'int16'),
    'count':               ('number of events starting', '1', 'int16'),
    'duration_mean':       ('mean event duration', 'day', 'float32'),
    'duration_max':        ('longest event duration', 'day', 'int16'),
    'intensity_mean':      ('mean intensity over event days', 'degC', 'float32'),
    'intensity_max':       ('peak intensity', 'degC', 'float32'),
    'intensity_cumulative': ('cumulative intensity over event days', 'degC day', 'float32'),
    'days_moderate':       ('days at category 1 (Moderate)', '1', 'int16'),
    'days_strong':         ('days at category 2 (Strong)', '1', 'int16'),
    'days_severe':         ('days at category 3 (Severe)', '1', 'int16'),
    'days_extreme':        ('days at category 4 (Extreme)', '1', 'int16'),
}


def _stats_one_sign(in_event, magnitude, intensity, year_of_day, n_years):
    """
    Annual event statistics for one sign (heatwaves or cold spells), for one
    vertical level.

    Parameters
    ----------
    in_event    : (time, eta, xi) bool - is this day part of an event
    magnitude   : (time, eta, xi) int8 - the Hobday category 1..4, 0 outside
    intensity   : (time, eta, xi) float - |T - climatology|, a positive
                  magnitude for both signs so that a cold spell 3 degC below the
                  climatology reads as 3, not -3
    year_of_day : (time,) int - index into the output year axis for each day
    n_years     : length of the output year axis

    Returns {suffix: (n_years, eta, xi) array} keyed as in EVENT_STAT_VARS.
    """
    n_time, n_eta, n_xi = in_event.shape
    shape = (n_years, n_eta, n_xi)

    out = {name: np.zeros(shape, dtype=('float64' if spec[2] == 'float32' else 'int64'))
           for name, spec in EVENT_STAT_VARS.items()}

    # --- day-based metrics: whole-year slices, since days are contiguous ------
    # a plain reduction per year, which is why these need no event logic at all
    bounds = np.searchsorted(year_of_day, np.arange(n_years + 1))
    for y in range(n_years):
        t0, t1 = bounds[y], bounds[y + 1]
        if t1 <= t0:
            continue
        ev = in_event[t0:t1]
        n_days = ev.sum(axis=0)
        out['days'][y] = n_days

        mag = magnitude[t0:t1]
        for cat, key in enumerate(('days_moderate', 'days_strong',
                                   'days_severe', 'days_extreme'), start=1):
            out[key][y] = (mag == cat).sum(axis=0)

        # Zero outside events rather than NaN, then reduce normally: within an
        # event the intensity is positive by construction for both signs, so a
        # plain sum and max are correct and there are no all-NaN slices to warn
        # about. Cells with no event days are set to NaN afterwards - 0 would
        # read as 'an event of zero intensity' rather than 'no event'.
        inten = np.where(ev, intensity[t0:t1], 0.0)
        none = n_days == 0
        total = inten.sum(axis=0)
        out['intensity_cumulative'][y] = np.where(none, np.nan, total)
        out['intensity_mean'][y] = np.where(none, np.nan,
                                            total / np.maximum(n_days, 1))
        out['intensity_max'][y] = np.where(none, np.nan, inten.max(axis=0))

    # --- event-based metrics: one pass, tracking the run length --------------
    # 'run' holds how many consecutive event days each cell is currently in. An
    # event ENDS on the last day it is still in one, and at that point run is its
    # full duration and t - run + 1 is the day it started - which is what gets
    # it attributed to its start year.
    run = np.zeros((n_eta, n_xi), dtype='int32')
    end_year, end_cell, end_dur = [], [], []
    for t in range(n_time):
        here = in_event[t]
        run = np.where(here, run + 1, 0)
        ends = here & (~in_event[t + 1] if t + 1 < n_time else here)
        if not ends.any():
            continue
        ii, jj = np.nonzero(ends)
        dur = run[ii, jj]
        end_year.append(year_of_day[t - dur + 1])
        end_cell.append(ii * n_xi + jj)
        end_dur.append(dur)

    if end_dur:
        ey = np.concatenate(end_year)
        ec = np.concatenate(end_cell)
        ed = np.concatenate(end_dur).astype('int64')
        flat = ey * (n_eta * n_xi) + ec
        size = n_years * n_eta * n_xi
        out['count'] = np.bincount(flat, minlength=size).reshape(shape)
        dur_sum = np.bincount(flat, weights=ed, minlength=size).reshape(shape)
        dur_max = np.zeros(size, dtype='int64')
        np.maximum.at(dur_max, flat, ed)
        out['duration_max'] = dur_max.reshape(shape)
        with np.errstate(invalid='ignore', divide='ignore'):
            out['duration_mean'] = np.where(out['count'] > 0,
                                            dur_sum / np.maximum(out['count'], 1),
                                            np.nan)
    else:
        out['duration_mean'][:] = np.nan

    return out


def annual_event_stats(category, temp_anom, times, mask=None):
    """
    Annual MHW and MCS event statistics for one vertical level.

    Parameters
    ----------
    category  : (time, eta, xi) signed Hobday category, as written by
                products.add_mhw_mcs. Land may be the int8 fill of -127 or NaN;
                either is handled, and `mask` blanks it either way
    temp_anom : (time, eta, xi) temperature anomaly T - climatology, as written
                by products.add_anomalies. This is the event intensity
    times     : (time,) datetime64 - the days the two arrays cover
    mask      : (eta, xi) land/sea mask, 1 = sea. Land is set to NaN (float
                outputs) or 0 (integer outputs) in the result

    Returns (stats, years) where stats is {variable_name: (n_years, eta, xi)}
    with names like 'mhw_days' and 'mcs_intensity_max', and years is the
    calendar years covered.
    """
    times = pd.DatetimeIndex(times)
    years = np.arange(times.year.min(), times.year.max() + 1)
    year_of_day = (times.year.values - years[0]).astype('int64')

    cat = np.asarray(category)
    # land may arrive as the int8 fill or as NaN depending on how it was read
    valid = np.isfinite(cat) if cat.dtype.kind == 'f' else (cat != -127)
    cat = np.where(valid, cat, 0).astype('int8')
    anom = np.asarray(temp_anom, dtype='float32')

    stats = {}
    for sign, prefix in ((1, 'mhw'), (-1, 'mcs')):
        in_event = (cat * sign) > 0
        magnitude = np.abs(cat) * in_event
        # intensity as a positive magnitude for both signs - a cold spell 3 degC
        # below the climatology reads as 3, matching how its category is |cat|
        intensity = anom * sign
        one = _stats_one_sign(in_event, magnitude, intensity, year_of_day, len(years))
        for suffix, arr in one.items():
            stats[f'{prefix}_{suffix}'] = arr

    if mask is not None:
        land = np.asarray(mask) == 0
        for name, arr in stats.items():
            dtype = EVENT_STAT_VARS[name.split('_', 1)[1]][2]
            arr[:, land] = np.nan if dtype == 'float32' else 0

    return stats, years


STATS_TITLE = 'SOMISANA annual marine heatwave / marine cold spell statistics'

# Heavy 4D fields in a products file that the statistics never look at. Dropped
# on open so that reading 33 years of monthly files does not pull them through.
STATS_DROP_VARS = ['temp', 'salt', 'u', 'v', 'zeta_anom', 'sst_front',
                   'stratification', 'density_deep', 'density_target_depth',
                   'temp_clim', 'temp_thresh_90', 'temp_thresh_10']


def mean_sigma_depths(ds):
    """
    Depth of each sigma level (metres, negative down) at the record-mean sea
    surface, as (s_rho, eta_rho, xi_rho).

    The statistics file has a year axis, not a time axis, so it cannot carry
    the time-varying depths a tier 1 file does. It does not need to: the sigma
    depths breathe with zeta, which on this grid is centimetres at the surface
    and millimetres at depth, against level spacings of metres to hundreds of
    metres. One static field computed at the mean sea surface describes the
    grid the statistics sit on, and is what makes the file independently
    usable without going back to the products files.
    """
    zeta_mean = ds['zeta'].mean(dim='time').values
    ds_one = ds.isel(time=[0]).copy()
    ds_one['zeta'] = (('time', 'eta_rho', 'xi_rho'), zeta_mean[np.newaxis, ...])
    return post.get_depths(ds_one).values[0]


def write_event_stats(fname_in, fname_out, Yorig=2000, doi_link=None,
                      compress=True):
    """
    Write annual MHW/MCS event statistics from a set of products files.

    The input is any number of products files carrying 'category' and
    'temp_anom' - typically the monthly files of a hindcast, given as a wildcard
    or a list. They are concatenated along time, so an event crossing a file
    boundary is one event, which is only sound because products.add_mhw_mcs was
    run with the neighbouring files as context (see products.detection_window).

    The output is NOT in CROCO format and cannot be fed to the regridding tiers:
    its vertical axis is still the model's sigma grid but its time axis is
    calendar years, so there is nothing for regrid_tier2 to interpolate onto
    depth levels. It is written to stand on its own instead, in the shape of a
    tier 1 file - lon_rho/lat_rho, h, mask and a 'depth' variable giving where
    each sigma level actually sits - so it can be read and plotted without
    reference to anything else.

    Parameters
    ----------
    fname_in  : products file(s) - a wildcard string or a list
    fname_out : statistics file to write
    Yorig     : origin year of the CROCO time axis of the input files
    doi_link  : bare DOI, written as a full https://doi.org/ URL
    compress  : deflate the output (default True; it is small either way)
    """
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    files = sorted(glob(fname_in)) if isinstance(fname_in, str) else sorted(fname_in)
    if not files:
        raise ValueError(f'no files matched {fname_in}')
    print(f'Reading {len(files)} products file(s)')

    # mask_and_scale=False keeps 'category' as int8: letting xarray decode its
    # -127 fill to NaN would promote 33 years of it to float and quadruple the
    # read for nothing, since the land mask says the same thing
    ds = xr.open_mfdataset(files, decode_times=False, mask_and_scale=False,
                           combine='by_coords', data_vars='minimal',
                           coords='minimal', compat='override',
                           drop_variables=STATS_DROP_VARS)
    ds = post.handle_time(ds, Yorig=Yorig)

    for needed in ('category', 'temp_anom'):
        if needed not in ds:
            raise ValueError(
                f"'{needed}' is not in the input files, so the statistics cannot "
                'be computed. They are written by products.add_mhw_mcs and '
                'products.add_anomalies respectively.')

    times = pd.DatetimeIndex(ds['time'].values)
    n_lev = ds.sizes['s_rho']
    n_eta, n_xi = ds.sizes['eta_rho'], ds.sizes['xi_rho']
    years = np.arange(times.year.min(), times.year.max() + 1)
    print(f'  {len(times)} days, {times[0].date()} to {times[-1].date()} '
          f'-> {len(years)} years, {n_lev} levels')

    gaps = np.diff(times.values).astype('timedelta64[D]').astype(int)
    if (gaps != 1).any():
        n_gap = int((gaps != 1).sum())
        print(f'  WARNING: the time axis is not continuous - {n_gap} gap(s). An '
              'event either side of a gap will be counted as two.')

    mask = ds['mask_rho'].values
    mask = mask[0] if mask.ndim > 2 else mask

    print('Computing the mean sigma level depths')
    depth = mean_sigma_depths(ds)

    out = {}
    for prefix in ('mhw', 'mcs'):
        for suffix, (_, _, dtype) in EVENT_STAT_VARS.items():
            out[f'{prefix}_{suffix}'] = np.zeros((len(years), n_lev, n_eta, n_xi),
                                                 dtype=dtype)

    for k in range(n_lev):
        print(f'  level {k + 1}/{n_lev}', end='\r')
        # one chunk per file per level, because products.default_encoding chunks
        # the 4D fields one sigma level at a time
        cat = ds['category'].isel(s_rho=k).values
        anom = ds['temp_anom'].isel(s_rho=k).values
        stats, _ = annual_event_stats(cat, anom, times, mask=mask)
        for name, arr in stats.items():
            out[name][:, k] = arr.astype(out[name].dtype)
        del cat, anom, stats
    print(f'  {n_lev} levels complete' + ' ' * 20)

    ds_out = xr.Dataset(coords={
        'year': ('year', years.astype('int32')),
        's_rho': ('s_rho', ds['s_rho'].values),
        'lon_rho': (('eta_rho', 'xi_rho'), ds['lon_rho'].values),
        'lat_rho': (('eta_rho', 'xi_rho'), ds['lat_rho'].values)})
    ds_out['year'].attrs = {'long_name': 'calendar year', 'units': '1'}

    ds_out['depth'] = (('s_rho', 'eta_rho', 'xi_rho'), depth.astype('float32'))
    ds_out['depth'].attrs = {
        'long_name': 'Depth of the sigma levels at the record-mean sea surface',
        'units': 'm', 'standard_name': 'depth', 'positive': 'up',
        'comment': 'Static, unlike the time-varying depth of a tier 1 file - see '
                   'marineheatwaves.mean_sigma_depths()'}
    ds_out['h'] = (('eta_rho', 'xi_rho'), ds['h'].values.astype('float32'))
    # the input is opened with mask_and_scale=False, which leaves _FillValue and
    # friends in attrs rather than moving them to encoding - carrying them over
    # would collide with the encoding set below
    ds_out['h'].attrs = {k: v for k, v in ds['h'].attrs.items()
                         if k not in ('_FillValue', 'scale_factor', 'add_offset',
                                      'missing_value')}
    ds_out['mask'] = (('eta_rho', 'xi_rho'), mask.astype('float32'))
    ds_out['mask'].attrs = {'long_name': 'mask on RHO-points', 'option_0': 'land',
                            'option_1': 'water', 'standard_name': 'land_binary_mask'}

    dims4 = ('year', 's_rho', 'eta_rho', 'xi_rho')
    for prefix, what in (('mhw', 'marine heatwave'), ('mcs', 'marine cold spell')):
        for suffix, (long_name, units, _) in EVENT_STAT_VARS.items():
            name = f'{prefix}_{suffix}'
            ds_out[name] = (dims4, out[name])
            by_event = suffix in ('count', 'duration_mean', 'duration_max')
            binned = ('the year the event started in, so an event spanning a year '
                      'boundary is counted once'
                      if by_event else 'the calendar year each day falls in')
            method = ('year: maximum' if suffix.endswith('_max') else
                      'year: mean' if suffix.endswith('_mean') else 'year: sum')
            attrs = {'long_name': f'Annual {what} {long_name}',
                     'units': units, 'cell_methods': method,
                     'comment': f'Binned by {binned}.'}
            if suffix.startswith('intensity'):
                attrs['comment'] += (' Intensity is |temperature - climatology|, '
                                     'a positive magnitude for both heatwaves and '
                                     'cold spells. NaN where there were no event '
                                     'days that year.')
            ds_out[name].attrs = attrs

    ds_out.attrs = cf_attrs.global_attrs(
        title=STATS_TITLE,
        source=f'{len(files)} products file(s): {files[0]} .. {files[-1]}',
        ds=ds_out, doi_link=doi_link,
        src_history=str(ds.attrs.get('history', '')),
        extra={'time_coverage_start': f'{times[0].date()}T00:00:00Z',
               'time_coverage_end': f'{times[-1].date()}T00:00:00Z',
               'time_coverage_resolution': 'P1Y'},
        action=f'write_event_stats: annual MHW/MCS event statistics over '
               f'{times[0].date()} to {times[-1].date()}, from {len(files)} '
               'products file(s)')

    encoding = {}
    for var in ds_out.data_vars:
        enc = {'zlib': True, 'complevel': 2} if compress else {}
        if ds_out[var].dtype.kind == 'f':
            enc['_FillValue'] = np.nan
        encoding[var] = enc

    out_file = Path(fname_out)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    if out_file.exists():
        out_file.unlink()
    print(f'Writing {fname_out}')
    ds_out.to_netcdf(fname_out, encoding=encoding)
    ds.close()
    print(f'Done: {fname_out} ({out_file.stat().st_size / (1024 ** 2):.1f} MB)')
