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
    Returns same structure as marineHeatWaves.detect
    
    Parameters:
    -----------
    temp_data : array (time,)
        Daily temperature time series
    clim_seas : array (time,)
        Seasonal climatology
    clim_thresh : array (time,)
        Threshold values (90th percentile for MHW, 10th for MCS)
    is_cold : bool
        True for cold spells, False for heat waves
    t_dates : array, optional
        Time vector in ordinal format
    
    Returns:
    --------
    mhw : dict
        Dictionary with all MHW properties (same as Hobday detect function)
    categories : array (time,)
        Event intensity categories (0=none, 1=Moderate, 2=Strong, 3=Severe, 4=Extreme)
    """
    
    n_time = len(temp_data)
    categories = np.zeros(n_time, dtype='int8')
    
    # Initialize MHW output structure
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
        'n_events': 0
    }
    
    # Handle NaN/missing data
    if np.all(np.isnan(temp_data)) or np.all(temp_data == 0):
        return mhw, categories
    
    # Replace NaN with climatology (as per original code)
    temp_clean = temp_data.copy()
    temp_clean[np.isnan(temp_clean)] = clim_seas[np.isnan(temp_clean)]
    
    # Flip for cold spells
    if is_cold:
        temp_clean = -1. * temp_clean
        clim_seas_use = -1. * clim_seas
        clim_thresh_use = -1. * clim_thresh
    else:
        clim_seas_use = clim_seas
        clim_thresh_use = clim_thresh
    
    # Find exceedances
    exceed_bool = temp_clean - clim_thresh_use
    exceed_bool[exceed_bool <= 0] = False
    exceed_bool[exceed_bool > 0] = True
    exceed_bool[np.isnan(exceed_bool)] = False
    
    # Find contiguous events
    events, n_events = ndimage.label(exceed_bool)
    
    # Process each event
    min_duration = 5
    categories_map = {'Moderate': 1, 'Strong': 2, 'Severe': 3, 'Extreme': 4}
    
    for ev in range(1, n_events + 1):
        event_indices = np.where(events == ev)[0]
        event_duration = len(event_indices)
        
        if event_duration < min_duration:
            continue
        
        tt_start = event_indices[0]
        tt_end = event_indices[-1]
        
        # Get event data
        temp_mhw = temp_clean[tt_start:tt_end+1]
        thresh_mhw = clim_thresh_use[tt_start:tt_end+1]
        seas_mhw = clim_seas_use[tt_start:tt_end+1]
        
        mhw_relSeas = temp_mhw - seas_mhw
        mhw_relThresh = temp_mhw - thresh_mhw
        mhw_relThreshNorm = (temp_mhw - thresh_mhw) / (thresh_mhw - seas_mhw)
        mhw_abs = temp_mhw
        
        # Find peak
        tt_peak = np.argmax(mhw_relSeas)
        
        # Store basic indices
        mhw['index_start'].append(int(tt_start))
        mhw['index_end'].append(int(tt_end))
        mhw['index_peak'].append(int(tt_start + tt_peak))
        
        # Time/date information
        if t_dates is not None:
            mhw['time_start'].append(int(t_dates[tt_start]))
            mhw['time_end'].append(int(t_dates[tt_end]))
            mhw['time_peak'].append(int(t_dates[tt_start + tt_peak]))
            mhw['date_start'].append(date.fromordinal(int(t_dates[tt_start])))
            mhw['date_end'].append(date.fromordinal(int(t_dates[tt_end])))
            mhw['date_peak'].append(date.fromordinal(int(t_dates[tt_start + tt_peak])))
        else:
            mhw['time_start'].append(tt_start)
            mhw['time_end'].append(tt_end)
            mhw['time_peak'].append(tt_start + tt_peak)
            mhw['date_start'].append(None)
            mhw['date_end'].append(None)
            mhw['date_peak'].append(None)
        
        # Duration
        mhw['duration'].append(event_duration)
        
        # Intensity metrics
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
        
        # Categories
        tt_peakCat = np.argmax(mhw_relThreshNorm)
        cats = np.floor(1. + mhw_relThreshNorm)
        cats = np.clip(cats, 1, 5)  # Clip to valid range (1-5, where 5+ = extreme)
        cat_names = ['Moderate', 'Strong', 'Severe', 'Extreme']
        # Ensure category is between 1 and 4
        peak_cat_value = np.clip(cats[tt_peakCat], 1, 4)
        cat_idx = int(peak_cat_value) - 1
        mhw['category'].append(cat_names[cat_idx])
        
        mhw['duration_moderate'].append(int(np.sum(cats == 1)))
        mhw['duration_strong'].append(int(np.sum(cats == 2)))
        mhw['duration_severe'].append(int(np.sum(cats == 3)))
        mhw['duration_extreme'].append(int(np.sum(cats >= 4)))
        
        # Rate of onset/decline
        if event_duration > 1:
            mhw['rate_onset'].append(float(mhw_relSeas[tt_peak] / (tt_peak + 1)))
            mhw['rate_decline'].append(float(mhw_relSeas[tt_peak] / (event_duration - tt_peak)))
        else:
            mhw['rate_onset'].append(0.)
            mhw['rate_decline'].append(0.)
        
        # Assign categories to full time series
        categories[tt_start:tt_end+1] = cats.astype('int8')
    
    mhw['n_events'] = len(mhw['duration'])
    
    return mhw, categories


def process_level_batch(temp_slice, clim_seas_slice, clim_thresh_slice, is_cold, t_dates):
    """
    Process a batch of grid cells (one level, subset of rows)
    
    Parameters:
    -----------
    temp_slice : array (time, eta, xi)
        Temperature data for this batch
    clim_seas_slice : array (time, eta, xi)
        Seasonal climatology aligned to temperature times
    clim_thresh_slice : array (time, eta, xi)
        Threshold values aligned to temperature times
    is_cold : bool
        True for MCS, False for MHW
    t_dates : array (time,)
        Time vector in ordinal format
    
    Returns:
    --------
    categories : array (time, eta, xi)
        Category values for all pixels in this batch
    mhw_dicts : list
        List of MHW dictionaries (one per pixel)
    """
    n_time, n_eta, n_xi = temp_slice.shape
    categories = np.zeros((n_time, n_eta, n_xi), dtype='int8')
    mhw_dicts = []
    
    for i in range(n_eta):
        for j in range(n_xi):
            temp_ts = temp_slice[:, i, j]
            clim_seas_ts = clim_seas_slice[:, i, j]
            clim_thresh_ts = clim_thresh_slice[:, i, j]
            
            # Skip if all NaN
            if np.all(np.isnan(temp_ts)):
                mhw_dicts.append(None)
                continue
            
            # Detect events
            mhw, cats = detect_events_with_climatology(
                temp_ts, clim_seas_ts, clim_thresh_ts, is_cold, t_dates
            )
            
            categories[:, i, j] = cats
            mhw_dicts.append(mhw)
    
    return categories, mhw_dicts


def align_climatology_to_temp(temp_time, doy_values, clim_data):
    """
    Align day-of-year climatology to temperature time series.
    
    Parameters:
    -----------
    temp_time : array
        Temperature time coordinates (datetime64)
    doy_values : array
        Day-of-year values from climatology (1-366)
    clim_data : array
        Climatology values corresponding to each day-of-year
    
    Returns:
    --------
    aligned_clim : array
        Climatology values aligned to temperature time series
    """
    n_time = len(temp_time)
    aligned_clim = np.zeros(n_time)
    
    # Create lookup dictionary for fast access
    clim_dict = {int(doy): clim_data[i] for i, doy in enumerate(doy_values)}
    
    for i, t in enumerate(temp_time):
        # Get day of year (1-366)
        doy = pd.Timestamp(t).dayofyear
        
        if doy in clim_dict:
            aligned_clim[i] = clim_dict[doy]
        else:
            # Handle leap year edge case (Feb 29 = day 60)
            # Use Feb 28 climatology for Feb 29 if not available
            if doy == 60 and 60 not in clim_dict:
                aligned_clim[i] = clim_dict.get(59, np.nan)
            else:
                print(f"Warning: DOY {doy} not found in climatology")
                aligned_clim[i] = np.nan
    
    return aligned_clim


def create_output_netcdf(output_file, n_time_daily, n_levels, n_eta, n_xi, 
                         ds_temp_daily, ds_temp, ds_clim, mode_name):
    """
    Create and initialize the output NetCDF file with proper structure.
    
    Parameters:
    -----------
    output_file : str or Path
        Path to output file
    n_time_daily : int
        Number of daily time steps
    n_levels : int
        Number of vertical levels
    n_eta : int
        Number of eta_rho grid points
    n_xi : int
        Number of xi_rho grid points
    ds_temp_daily : xarray.Dataset
        Daily-averaged temperature dataset
    ds_temp : xarray.Dataset
        Original temperature dataset
    ds_clim : xarray.Dataset
        Climatology dataset
    mode_name : str
        "MHW" or "MCS"
    
    Returns:
    --------
    nc_out : netCDF4.Dataset
        Opened NetCDF file ready for writing
    """
    nc_out = Dataset(output_file, mode='w', format='NETCDF4')
    
    # Create dimensions
    nc_out.createDimension('time', n_time_daily)
    nc_out.createDimension('s_rho', n_levels)
    nc_out.createDimension('eta_rho', n_eta)
    nc_out.createDimension('xi_rho', n_xi)
    
    # Create coordinate variables
    time_var = nc_out.createVariable('time', 'f8', ('time',))
    time_var.units = 'days since 1993-01-01'
    time_var.calendar = 'standard'
    time_values = (pd.to_datetime(ds_temp_daily.time.values) - 
                   pd.Timestamp('1993-01-01')).total_seconds() / 86400
    time_var[:] = time_values
    
    s_rho_var = nc_out.createVariable('s_rho', 'i4', ('s_rho',))
    s_rho_var[:] = np.arange(n_levels, dtype='int32')
    
    # Copy lat/lon from temperature file or climatology file
    if 'lat_rho' in ds_temp:
        lat_var = nc_out.createVariable('lat_rho', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
        lat_var[:] = ds_temp['lat_rho'].values
    elif 'lat_rho' in ds_clim:
        lat_var = nc_out.createVariable('lat_rho', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
        lat_var[:] = ds_clim['lat_rho'].values
        
    if 'lon_rho' in ds_temp:
        lon_var = nc_out.createVariable('lon_rho', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
        lon_var[:] = ds_temp['lon_rho'].values
    elif 'lon_rho' in ds_clim:
        lon_var = nc_out.createVariable('lon_rho', 'f4', ('eta_rho', 'xi_rho'), zlib=True)
        lon_var[:] = ds_clim['lon_rho'].values
    
    # Create category variable
    cat_var = nc_out.createVariable('category', 'i1', 
                                   ('time', 's_rho', 'eta_rho', 'xi_rho'),
                                   zlib=True, complevel=4, fill_value=0)
    cat_var.long_name = f'{mode_name} Category'
    cat_var.description = 'Event intensity: 0=None, 1=Moderate, 2=Strong, 3=Severe, 4=Extreme'
    
    # Add global attributes
    nc_out.title = f'{mode_name} Event Detection'
    nc_out.description = 'Detected using pre-computed climatology from 1993-2019'
    nc_out.created = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
    
    return nc_out


def process_single_level(level, n_levels, ds_temp_daily, ds_clim, temp_var_name,
                        doy_values, is_cold, t_dates, batch_size, nc_out):
    """
    Process a single vertical level and write results to NetCDF.
    
    Parameters:
    -----------
    level : int
        Vertical level index to process
    n_levels : int
        Total number of levels
    ds_temp_daily : xarray.Dataset
        Daily temperature data
    ds_clim : xarray.Dataset
        Climatology data
    temp_var_name : str
        Name of temperature variable
    doy_values : array
        Day-of-year values from climatology
    is_cold : bool
        True for MCS, False for MHW
    t_dates : array
        Time vector in ordinal format
    batch_size : int
        Number of rows to process at once
    nc_out : netCDF4.Dataset
        Output file to write results to
    
    Returns:
    --------
    None (writes directly to nc_out)
    """
    n_time_daily = ds_temp_daily.sizes['time']
    n_eta = ds_temp_daily.sizes['eta_rho']
    n_xi = ds_temp_daily.sizes['xi_rho']
    
    print(f"\n   Level {level}/{n_levels-1}:")
    
    # Get climatology for this level (day_of_year, eta_rho, xi_rho)
    clim_seas_level = ds_clim['climatology'].isel(s_rho=level).values
    clim_thresh_level = ds_clim['threshold'].isel(s_rho=level).values
    
    # Get the category variable from the file
    cat_var = nc_out.variables['category']
    
    # Process in spatial batches
    for i in range(0, n_eta, batch_size):
        end_i = min(i + batch_size, n_eta)
        
        # Get temperature slice (time, eta_rho, xi_rho)
        temp_slice = ds_temp_daily[temp_var_name].isel(
            s_rho=level, 
            eta_rho=slice(i, end_i)
        ).values
        
        # Align climatology to temperature time period
        clim_seas_slice = np.zeros((n_time_daily, end_i - i, n_xi))
        clim_thresh_slice = np.zeros((n_time_daily, end_i - i, n_xi))
        
        for lat_idx in range(end_i - i):
            for lon_idx in range(n_xi):
                # Extract climatology for this pixel (366 days)
                pixel_clim_seas = clim_seas_level[:, i + lat_idx, lon_idx]
                pixel_clim_thresh = clim_thresh_level[:, i + lat_idx, lon_idx]
                
                # Align to temperature dates
                clim_seas_slice[:, lat_idx, lon_idx] = align_climatology_to_temp(
                    ds_temp_daily.time.values, doy_values, pixel_clim_seas
                )
                clim_thresh_slice[:, lat_idx, lon_idx] = align_climatology_to_temp(
                    ds_temp_daily.time.values, doy_values, pixel_clim_thresh
                )
        
        # Detect events
        categories, mhw_dicts = process_level_batch(
            temp_slice, clim_seas_slice, clim_thresh_slice, is_cold, t_dates
        )
        
        # Write to file
        cat_var[:, level, i:end_i, :] = categories
        
        print(f"      Rows {i:3d}-{end_i:3d} complete", end='\r')
    
    print(f"      Level {level} complete" + " "*20)
    
    # Sync periodically
    if level % 5 == 0:
        nc_out.sync()
    
    gc.collect()


def prepare_temperature_data(ds_temp, temp_var_name, time_origin):
    """
    Prepare temperature data by fixing time coordinates and resampling to daily.
    
    Parameters:
    -----------
    ds_temp : xarray.Dataset
        Temperature dataset
    temp_var_name : str
        Name of temperature variable
    time_origin : str
        Origin time string for conversion
    
    Returns:
    --------
    ds_temp_daily : xarray.Dataset
        Daily-averaged temperature dataset
    ds_temp : xarray.Dataset
        Original dataset with fixed time coordinates
    """
    # Fix time coordinate if needed
    if not np.issubdtype(ds_temp.time.dtype, np.datetime64):
        try:
            seconds = ds_temp.time.values.astype(float)
            delta = seconds.astype("timedelta64[s]")
            origin = np.datetime64(time_origin)
            real_time = origin + delta
            ds_temp = ds_temp.assign_coords(time=("time", real_time))
        except Exception as e:
            print(f"   Warning: Could not convert time: {e}")
    
    # Resample to daily if needed
    ds_temp_daily = ds_temp.resample(time="1D").mean()
    ds_temp_daily = ds_temp_daily.interpolate_na(dim="time", limit=7)
    
    return ds_temp_daily, ds_temp