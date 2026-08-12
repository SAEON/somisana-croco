"""
Operational derived products from CROCO output.

This module assembles the daily-averaged 'products' file that accompanies each
operational forecast: MHW/MCS event categories, temperature/sea-level
anomalies, SST fronts and stratification, plus the figures and animations made
from them.

The MHW/MCS detection algorithm itself lives in crocotools_py.marineheatwaves,
which is kept free of any operational plumbing so it can be imported on its own
for hindcast analysis.
"""

import os
import gc
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

# Plotting & GIS Libraries
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.patheffects as pe
from matplotlib.animation import FuncAnimation
from matplotlib.patches import FancyBboxPatch, Wedge
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Internal Dependencies
import crocotools_py.postprocess as post
import crocotools_py.marineheatwaves as mhw


# The products file
#
# All the operational products land in a single netcdf file (conventionally
# croco_avg_products.nc) which is itself a valid CROCO output file, just daily
# averaged rather than hourly. Keeping it in CROCO format means everything in
# postprocess.py (get_var, get_depths, regrid_tier*, ...) works on it directly,
# with no special cases for the products.
#
# make_products_base() creates that file from the raw CROCO output, and each
# add_* step then appends its own variables to it via append_to_products().

# Static grid variables carried through from the raw CROCO output. This is the
# set postprocess.get_ds() treats as essential, so that anything get_var() can
# do on croco_avg.nc it can also do on the products file. They are all small.
# 'theta_s', 'theta_b' and 'hc' are global attributes in CROCO output (and 'hc'
# is additionally a variable), so they are handled in both forms.
GRID_VARS = ['s_rho', 's_w', 'sc_r', 'sc_w', 'Cs_r', 'Cs_w',
             'hc', 'angle', 'h', 'f', 'pm', 'pn',
             'Vtransform', 'theta_s', 'theta_b',
             'lon_rho', 'lat_rho', 'mask_rho',
             'lon_u', 'lat_u', 'lon_v', 'lat_v']

# Variables daily-averaged into the products file by default. 'zeta' is added
# to whatever is asked for rather than listed here, because it is not optional:
# postprocess.get_depths() reads ds.zeta to work out the depths of the sigma
# levels, so without it nothing in the file can be interpolated to a z-level.
BASE_VARS = ['temp', 'salt']


def time_encoding(Yorig):
    """
    netcdf encoding for the time axis of the products file.

    Written as 'seconds since Yorig-01-01', as the rest of the repo does, so
    that reading the file back recovers the correct dates. The products file
    gets re-read by regrid_tier2 -> get_var -> get_ds(decode_times=False) +
    handle_time(), so it has to be readable by that path:
      - without any time encoding, xarray defaults to
        'days since <first timestamp>' and handle_time() reads those small
        integers as seconds, collapsing every date to 1970-01-01.
      - dtype must be float64 (as in native CROCO output), because
        handle_time() only applies Yorig when the raw values are float64
        (isinstance(np.float64, float) is True, but np.int32/np.float32 are
        not); an integer time axis silently falls through to being read as
        seconds since the 1970 epoch, shifting every date by Yorig-1970.
    """
    return {'units': f'seconds since {Yorig}-01-01 00:00:00',
            'calendar': 'standard',
            'dtype': 'f8'}


def default_encoding(ds):
    """
    Sensible netcdf encoding for the data variables of ds: deflate everything,
    and chunk 3D fields one sigma level at a time.

    The chunking matters. If a chunk spans the whole s_rho axis then reading a
    single level forces the whole file to be decompressed - that is exactly
    what made reading the day-of-year climatology files so slow.
    """
    encoding = {}
    for var in ds.data_vars:
        enc = {'zlib': True, 'complevel': 2}
        if ds[var].dtype.kind == 'f':
            enc['_FillValue'] = np.nan
        if 's_rho' in ds[var].dims:
            enc['chunksizes'] = tuple(1 if d == 's_rho' else ds.sizes[d]
                                      for d in ds[var].dims)
        encoding[var] = enc
    return encoding


def get_grid_vars(ds_raw):
    """
    Pull the static grid variables and vertical-coordinate attributes out of a
    raw CROCO dataset, as a dataset ready to be merged into the products file.
    """
    ds_grid = xr.Dataset()
    for var in GRID_VARS:
        if var in ds_raw.variables:
            # take dims/values/attrs only - carrying the DataArray across with
            # its own coords attached leaves xarray unable to work out whether
            # things like lon_u are coordinates or data variables on merge
            da = ds_raw[var]
            ds_grid[var] = (da.dims, da.values, dict(da.attrs))
    # keep whatever the raw file treats as a coordinate (lon_rho/lat_rho are
    # promoted to coordinates via the 'coordinates' attribute on the data
    # variables) as a coordinate here too, otherwise merging into a dataset
    # that already has them as coordinates is ambiguous
    ds_grid = ds_grid.set_coords([v for v in ds_grid.data_vars if v in ds_raw.coords])
    for attr in ['theta_s', 'theta_b', 'hc', 'Vtransform']:
        if attr in ds_raw.attrs:
            ds_grid.attrs[attr] = ds_raw.attrs[attr]
    return ds_grid


def make_products_base(fname, fname_out, Yorig=2000, varList=None):
    """
    Create the products file: a daily average of the raw CROCO output, carrying
    through the full grid so that the result is itself a valid CROCO file.

    Parameters
    ----------
    fname     : raw CROCO output file(s) - anything postprocess.get_ds accepts
    fname_out : products file to create (conventionally croco_avg_products.nc)
    Yorig     : origin year of the CROCO time axis
    varList   : variables to daily-average, default ['temp','salt']. 'zeta' is
                always included whether or not it is asked for, because
                get_depths() needs it.

    Any existing fname_out is overwritten - this is the step that creates the
    file, and the add_* steps then append to it.
    """
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    out_file = Path(fname_out)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    varList = list(BASE_VARS if varList is None else varList)
    if 'zeta' not in varList:
        varList.append('zeta')

    print(f'Loading raw CROCO output: {fname}')
    ds_raw = post.handle_time(post.get_ds(fname, 'temp'), Yorig=Yorig)

    missing = [v for v in varList if v not in ds_raw]
    if missing:
        raise ValueError(f'variable(s) {missing} not found in {fname}')

    # resample() already produces one bin per day over the full span, leaving
    # NaN for any day with no data - so no reindexing or gap filling is needed
    # (and a 'nearest' reindex would quietly duplicate a neighbouring day into
    # a real gap rather than showing it as missing).
    print(f'Computing daily means of: {", ".join(varList)}')
    ds_out = ds_raw[varList].resample(time='1D').mean().compute()

    print('Carrying through the grid variables')
    ds_grid = get_grid_vars(ds_raw)
    ds_out = ds_out.merge(ds_grid)
    ds_out.attrs.update(ds_grid.attrs)

    # resample() drops the variable attributes, so put them back
    for var in varList:
        ds_out[var].attrs = dict(ds_raw[var].attrs)

    encoding = default_encoding(ds_out)
    encoding['time'] = time_encoding(Yorig)

    print(f'Writing {fname_out}')
    if out_file.exists():
        out_file.unlink()
    ds_out.to_netcdf(fname_out, encoding=encoding)
    ds_raw.close()
    print(f'Done: {fname_out} ({out_file.stat().st_size / (1024 ** 2):.1f} MB, '
          f'{ds_out.sizes["time"]} days)')


def open_products(fname_out, fname_raw=None, Yorig=2000):
    """
    Open the products file, ready to be read by an add_* step.

    If the file does not exist yet it is created from the raw CROCO output
    first, so that a step can be run standalone without make_products_base
    having been run - this is the 'either write it or add to it' behaviour that
    lets the add_* steps run in any order or be re-run individually.

    The contents are read into memory and the file is then closed, because the
    step is going to append to that same file and netcdf4 raises an HDF error
    if it is still open for reading. Closing alone is not enough to rely on:
    handle_time() returns datasets derived via assign_coords/sel, so the object
    the caller holds is not the one that owns the file handle. The file is
    daily and only a few hundred MB, so reading it up front is cheap.
    """
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    if not Path(fname_out).exists():
        if fname_raw is None:
            raise ValueError(
                f'{fname_out} does not exist, and no raw CROCO file was given to '
                'create it from. Either run make_products_base first, or pass '
                'the raw file so this step can create the products file itself.')
        print(f'{fname_out} does not exist yet - creating it from {fname_raw}')
        make_products_base(fname_raw, fname_out, Yorig=Yorig)
    ds = post.handle_time(post.get_ds(fname_out, 'temp'), Yorig=Yorig).load()
    ds.close()
    return ds


def append_to_products(ds_new, fname_out, Yorig=2000, encoding=None):
    """
    Add the variables in ds_new to an existing products file (see
    open_products, which is what creates it if it is not there yet).

    Appending is done in place, so it costs the size of ds_new rather than the
    size of the file. Re-running a step overwrites its own variables.
    """
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    out_file = Path(fname_out)
    enc = dict(encoding) if encoding else default_encoding(ds_new)

    if not out_file.exists():
        raise ValueError(f'{fname_out} does not exist - open_products() should '
                         'have created it before anything is appended')

    # The file must not be held open while we append to it, or netcdf4 raises
    # an HDF error, so read what we need from it into memory and close it first.
    ds_existing = xr.open_dataset(fname_out, decode_times=False)
    existing_vars = set(ds_existing.data_vars)
    existing_time = np.asarray(post.handle_time(ds_existing, Yorig=Yorig).time.values)
    ds_existing.close()

    if 'time' in ds_new.dims:
        new_time = ds_new.time.values
        if len(new_time) != len(existing_time) or not np.array_equal(
                pd.to_datetime(new_time), pd.to_datetime(existing_time)):
            raise ValueError(
                f'time axis of the new variables does not match {fname_out}\n'
                f'  existing: {len(existing_time)} steps, '
                f'{pd.Timestamp(existing_time[0])} to {pd.Timestamp(existing_time[-1])}\n'
                f'  new     : {len(new_time)} steps, '
                f'{pd.Timestamp(new_time[0])} to {pd.Timestamp(new_time[-1])}\n'
                'The products file is built from a single forecast - re-run '
                'make_products_base if the forecast window has changed.')

    # coordinates already on disk must not be rewritten on append
    ds_new = ds_new.drop_vars([c for c in ds_new.coords if c in existing_vars
                               or c in ('time', 's_rho', 'eta_rho', 'xi_rho')],
                              errors='ignore')
    enc = {v: e for v, e in enc.items() if v in ds_new.data_vars}

    overwriting = sorted(set(ds_new.data_vars) & existing_vars)
    if overwriting:
        print(f'Overwriting existing variable(s): {", ".join(overwriting)}')

    ds_new.to_netcdf(fname_out, mode='a', encoding=enc)
    print(f'Appended to {fname_out}: {", ".join(ds_new.data_vars)}')


# Stratification (EOS-80 density)
#
# Density is computed using the full UNESCO EOS-80 equation of state
# (Millero & Poisson 1981 for density at 1 atm; Fofonoff & Millard 1983 for
# the pressure/compressibility correction) - no external dependency beyond
# numpy, and stays accurate at depth (important since "bottom" can be
# hundreds of meters deep here, not just a few meters).

def eos80_density(temp, salt, pressure):
    """
    Compute in-situ seawater density (kg/m^3) using the UNESCO EOS-80
    equation of state.

    Parameters
    ----------
    temp     : in-situ temperature (degC), array-like
    salt     : practical salinity (PSU), array-like
    pressure : pressure (dbar, positive downward). pressure (dbar) ~= depth (m)
               is used as a standard approximation for ocean-depth use.

    Returns
    -------
    rho : in-situ density (kg/m^3), same shape as inputs (broadcast-compatible)
    """
    T = np.asarray(temp, dtype='float64')
    S = np.asarray(salt, dtype='float64')
    P = np.asarray(pressure, dtype='float64') / 10.0  # dbar -> bar

    # Density of pure water at 1 atm (Bigg 1967 / Millero & Poisson 1981)
    rho_w = (999.842594 + 6.793952e-2 * T - 9.095290e-3 * T**2
             + 1.001685e-4 * T**3 - 1.120083e-6 * T**4 + 6.536332e-9 * T**5)

    # Density of seawater at 1 atm (Millero & Poisson 1981)
    B = (8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * T**2
         - 8.2467e-7 * T**3 + 5.3875e-9 * T**4)
    C = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2
    d0 = 4.8314e-4
    rho_st0 = rho_w + B * S + C * S**1.5 + d0 * S**2

    # Secant bulk modulus K(S,T,P) (Fofonoff & Millard 1983)
    Kw = (19652.21 + 148.4206 * T - 2.327105 * T**2
          + 1.360477e-2 * T**3 - 5.155288e-5 * T**4)
    F1 = 54.6746 - 0.603459 * T + 1.09987e-2 * T**2 - 6.1670e-5 * T**3
    F2 = 7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * T**2
    K_st0 = Kw + F1 * S + F2 * S**1.5

    Aw = (3.239908 + 1.43713e-3 * T + 1.16092e-4 * T**2 - 5.77905e-7 * T**3)
    G1 = 2.2838e-3 - 1.0981e-5 * T - 1.6078e-6 * T**2
    G2 = 1.91075e-4
    A = Aw + G1 * S + G2 * S**1.5

    Bw = 8.50935e-5 - 6.12293e-6 * T + 5.2787e-8 * T**2
    H1 = -9.9348e-7 + 2.0816e-8 * T + 9.1697e-10 * T**2
    B_coef = Bw + H1 * S

    K = K_st0 + A * P + B_coef * P**2

    rho = rho_st0 / (1.0 - P / K)
    return rho


def _compute_stratification_arrays(ds_prod, temp_var, salt_var, target_depth, batch_size):
    """
    Compute bottom density, density at target_depth, and the difference between
    them, from the daily fields in the products file. Returns three
    (time, eta_rho, xi_rho) arrays.

    The products file is already a valid CROCO file on the daily timeline, so
    get_depths() can be handed it directly - no need to assemble a stand-in
    dataset with the vertical grid copied across by hand.
    """
    print(f"\nComputing stratification (bottom vs {target_depth}m density)")
    T_daily, num_levels = ds_prod.sizes['time'], ds_prod.sizes['s_rho']
    n_eta, n_xi = ds_prod.sizes['eta_rho'], ds_prod.sizes['xi_rho']

    z_daily = post.get_depths(ds_prod)  # (time, s_rho, eta_rho, xi_rho), negative down

    out_density_bottom = np.full((T_daily, n_eta, n_xi), np.nan, dtype='float32')
    out_density_target = np.full((T_daily, n_eta, n_xi), np.nan, dtype='float32')

    for i in range(0, n_eta, batch_size):
        end_i = min(i + batch_size, n_eta)
        print(f"   Stratification rows {i:3d}-{end_i:3d}", end='\r', flush=True)

        temp_slice = ds_prod[temp_var].isel(eta_rho=slice(i, end_i)).values
        salt_slice = ds_prod[salt_var].isel(eta_rho=slice(i, end_i)).values
        z_slice = z_daily.isel(eta_rho=slice(i, end_i)).values

        # Bottom: CROCO convention s_rho index 0 = bottom-most sigma layer
        temp_bottom = temp_slice[:, 0, :, :]
        salt_bottom = salt_slice[:, 0, :, :]
        depth_bottom = -z_slice[:, 0, :, :]
        rho_bottom = eos80_density(temp_bottom, salt_bottom, depth_bottom)

        # Target depth: linear-interpolate between bracketing sigma levels
        # (vectorized: replaces per-column apply_along_axis/searchsorted
        # plus a triple nested Python loop with a boolean comparison +
        # np.take_along_axis gather)
        target_z = -abs(target_depth)
        above = (z_slice >= target_z)
        levs = np.argmax(above, axis=1)
        levs = np.clip(levs, 1, num_levels - 1)

        levs_up = levs[:, None, :, :]
        levs_down = levs_up - 1
        z_up      = np.take_along_axis(z_slice, levs_up, axis=1)[:, 0]
        z_down    = np.take_along_axis(z_slice, levs_down, axis=1)[:, 0]
        temp_up   = np.take_along_axis(temp_slice, levs_up, axis=1)[:, 0]
        temp_down = np.take_along_axis(temp_slice, levs_down, axis=1)[:, 0]
        salt_up   = np.take_along_axis(salt_slice, levs_up, axis=1)[:, 0]
        salt_down = np.take_along_axis(salt_slice, levs_down, axis=1)[:, 0]

        denom = z_up - z_down
        with np.errstate(invalid='ignore', divide='ignore'):
            frac = np.where(denom != 0, (target_z - z_down) / denom, np.nan)
        temp_at_target = (temp_down + frac * (temp_up - temp_down)).astype('float32')
        salt_at_target = (salt_down + frac * (salt_up - salt_down)).astype('float32')

        rho_target = eos80_density(temp_at_target, salt_at_target, abs(target_depth))

        out_density_bottom[:, i:end_i, :] = rho_bottom
        out_density_target[:, i:end_i, :] = rho_target

        del temp_slice, salt_slice, z_slice
        gc.collect()

    print("\n   Stratification complete" + " " * 20)
    out_stratification = out_density_bottom - out_density_target
    return out_density_bottom, out_density_target, out_stratification


def _mask_2d(ds_prod):
    """The 2D land/sea mask of the products file (1 = sea, 0 = land)."""
    mask = ds_prod['mask_rho'].values
    return mask[0] if mask.ndim > 2 else mask


def add_anomalies(fname_out, clim_file, fname=None, Yorig=2000):
    """
    Add daily temperature and sea-level anomalies to the products file.

    Anomalies are taken against the day-of-year climatology - the same baseline
    the MHW/MCS categories are computed from - so that the anomaly and the
    category for a given day are consistent with each other.

    Parameters
    ----------
    fname_out : products file to add to (created first if it does not exist)
    clim_file : pre-built day-of-year climatology file
    fname     : raw CROCO file, only needed if fname_out does not exist yet
    Yorig     : origin year of the CROCO time axis
    """
    ds_prod = open_products(fname_out, fname, Yorig)
    ds_clim = mhw.load_and_harmonize_baselines(clim_file)

    print('Computing daily temperature anomalies')
    clim_daily, = mhw.load_daily_baselines(ds_clim, ds_prod.time.values,
                                           variables=('climatology',))
    ds_new = xr.Dataset(
        {'temp_anom': (['time', 's_rho', 'eta_rho', 'xi_rho'],
                       (ds_prod['temp'].values - clim_daily).astype('float32'))},
        coords={'time': ds_prod.time.values})
    ds_new['temp_anom'].attrs = {'long_name': 'Sea Water Temperature Daily Anomaly',
                                 'units': 'degC'}

    if 'zeta' in ds_clim and 'zeta' in ds_prod:
        print('Computing daily sea level anomalies')
        # positional selection on the same alignment index the temperature
        # anomalies use, so both handle a climatology without day 366 the same
        idx, valid = mhw.build_doy_alignment_index(ds_prod.time.values,
                                                   ds_clim['dayofyear'].values)
        zeta_clim = ds_clim['zeta'].isel(dayofyear=idx).values.astype('float32')
        if not valid.all():
            zeta_clim[~valid, ...] = np.nan
        ds_new['zeta_anom'] = (['time', 'eta_rho', 'xi_rho'],
                               (ds_prod['zeta'].values - zeta_clim).astype('float32'))
        ds_new['zeta_anom'].attrs = {'long_name': 'Sea Surface Elevation Daily Anomaly',
                                     'units': 'm'}

    ds_prod.close()
    ds_clim.close()
    append_to_products(ds_new, fname_out, Yorig=Yorig)


def add_mhw_mcs(fname_out, clim_file, thresh_file, fname=None, Yorig=2000, batch_size=5):
    """
    Add signed Marine Heatwave / Marine Cold Spell event categories to the
    products file.

    'category' is +1..+4 for a heatwave, -1..-4 for a cold spell and 0 for
    neither, following Hobday et al. (2016). Where a cell is in neither state
    the MCS category is used, so the two never overlap. Land is set to -127.

    Parameters
    ----------
    fname_out   : products file to add to (created first if it does not exist)
    clim_file   : pre-built day-of-year climatology file
    thresh_file : pre-built day-of-year percentile threshold file
    fname       : raw CROCO file, only needed if fname_out does not exist yet
    Yorig       : origin year of the CROCO time axis
    batch_size  : number of eta_rho rows detected at a time (memory vs speed)
    """
    ds_prod = open_products(fname_out, fname, Yorig)
    ds_clim = mhw.load_and_harmonize_baselines(clim_file, thresh_file)

    T_daily, num_levels = ds_prod.sizes['time'], ds_prod.sizes['s_rho']
    n_eta, n_xi = ds_prod.sizes['eta_rho'], ds_prod.sizes['xi_rho']
    t_dates = np.array([d.toordinal() for d in pd.to_datetime(ds_prod.time.values)],
                       dtype=int)
    mask_rho_2d = _mask_2d(ds_prod)

    # Read the climatology/thresholds once, for the forecast days only, across
    # all levels - see marineheatwaves.load_daily_baselines() for why this
    # matters so much.
    print('Loading climatology and thresholds for the forecast days')
    clim_daily, thresh90_daily, thresh10_daily = mhw.load_daily_baselines(
        ds_clim, ds_prod.time.values)

    out_category = np.zeros((T_daily, num_levels, n_eta, n_xi), dtype='int8')

    print('Processing vertical planes')
    for k in range(num_levels - 1, -1, -1):
        temp_level = ds_prod['temp'].isel(s_rho=k).values
        clim_level = clim_daily[:, k, :, :]

        mhw_layer = np.zeros((T_daily, n_eta, n_xi), dtype='int8')
        mhw.process_single_level(k, num_levels, temp_level, clim_level,
                                 thresh90_daily[:, k, :, :], False, t_dates,
                                 batch_size, mhw_layer)

        mcs_layer = np.zeros((T_daily, n_eta, n_xi), dtype='int8')
        mhw.process_single_level(k, num_levels, temp_level, clim_level,
                                 thresh10_daily[:, k, :, :], True, t_dates,
                                 batch_size, mcs_layer)

        combined = mhw_layer.copy()
        combined[combined == 0] = mcs_layer[combined == 0]
        combined[:, mask_rho_2d == 0] = -127
        out_category[:, k, :, :] = combined
        gc.collect()

    ds_new = xr.Dataset(
        {'category': (['time', 's_rho', 'eta_rho', 'xi_rho'], out_category)},
        coords={'time': ds_prod.time.values})
    ds_new['category'].attrs = {
        'long_name': 'MHW_MCS Combined Event Categories',
        'description': 'Positive = Heatwave, Negative = Cold Spell, 0 = Neutral',
        'valid_range': np.array([-4, 4], dtype='int8'),
        'reference': ('Hobday et al. (2016), categories 1-4 = '
                      'Moderate/Strong/Severe/Extreme, where Extreme is >= 4x '
                      'the climatology-to-threshold difference')}

    ds_prod.close()
    ds_clim.close()
    append_to_products(ds_new, fname_out, Yorig=Yorig,
                       encoding={'category': {'zlib': True, 'complevel': 2,
                                              '_FillValue': -127,
                                              'chunksizes': (T_daily, 1, n_eta, n_xi)}})


def add_sst_front(fname_out, fname=None, Yorig=2000):
    """
    Add the daily surface thermal front magnitude to the products file.

    This is the horizontal gradient magnitude of the daily mean sea surface
    temperature, in degC/km, using the grid metrics pm/pn to convert the
    per-cell differences to physical distances.

    Parameters
    ----------
    fname_out : products file to add to (created first if it does not exist)
    fname     : raw CROCO file, only needed if fname_out does not exist yet
    Yorig     : origin year of the CROCO time axis
    """
    ds_prod = open_products(fname_out, fname, Yorig)

    T_daily = ds_prod.sizes['time']
    n_eta, n_xi = ds_prod.sizes['eta_rho'], ds_prod.sizes['xi_rho']
    mask_rho_2d = _mask_2d(ds_prod)

    print('Calculating daily surface thermal fronts (SST fronts)')
    sst = ds_prod['temp'].isel(s_rho=-1).values  # s_rho -1 = surface in CROCO
    pm = ds_prod['pm'].values if 'pm' in ds_prod else np.ones((n_eta, n_xi))
    pn = ds_prod['pn'].values if 'pn' in ds_prod else np.ones((n_eta, n_xi))
    if pm.ndim > 2:
        pm, pn = pm[0], pn[0]

    out_sst_front = np.zeros((T_daily, n_eta, n_xi), dtype='float32')
    for t_idx in range(T_daily):
        d_eta, d_xi = np.gradient(sst[t_idx])
        out_sst_front[t_idx] = np.hypot(d_xi * pm, d_eta * pn) * 1000.0
    out_sst_front = np.where(mask_rho_2d[np.newaxis, :, :] == 1, out_sst_front, np.nan)

    ds_new = xr.Dataset(
        {'sst_front': (['time', 'eta_rho', 'xi_rho'], out_sst_front.astype('float32'))},
        coords={'time': ds_prod.time.values})
    ds_new['sst_front'].attrs = {
        'long_name': 'Sea Surface Temperature Horizontal Front Magnitude',
        'units': 'degC / km'}

    ds_prod.close()
    append_to_products(ds_new, fname_out, Yorig=Yorig)


def add_stratification(fname_out, fname=None, Yorig=2000, target_depth=5.0,
                       batch_size=5, temp_var='temp', salt_var='salt'):
    """
    Add daily bottom density, density at target_depth, and the difference
    between them (the stratification) to the products file.

    Parameters
    ----------
    fname_out    : products file to add to (created first if it does not exist)
    fname        : raw CROCO file, only needed if fname_out does not exist yet
    Yorig        : origin year of the CROCO time axis
    target_depth : depth in metres (positive down) whose density is compared
                   against the bottom density
    batch_size   : number of eta_rho rows processed at a time
    """
    ds_prod = open_products(fname_out, fname, Yorig)

    if salt_var not in ds_prod:
        raise ValueError(
            f"'{salt_var}' is not in {fname_out}, so the stratification cannot be "
            'computed. Re-run make_products_base including salinity in --varList.')

    density_bottom, density_target, stratification = _compute_stratification_arrays(
        ds_prod, temp_var, salt_var, target_depth, batch_size)

    mask_rho_2d = _mask_2d(ds_prod)
    for arr in (density_bottom, density_target, stratification):
        arr[:, mask_rho_2d == 0] = np.nan

    dims = ['time', 'eta_rho', 'xi_rho']
    ds_new = xr.Dataset({'density_bottom': (dims, density_bottom),
                         'density_target_depth': (dims, density_target),
                         'stratification': (dims, stratification.astype('float32'))},
                        coords={'time': ds_prod.time.values})
    ds_new['density_bottom'].attrs = {'long_name': 'Bottom in-situ density',
                                      'units': 'kg m-3'}
    ds_new['density_target_depth'].attrs = {
        'long_name': f'In-situ density at {target_depth} m depth', 'units': 'kg m-3'}
    ds_new['stratification'].attrs = {
        'long_name': 'Water column stratification (bottom density minus density at target depth)',
        'units': 'kg m-3',
        'description': 'Positive = stably stratified, near-zero/negative = well-mixed',
        'target_depth_m': target_depth}

    ds_prod.close()
    append_to_products(ds_new, fname_out, Yorig=Yorig)


def detect_mhw_forecast(temp_file, clim_file, thresh_file, fname_out, temp_var='temp',
                        Yorig=2000, batch_size=5, salt_var='salt', target_depth=5.0,
                        compute_stratification=True):
    """
    Run the whole products pipeline in one call.

    Kept so that the existing operational workflow keeps working while it is
    moved onto the individual steps; new callers should use make_products_base
    followed by whichever add_* steps they want, which is what this does.

    Note that the output now also contains the daily mean temp/salt/zeta and
    the full CROCO grid, so it is a valid CROCO file rather than a bare
    collection of product fields.
    """
    make_products_base(temp_file, fname_out, Yorig=Yorig,
                       varList=[temp_var, salt_var] if compute_stratification else [temp_var])
    add_anomalies(fname_out, clim_file, Yorig=Yorig)
    add_mhw_mcs(fname_out, clim_file, thresh_file, Yorig=Yorig, batch_size=batch_size)
    add_sst_front(fname_out, Yorig=Yorig)
    if compute_stratification:
        add_stratification(fname_out, Yorig=Yorig, target_depth=target_depth,
                           batch_size=batch_size, temp_var=temp_var, salt_var=salt_var)
    print(f'Done: {fname_out} '
          f'({Path(fname_out).stat().st_size / (1024 ** 2):.1f} MB)')


# Operational Plotting & Animation Routines
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
        # land/fill -> no event. xarray decodes the file's _FillValue (-127) to
        # NaN on read, so test for that rather than for the raw fill value
        cat[~np.isfinite(cat)] = 0
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

        # Heatwave / coldspell shading — filled wherever SST crosses the threshold
        ax.fill_between(all_dates, all_temp, all_h_thr, where=(all_temp > all_h_thr),
                         interpolate=True, color="#f0ad4e", alpha=0.85, zorder=1)
        ax.fill_between(all_dates, all_temp, all_c_thr, where=(all_temp < all_c_thr),
                         interpolate=True, color="#5bc0de", alpha=0.85, zorder=1)

        ax.plot(all_dates, all_seas, ":", color="gray", label="Climatology", lw=1.5, zorder=2)
        ax.plot(all_dates, all_h_thr, "--", color="#d9534f", label="MHW threshold", lw=1.2, zorder=2)
        ax.plot(all_dates, all_c_thr, "--", color="#337ab7", label="MCS threshold", lw=1.2, zorder=2)

        if len(obs_dates) > 0:
            ax.plot(obs_dates, obs_temp, color="#777777", lw=2.5, label="SST observed", zorder=5)
        if len(fct_dates) > 0:
            ax.plot(fct_dates, fct_temp, color="black", lw=2.5, label="SST forecast", zorder=5)

        ax.axvline(today, color="black", lw=1.0, zorder=6)
        ax.text(today + pd.Timedelta(hours=4), ax.get_ylim()[1], "Today", va="bottom", ha="left", fontsize=10)

        ax.set_title(f"{site_name}  ({abs(data['lat']):.3f}°S, {data['lon']:.3f}°E)", fontsize=14, fontweight="bold", pad=10, color="#1a3a5c")
        ax.set_ylabel("Temperature [°C]", fontsize=11, fontweight="bold", color="#1a3a5c")
        ax.set_xlim(all_dates[0], all_dates[-1])

        ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d'))
        for spine in ("top", "right"): ax.spines[spine].set_visible(False)
        ax.spines['left'].set_color('#1a3a5c'); ax.spines['bottom'].set_color('#1a3a5c')

        # Explicit legend order to match: row1 = observed/MHW/climatology, row2 = forecast/MCS
        handles, labels = ax.get_legend_handles_labels()
        order = ["SST observed", "MHW threshold", "Climatology", "SST forecast", "MCS threshold"]
        ordered = [(h, l) for l in order for h, l2 in zip(handles, labels) if l2 == l]
        ax.legend([h for h, l in ordered], [l for h, l in ordered],
                   loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=9, frameon=False)

        plt.savefig(out_dir / f"{site_name.replace(' ', '_')}_{depth_name}_{today.strftime('%Y%m%d')}.png", dpi=150, bbox_inches="tight")
        plt.close()

# sites of interest
TARGETS = {"Kleinsee":       (17.030382, -29.680623), "Hondeklipbaai":  (17.252461, -30.315292), "Doringbaai":     (18.213554, -31.814509),
    "Elandsbaai":     (18.30165,  -32.312317), "Laaiplek":       (18.125354, -32.742041), "Paternoster":    (17.870305, -32.777566),
    "Saldanha":       (17.929861, -33.074807), "Yzerfontein":    (18.13382,  -33.361876), "Bloubergstrand": (18.443896, -33.803906),
    "Oudekraal":      (18.342541, -33.980098), "Cape Point":     (18.46024,  -34.358313), "Simonstown":     (18.442294, -34.176514),
    "Strand":         (18.810174, -34.120553), "Hangklip":       (18.803882, -34.374716), "Kleinmond":      (19.026591, -34.355882),
    "Hermanus":       (19.256989, -34.425957), "Gansbaai":       (19.323381, -34.576985),}

WINDOW_DAYS = 10

# Hobday category Colors
FILL_MOD   = "#ffc73e";  FILL_STR   = "#f77819"
FILL_SEV   = "#bf460c";  FILL_EXT   = "#4e1909"
FILL_C_MOD = "#a6d3e8";  FILL_C_STR = "#5da6c9"
FILL_C_SEV = "#2074a3";  FILL_C_EXT = "#103c68"

MHW_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_MOD,   2: FILL_STR,   3: FILL_SEV,   4: FILL_EXT}
MCS_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_C_MOD, 2: FILL_C_STR, 3: FILL_C_SEV, 4: FILL_C_EXT}

CMAP_9 = mplc.ListedColormap([FILL_C_EXT, FILL_C_SEV, FILL_C_STR, FILL_C_MOD,"#ffffff", FILL_MOD, FILL_STR, FILL_SEV, FILL_EXT])
CMAP_9.set_bad("white")
BNORM_9 = mplc.BoundaryNorm([-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5], CMAP_9.N)

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
        ax_g.text(0,  r_out + 0.10, "MHW", ha="center", va="bottom", fontsize=9, fontweight="bold", color=MHW_FLAG_COLOURS[2])
        ax_g.text(0, -(r_out + 0.10), "MCS", ha="center", va="top", fontsize=9, fontweight="bold", color=MCS_FLAG_COLOURS[3])
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

def animate_spatial_categories(cat_ds, ds_fcst, lat, lon, depth_name, lev, out_path):
    out_path = Path(out_path); times = pd.to_datetime(cat_ds.time.values)
    cat = cat_ds["category"].isel(s_rho=lev).values.astype(float)
    # land/fill is already NaN here - xarray decodes the file's _FillValue
    # (-127) on read - and the mask_rho step below covers it in any case
    mask = ds_fcst["mask_rho"].values if "mask_rho" in ds_fcst else np.ones_like(lat)
    if mask.ndim > 2: mask = mask[0]
    cat = np.where(mask[np.newaxis, :, :] == 1, cat, np.nan)

    fig = plt.figure(figsize=(10, 13)); ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(lon, lat, cat[0], transform=ccrs.PlateCarree(), cmap=CMAP_9, norm=BNORM_9, shading="auto")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6); ax.add_feature(cfeature.LAND, facecolor="white", edgecolor='black', zorder=2); ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=":")
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="#aaaaaa", alpha=0.8, linestyle="--", zorder=2); gl.top_labels = gl.right_labels = False

    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.04, ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4])
    cbar.ax.set_yticklabels(["Ext", "Sev", "Str", "Mod", "Neut", "Mod", "Str", "Sev", "Ext"])
    cbar.set_label("MCS (Cold)  ←  Intensity  →  MHW (Heat)")

    title = ax.set_title(f"MHW & MCS Categories ({depth_name})\nDate: {str(times[0])[:10]}")
    ani = FuncAnimation(fig, _update_spatial_frame, frames=len(times), fargs=(cat, times, mesh, title, depth_name), blit=False)
    ani.save(out_path, writer='ffmpeg', fps=1, dpi=120); plt.close(fig)
    
def _update_anomaly_frame(frame, anom_data, time_data, mesh_obj, title_obj):
    mesh_obj.set_array(anom_data[frame].ravel())
    title_obj.set_text(f"Sea Water Temperature Daily Anomaly (Surface)\nDate: {str(time_data[frame])[:10]}")
    return mesh_obj, title_obj

def animate_surface_anomalies(cat_ds, lat, lon, out_path):
    out_path = Path(out_path); times = pd.to_datetime(cat_ds.time.values)
    anom = cat_ds["temp_anom"].isel(s_rho=-1).values.astype(float)

    fig = plt.figure(figsize=(10, 13)); ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(lon, lat, anom[0], transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=-3.0, vmax=3.0, shading="auto")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6); ax.add_feature(cfeature.LAND, facecolor="lightgray", edgecolor='black', zorder=2)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="#aaaaaa", alpha=0.8, linestyle="--", zorder=2); gl.top_labels = gl.right_labels = False
    
    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.04); cbar.set_label("Temperature Anomaly (°C)")
    title = ax.set_title(f"Sea Water Temperature Daily Anomaly (Surface)\nDate: {str(times[0])[:10]}")
    ani = FuncAnimation(fig, _update_anomaly_frame, frames=len(times), fargs=(anom, times, mesh, title), blit=False)
    ani.save(out_path, writer='ffmpeg', fps=1, dpi=120); plt.close(fig)

def _update_frontal_frame(frame, front_data, time_data, mesh_obj, title_obj):
    mesh_obj.set_array(front_data[frame].ravel())
    title_obj.set_text(f"Horizontal Thermal Front Magnitude (Surface)\nDate: {str(time_data[frame])[:10]}")
    return mesh_obj, title_obj

def animate_surface_fronts(cat_ds, lat, lon, out_path):
    out_path = Path(out_path); times = pd.to_datetime(cat_ds.time.values)
    front = cat_ds["sst_front"].values.astype(float)

    fig = plt.figure(figsize=(10, 13)); ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(lon, lat, front[0], transform=ccrs.PlateCarree(), cmap='inferno', vmin=0.05, vmax=0.50, shading="auto")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6); ax.add_feature(cfeature.LAND, facecolor="lightgray", edgecolor='black', zorder=2)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="#aaaaaa", alpha=0.8, linestyle="--", zorder=2); gl.top_labels = gl.right_labels = False
    
    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.04); cbar.set_label("SST Gradient Magnitude (°C / km)")
    title = ax.set_title(f"Horizontal Thermal Front Magnitude (Surface)\nDate: {str(times[0])[:10]}")
    ani = FuncAnimation(fig, _update_frontal_frame, frames=len(times), fargs=(front, times, mesh, title), blit=False)
    ani.save(out_path, writer='ffmpeg', fps=1, dpi=120); plt.close(fig)


def plot_timeseries_stratification(ds_cat, today, out_dir):
    """ds_cat is the same already-open detect_mhw_forecast output dataset -
    reads 'stratification' directly from it, no separate file needed."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    today = pd.Timestamp(today)

    lon2d = ds_cat['lon_rho'].values
    lat2d = ds_cat['lat_rho'].values
    all_dates = pd.to_datetime(ds_cat['time'].values)

    for site_name, (site_lon, site_lat) in TARGETS.items():
        pj, pi = nearest(lon2d, lat2d, site_lon, site_lat)
        strat_ts = ds_cat['stratification'].isel(eta_rho=pj, xi_rho=pi).values

        obs_m = all_dates <= today
        fct_m = all_dates >= today

        fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
        ax.yaxis.grid(True, color='#cccccc', linewidth=0.7, zorder=0)
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')

        ax.axhline(0, color='#999999', lw=1.0, ls=':', zorder=1)
        if obs_m.any():
            ax.plot(all_dates[obs_m], strat_ts[obs_m], color='#777777', lw=2.5, label='Stratification (observed)', zorder=5)
        if fct_m.any():
            ax.plot(all_dates[fct_m], strat_ts[fct_m], color='black', lw=2.5, label='Stratification (forecast)', zorder=5)

        ax.axvline(today, color='black', lw=1.0, zorder=6)
        ax.text(today + pd.Timedelta(hours=4), ax.get_ylim()[1], 'Today', va='bottom', ha='left', fontsize=10)

        ax.set_title(f"{site_name}  ({abs(site_lat):.3f}°S, {site_lon:.3f}°E)", fontsize=14, fontweight='bold', pad=10, color='#1a3a5c')
        ax.set_ylabel('Stratification (bottom − 5m density) [kg m$^{-3}$]', fontsize=11, fontweight='bold', color='#1a3a5c')
        ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d'))
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9, frameon=False)

        plt.savefig(out_dir / f"{site_name.replace(' ', '_')}_Stratification_{today.strftime('%Y%m%d')}.png",
                    dpi=150, bbox_inches='tight')
        plt.close()


def _update_strat_frame(frame, data, time_data, mesh_obj, title_obj):
    mesh_obj.set_array(data[frame].ravel())
    title_obj.set_text(f"Water Column Stratification\nDate: {str(time_data[frame])[:10]}")
    return mesh_obj, title_obj


def animate_stratification(ds_cat, lat, lon, out_path):
    """ds_cat is the same already-open detect_mhw_forecast output dataset;
    lat/lon passed in already match plot_operational_mhw_mcs's own lat/lon
    arrays, so no separate lon_rho/lat_rho lookup needed here."""
    out_path = Path(out_path)
    times = pd.to_datetime(ds_cat['time'].values)
    data = ds_cat['stratification'].values

    fig = plt.figure(figsize=(10, 13))
    ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(lon, lat, data[0], transform=ccrs.PlateCarree(),
                          cmap='viridis', vmin=0, vmax=2.0, shading='auto')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.LAND, facecolor='white', edgecolor='black', zorder=2)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='#aaaaaa', alpha=0.8, linestyle='--', zorder=2)
    gl.top_labels = gl.right_labels = False

    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label('Stratification (bottom − 5m density) [kg m$^{-3}$]')
    title = ax.set_title(f"Water Column Stratification\nDate: {str(times[0])[:10]}")
    ani = FuncAnimation(fig, _update_strat_frame, frames=len(times),
                         fargs=(data, times, mesh, title), blit=False)
    ani.save(out_path, writer='ffmpeg', fps=1, dpi=120)
    plt.close(fig)


def plot_operational_mhw_mcs(forecast_file, cat_file, clim_file, thresh_file, out_dir, today, Yorig=2000):
    """
    Render the operational MHW/MCS figures and animations from the output of
    detect_mhw_forecast().

    'today' is the hindcast/forecast transition date used to split the
    time-series plots. The window shown in the flag map titles is taken from
    cat_file itself rather than being passed in, so the label cannot disagree
    with the data (it previously came from hardcoded +/-4 day offsets in the
    workflow, which did not match HDAYS/FDAYS).
    """
    print("Rendering Operational MHW/MCS Visuals")
    out_dir = Path(out_dir)

    ds_clim = mhw.load_and_harmonize_baselines(clim_file, thresh_file)
    ds_cat  = xr.open_dataset(cat_file)
    ds_fcst = post.handle_time(post.get_ds(forecast_file, "temp"), Yorig=Yorig)

    cat_dates = pd.to_datetime(ds_cat.time.values)
    start_date, end_date = cat_dates[0], cat_dates[-1]

    lat = ds_fcst.lat_rho.values if "lat_rho" in ds_fcst else ds_fcst.lat.values
    if lat.ndim > 2: lat, lon = lat[0], ds_fcst.lon_rho.values[0] if "lon_rho" in ds_fcst else ds_fcst.lon.values[0]
    else: lon = ds_fcst.lon_rho.values if "lon_rho" in ds_fcst else ds_fcst.lon.values
    
    h = ds_fcst.h.values if "h" in ds_fcst else np.zeros_like(lat)
    if h.ndim > 2: h = h[0]
    nlev = len(ds_fcst.s_rho) if "s_rho" in ds_fcst else ds_fcst.dims.get("s_rho", 32)
    today = pd.Timestamp(today).normalize()

    depth_levels = {"Surface": {"type": "fixed", "lev": nlev - 1},
        "Bottom":  {"type": "fixed", "lev": 0},}

    for depth_name, depth_info in depth_levels.items():
        print(f"\nProcessing Depth: {depth_name}...")
        lev_site = depth_info["lev"]
        # resample the whole level to daily means once per depth, not once per
        # site - every site reads a single point out of the same field
        ts = ds_fcst["temp"].isel(s_rho=lev_site).resample(time="1D").mean().load()
        all_dates = pd.to_datetime(ts.time.values)
        doy_all = doy_index(all_dates)
        obs_m, fct_m = all_dates <= today, all_dates >= today

        sites = {}
        for site_name, (site_lon, site_lat) in TARGETS.items():
            pj, pi = nearest(lon, lat, site_lon, site_lat)
            all_temps = ts.isel(eta_rho=pj, xi_rho=pi).values

            clim_profile   = ds_clim["climatology"].isel(s_rho=lev_site, eta_rho=pj, xi_rho=pi).values
            thresh90_profile = ds_clim["threshold_90"].isel(s_rho=lev_site, eta_rho=pj, xi_rho=pi).values
            thresh10_profile = ds_clim["threshold_10"].isel(s_rho=lev_site, eta_rho=pj, xi_rho=pi).values
            
            sites[site_name] = dict(
                pj=int(pj), pi=int(pi), lon=float(lon[pj, pi]), lat=float(lat[pj, pi]),
                obs_dates=pd.DatetimeIndex(all_dates[obs_m]), obs_temp=all_temps[obs_m],
                obs_seas=clim_profile[doy_all][obs_m],
                obs_h_thr=thresh90_profile[doy_all][obs_m],
                obs_c_thr=thresh10_profile[doy_all][obs_m],
                fct_dates=pd.DatetimeIndex(all_dates[fct_m]), fct_temp=all_temps[fct_m],
                fct_seas=clim_profile[doy_all][fct_m],
                fct_h_thr=thresh90_profile[doy_all][fct_m],
                fct_c_thr=thresh10_profile[doy_all][fct_m],)
            
        print("Time Series")
        plot_timeseries_multisite(sites, today, out_dir / depth_name, depth_name)

        print("Flag Maps")
        plot_flag_map(compute_site_flag_data(sites, ds_cat, depth_info["lev"]), today, start_date, end_date, out_dir / f"FlagMap_{depth_name}_{today.strftime('%Y%m%d')}.png", lat, lon, depth_name)

        print("Spatial Category Animation")
        animate_spatial_categories(ds_cat, ds_fcst, lat, lon, depth_name, depth_info["lev"], out_dir / f"Categories_Animation_{depth_name}.mp4")

        if depth_name == "Surface":
            print("Temperature Anomaly Animation")
            animate_surface_anomalies(ds_cat, lat, lon, out_dir / "Temperature_Anomaly_Animation_Surface.mp4")
            print("Thermal Front Animation")
            animate_surface_fronts(ds_cat, lat, lon, out_dir / "Thermal_Front_Animation_Surface.mp4")

    # Stratification plots - only if stratification is present in the categories/detection file)
    if 'stratification' in ds_cat:
        print("\nStratification Time Series")
        plot_timeseries_stratification(ds_cat, today, out_dir)

        print("Stratification Spatial Animation")
        animate_stratification(ds_cat, lat, lon, out_dir / "Stratification_Animation.mp4")

    ds_fcst.close(); ds_clim.close(); ds_cat.close()
    print(f"\nAll operational visuals saved cleanly to: {out_dir}")