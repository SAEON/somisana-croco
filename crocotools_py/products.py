"""
Operational derived products from CROCO output.

This module assembles the daily-averaged 'products' file that accompanies each
operational forecast: MHW/MCS event categories, temperature/sea-level
anomalies, SST fronts and stratification.

The MHW/MCS detection algorithm itself lives in crocotools_py.marineheatwaves,
which is kept free of any operational plumbing so it can be imported on its own
for hindcast analysis. The figures made from the products live in
crocotools_py.plotting_products, so this module is computation only and does
not pull in matplotlib or cartopy.
"""

import os
import gc
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

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


def make_products_base(fname_in, fname_out, Yorig=2000, varList=None):
    """
    Create the products file: a daily average of the raw CROCO output, carrying
    through the full grid so that the result is itself a valid CROCO file.

    Parameters
    ----------
    fname_in  : raw CROCO output file(s) - anything postprocess.get_ds accepts
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

    print(f'Loading raw CROCO output: {fname_in}')
    ds_raw = post.handle_time(post.get_ds(fname_in, 'temp'), Yorig=Yorig)

    missing = [v for v in varList if v not in ds_raw]
    if missing:
        raise ValueError(f'variable(s) {missing} not found in {fname_in}')

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


def open_products(fname_out, fname_in=None, Yorig=2000):
    """
    Open the daily fields an add_* step needs, ready for it to append to
    fname_out.

    fname_in is where the data comes from, and can be any CROCO file(s) that
    postprocess.get_ds accepts - the raw model output, a set of files matched
    by a wildcard, or an existing products file. If it is left as None it
    defaults to fname_out, which is the usual operational case of a step adding
    to a products file that is already there.

    If fname_out does not exist it is created from fname_in first, so a step
    can be run standalone without make_products_base having been run. That is
    the 'either write it or add to it' behaviour which lets the add_* steps run
    in any order, be re-run individually, or be used one at a time on a
    hindcast.

    The contents are read into memory and the file is then closed, because the
    step is going to append to that same file and netcdf4 raises an HDF error
    if it is still open for reading. Closing alone is not enough to rely on:
    handle_time() returns datasets derived via assign_coords/sel, so the object
    the caller holds is not the one that owns the file handle. The file is
    daily and only a few hundred MB, so reading it up front is cheap.
    """
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    if fname_in is None:
        fname_in = fname_out
    if not Path(fname_out).exists():
        if fname_in == fname_out:
            raise ValueError(
                f'{fname_out} does not exist, and no input file was given to create '
                'it from. Either run make_products_base first, or pass fname_in so '
                'this step can create the products file itself.')
        print(f'{fname_out} does not exist yet - creating it from {fname_in}')
        make_products_base(fname_in, fname_out, Yorig=Yorig)
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


def add_stratification(fname_out, fname_in=None, Yorig=2000, target_depth=5.0,
                       temp_var='temp', salt_var='salt'):
    """
    Add daily bottom density, density at target_depth, and the difference
    between them (the stratification) to the products file.

    Density is in-situ EOS-80 density; the bottom is the bottom-most sigma
    layer (s_rho index 0 in CROCO), and the value at target_depth comes from
    the repo's own sigma-to-z interpolation (postprocess.hlev_xarray, the same
    routine get_var uses for level=-N). Cells where the seabed is shallower
    than target_depth come out as NaN rather than being extrapolated.

    Parameters
    ----------
    fname_out    : products file to add to (created first if it does not exist)
    fname_in     : CROCO file(s) to read the daily fields from; defaults to
                   fname_out, and is used to create it if it is not there yet
    Yorig        : origin year of the CROCO time axis
    target_depth : depth in metres (positive down) whose density is compared
                   against the bottom density
    """
    ds_prod = open_products(fname_out, fname_in, Yorig)

    if salt_var not in ds_prod:
        raise ValueError(
            f"'{salt_var}' is not in {fname_out}, so the stratification cannot be "
            'computed. Re-run make_products_base including salinity in --varList.')

    print(f'Computing stratification (bottom vs {target_depth}m density)')
    z = post.get_depths(ds_prod)  # (time, s_rho, eta_rho, xi_rho), negative down

    # Bottom: CROCO convention s_rho index 0 = bottom-most sigma layer
    density_bottom = eos80_density(ds_prod[temp_var].isel(s_rho=0),
                                   ds_prod[salt_var].isel(s_rho=0),
                                   -z.isel(s_rho=0))

    # hlev_xarray drops the depth dimension itself for a single depth
    target_z = -abs(target_depth)
    density_target = eos80_density(post.hlev_xarray(ds_prod[temp_var], z, target_z),
                                   post.hlev_xarray(ds_prod[salt_var], z, target_z),
                                   abs(target_depth))

    density_bottom = density_bottom.astype('float32')
    density_target = density_target.astype('float32')
    stratification = density_bottom - density_target

    mask_rho_2d = _mask_2d(ds_prod)
    for arr in (density_bottom, density_target, stratification):
        arr[:, mask_rho_2d == 0] = np.nan

    dims = ['time', 'eta_rho', 'xi_rho']
    ds_new = xr.Dataset({'density_bottom': (dims, density_bottom),
                         'density_target_depth': (dims, density_target),
                         'stratification': (dims, stratification)},
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


def _mask_2d(ds_prod):
    """The 2D land/sea mask of the products file (1 = sea, 0 = land)."""
    mask = ds_prod['mask_rho'].values
    return mask[0] if mask.ndim > 2 else mask


def add_anomalies(fname_out, clim_file, fname_in=None, Yorig=2000):
    """
    Add daily temperature and sea-level anomalies to the products file.

    Anomalies are taken against the day-of-year climatology - the same baseline
    the MHW/MCS categories are computed from - so that the anomaly and the
    category for a given day are consistent with each other.

    Parameters
    ----------
    fname_out : products file to add to (created first if it does not exist)
    clim_file : pre-built day-of-year climatology file
    fname_in  : CROCO file(s) to read the daily fields from; defaults to
                fname_out, and is used to create it if it is not there yet
    Yorig     : origin year of the CROCO time axis
    """
    ds_prod = open_products(fname_out, fname_in, Yorig)
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


def add_mhw_mcs(fname_out, clim_file, thresh_file, fname_in=None, Yorig=2000, batch_size=5):
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
    fname_in    : CROCO file(s) to read the daily fields from; defaults to
                  fname_out, and is used to create it if it is not there yet
    Yorig       : origin year of the CROCO time axis
    batch_size  : number of eta_rho rows detected at a time (memory vs speed)
    """
    ds_prod = open_products(fname_out, fname_in, Yorig)
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


def _masked_diff(field, axis):
    """
    Difference `field` along `axis` in index space, using only valid (non-NaN)
    neighbours: a centred difference where both neighbours are valid, a
    one-sided difference where only one is, and NaN where neither is.

    np.gradient cannot be used directly here because CROCO writes 0.0 over land
    rather than a fill value, and differencing sea against land would put an
    ~18 degC step across a single cell.
    """
    fwd = np.full(field.shape, np.nan)
    bwd = np.full(field.shape, np.nan)
    lo = [slice(None)] * field.ndim
    hi = [slice(None)] * field.ndim
    lo[axis] = slice(None, -1)
    hi[axis] = slice(1, None)
    diff = field[tuple(hi)] - field[tuple(lo)]
    fwd[tuple(lo)] = diff
    bwd[tuple(hi)] = diff
    both = np.isfinite(fwd) & np.isfinite(bwd)
    return np.where(both, (fwd + bwd) / 2,
                    np.where(np.isfinite(fwd), fwd, bwd))


def add_sst_front(fname_out, fname_in=None, Yorig=2000):
    """
    Add the daily surface thermal front magnitude to the products file.

    This is the horizontal gradient magnitude of the daily mean sea surface
    temperature, in degC/km, using the 2D grid metrics pm (1/dx, along xi) and
    pn (1/dy, along eta) to convert the per-cell differences to physical
    distances.

    Land is excluded before differencing, and the cells against it use a
    one-sided difference. Without that the coastal cells difference sea surface
    temperature against the 0.0 CROCO writes over land, which put a spurious
    front of ~6 degC/km along the whole coast - 200x the open-ocean value, on
    exactly the cells an upwelling product is about.

    Parameters
    ----------
    fname_out : products file to add to (created first if it does not exist)
    fname_in  : CROCO file(s) to read the daily fields from; defaults to
                fname_out, and is used to create it if it is not there yet
    Yorig     : origin year of the CROCO time axis
    """
    ds_prod = open_products(fname_out, fname_in, Yorig)

    T_daily = ds_prod.sizes['time']
    n_eta, n_xi = ds_prod.sizes['eta_rho'], ds_prod.sizes['xi_rho']
    mask_rho_2d = _mask_2d(ds_prod)

    print('Calculating daily surface thermal fronts (SST fronts)')
    sst = ds_prod['temp'].isel(s_rho=-1).values.astype('float64')
    if 'pm' not in ds_prod or 'pn' not in ds_prod:
        raise ValueError(
            f'{fname_out} has no pm/pn grid metrics, so the temperature gradient '
            'cannot be converted to degC/km. Re-run make_products_base, which '
            'carries them through from the raw CROCO output.')
    pm = ds_prod['pm'].values  # 1/dx, along xi  (eta_rho, xi_rho)
    pn = ds_prod['pn'].values  # 1/dy, along eta (eta_rho, xi_rho)
    if pm.ndim > 2:
        pm, pn = pm[0], pn[0]

    # blank out land before differencing - see _masked_diff()
    sst = np.where(mask_rho_2d[np.newaxis, :, :] == 1, sst, np.nan)

    out_sst_front = np.zeros((T_daily, n_eta, n_xi), dtype='float32')
    for t_idx in range(T_daily):
        d_eta = _masked_diff(sst[t_idx], axis=0)
        d_xi = _masked_diff(sst[t_idx], axis=1)
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




# Operational Plotting & Animation Routines




# sites of interest





    










