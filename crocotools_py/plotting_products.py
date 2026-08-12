"""
Figures for the operational MHW/MCS products.

What is here is the plotting that crocotools_py.plotting cannot do generically:
per-site time-series against the day-of-year climatology, and the coastal flag
map. The spatial animations are plain map animations of a single variable, so
the workflow makes them by calling crocplot on the products file rather than
duplicating animation machinery here.

Everything reads the time-series file written by get_ts_multivar rather than
the model output, so the sites are whatever is listed in the domain's
timeseries_coords.csv and nothing here has any site coordinates baked into it.

The products themselves are built by crocotools_py.products.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, Wedge
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import crocotools_py.marineheatwaves as mhw

# Category colours, deepening with intensity. These match the 'mhw_mcs'
# colormap registered in crocotools_py.plotting, which the spatial animations
# use, so a category reads the same colour on a map as on a flag.
FILL_MOD   = "#ffc73e";  FILL_STR   = "#f77819"
FILL_SEV   = "#bf460c";  FILL_EXT   = "#4e1909"
FILL_C_MOD = "#a6d3e8";  FILL_C_STR = "#5da6c9"
FILL_C_SEV = "#2074a3";  FILL_C_EXT = "#103c68"

MHW_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_MOD,   2: FILL_STR,   3: FILL_SEV,   4: FILL_EXT}
MCS_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_C_MOD, 2: FILL_C_STR, 3: FILL_C_SEV, 4: FILL_C_EXT}

CAT_NAMES = ['None', 'Moderate', 'Strong', 'Severe', 'Extreme']


def open_ts(ts_file):
    """
    Open a time-series file written by get_ts_multivar, checking it is one.
    """
    ds = xr.open_dataset(ts_file) if not isinstance(ts_file, xr.Dataset) else ts_file
    if 'site' not in ds.dims:
        raise ValueError(
            f"{ts_file} has no 'site' dimension - these figures are made from the "
            'time-series file written by get_ts_multivar, not from the model output')
    return ds


def site_grid_indices(ds_ts):
    """
    The 0-based eta,xi indices of the grid cells the time-series were taken
    from, for indexing another field on the same grid (a climatology, say).

    CROCO's eta_rho/xi_rho coordinates are 1-based, so the index is the
    coordinate value minus one - see find_nearest_point().
    """
    return (ds_ts['eta_rho'].values.astype(int) - 1,
            ds_ts['xi_rho'].values.astype(int) - 1)


def site_flags(ds_ts, level=-1):
    """
    Summarise each site's worst event over the window.

    Returns {site_name: {'mode': 'MHW'|'MCS'|'None', 'max_cat': float}}, where
    mode is whichever of the two reached the higher category.
    """
    flags = {}
    for site in ds_ts['site'].values:
        cat = ds_ts['category'].sel(site=site).isel(s_rho=level).values.astype(float)
        # land/fill -> no event. xarray decodes the file's _FillValue (-127) to
        # NaN on read, so test for that rather than for the raw fill value
        cat[~np.isfinite(cat)] = 0
        mhw_days, mcs_days = cat[cat > 0], cat[cat < 0]
        max_mhw = float(np.max(mhw_days)) if len(mhw_days) else 0.0
        max_mcs = float(np.abs(np.min(mcs_days))) if len(mcs_days) else 0.0

        if max_mhw == 0.0 and max_mcs == 0.0:
            flags[str(site)] = {'mode': 'None', 'max_cat': 0}
        elif max_mhw >= max_mcs:
            flags[str(site)] = {'mode': 'MHW', 'max_cat': max_mhw}
        else:
            flags[str(site)] = {'mode': 'MCS', 'max_cat': max_mcs}
    return flags


def _flag_colour(mode, cat):
    if mode == 'None' or cat == 0:
        return MHW_FLAG_COLOURS[0]
    colours = MHW_FLAG_COLOURS if mode == 'MHW' else MCS_FLAG_COLOURS
    return colours[max(0, min(4, int(round(cat))))]


def _draw_gauge(ax_g):
    """The MHW/MCS category wheel used as the flag map legend."""
    cat_labels = {1: 'Mod', 2: 'Str', 3: 'Sev', 4: 'Ext'}
    ax_g.set_xlim(-1.55, 1.55); ax_g.set_ylim(-1.55, 1.55)
    ax_g.set_aspect('equal'); ax_g.axis('off')
    n, r_out, r_in = 4, 1.30, 0.48

    upper = np.degrees(np.linspace(0, np.pi, n + 1))
    for k, (th1, th2) in enumerate(zip(upper, upper[1:])):
        ax_g.add_patch(Wedge((0, 0), r_out, th1, th2, width=r_out - r_in,
                             fc=MHW_FLAG_COLOURS[k+1], ec='white', lw=0.8, zorder=1))
        mid = np.radians((th1 + th2) / 2); rl = (r_out + r_in) / 2
        ax_g.text(rl * np.cos(mid), rl * np.sin(mid), cat_labels[k+1], ha='center',
                  va='center', fontsize=7.0, fontweight='bold', color='white',
                  rotation=np.degrees(mid) - 90, zorder=3)

    lower = np.degrees(np.linspace(np.pi, 2 * np.pi, n + 1))
    for k, (th1, th2) in enumerate(zip(lower, lower[1:])):
        ax_g.add_patch(Wedge((0, 0), r_out, th1, th2, width=r_out - r_in,
                             fc=MCS_FLAG_COLOURS[n-k], ec='white', lw=0.8, zorder=1))
        mid = np.radians((th1 + th2) / 2); rl = (r_out + r_in) / 2
        ax_g.text(rl * np.cos(mid), rl * np.sin(mid), cat_labels[n-k], ha='center',
                  va='center', fontsize=7.0, fontweight='bold', color='white',
                  rotation=np.degrees(mid) - 90, zorder=3)

    ax_g.add_patch(plt.Circle((0, 0), r_in, fc=MHW_FLAG_COLOURS[0], ec='white', lw=1.0, zorder=2))
    ax_g.text(0, 0, 'None', ha='center', va='center', fontsize=8, fontweight='bold',
              color='white', zorder=4)
    ax_g.text(0,  r_out + 0.10, 'MHW', ha='center', va='bottom', fontsize=9,
              fontweight='bold', color=MHW_FLAG_COLOURS[2])
    ax_g.text(0, -(r_out + 0.10), 'MCS', ha='center', va='top', fontsize=9,
              fontweight='bold', color=MCS_FLAG_COLOURS[3])
    ax_g.set_title('Max Intensity\n(Discrete Flags)', fontsize=7, fontweight='bold',
                   pad=3, color='#1a3a5c')


def plot_coastal_flags(ts_file, out_path, level=-1, depth_name='Surface',
                       title=None, extent=None, margin=1.2):
    """
    Map of the worst MHW/MCS category reached at each site over the window.

    A coloured box is drawn along the coast for each stretch, coloured by the
    nearest site's peak category, with each site marked and labelled.

    Parameters
    ----------
    ts_file    : time-series file from get_ts_multivar, containing 'category'
    out_path   : png to write
    level      : s_rho index to summarise; -1 (default) is the surface, 0 the bottom
    depth_name : label for that level, used in the title and nothing else
    title      : first line of the title. Defaults to a generic description
    extent     : [lon0,lon1,lat0,lat1] map extent. Derived from the sites if None
    margin     : degrees of padding around the sites when deriving the extent
    """
    ds = open_ts(ts_file)
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)

    flags = site_flags(ds, level=level)
    names = [str(s) for s in ds['site'].values]
    lons = ds['lon_site'].values
    lats = ds['lat_site'].values
    dates = pd.to_datetime(ds['time'].values)

    # The sites are listed along the coast in the domain's coords file, so that
    # ordering is the path the flag boxes follow - no separate ordering needed.
    BOX_SIZE, OFFSHORE, BOX_STEP_DIST = 0.45, -0.10, 0.50
    dense_lons, dense_lats = [], []
    for k in range(len(names) - 1):
        t = np.linspace(0, 1, 100)
        dense_lons.extend(lons[k] + t * (lons[k+1] - lons[k]))
        dense_lats.extend(lats[k] + t * (lats[k+1] - lats[k]))
    dense_lons, dense_lats = np.array(dense_lons), np.array(dense_lats)

    dists = np.zeros(len(dense_lons))
    dists[1:] = np.cumsum(np.hypot(np.diff(dense_lons), np.diff(dense_lats)))

    boxes = []
    for bd in np.arange(0, dists[-1], BOX_STEP_DIST):
        cx = np.interp(bd, dists, dense_lons)
        cy = np.interp(bd, dists, dense_lats)
        nearest_site = names[int(np.argmin(np.hypot(cx - lons, cy - lats)))]
        info = flags.get(nearest_site, {'mode': 'None', 'max_cat': 0})
        idx = max(1, min(np.searchsorted(dists, bd), len(dists) - 1))
        slen = np.hypot(dense_lons[idx] - dense_lons[idx-1],
                        dense_lats[idx] - dense_lats[idx-1]) or 1.0
        # unit normal to the coast path, to push the boxes offshore
        px = -(dense_lats[idx] - dense_lats[idx-1]) / slen
        py = (dense_lons[idx] - dense_lons[idx-1]) / slen
        boxes.append((cx + px * OFFSHORE, cy + py * OFFSHORE,
                      _flag_colour(info['mode'], info['max_cat'])))

    if extent is None:
        extent = [lons.min() - margin, lons.max() + margin,
                  lats.min() - margin, lats.max() + margin]

    fig = plt.figure(figsize=(10, 13), dpi=150)
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    for cx, cy, col in boxes:
        ax.add_patch(FancyBboxPatch((cx - BOX_SIZE/2, cy - BOX_SIZE/2), BOX_SIZE, BOX_SIZE,
                                    boxstyle='round,pad=0.04', facecolor=col,
                                    edgecolor='white', linewidth=0.6, zorder=3,
                                    transform=ccrs.PlateCarree()))

    for name, site_lon, site_lat in zip(names, lons, lats):
        info = flags.get(name, {'mode': 'None', 'max_cat': 0})
        cat_int = max(0, min(4, int(round(info['max_cat']))))
        ax.plot(site_lon, site_lat, 'o', ms=4, color='white', zorder=8, mec='black',
                mew=0.8, transform=ccrs.PlateCarree())
        label = (f'{name}\nNone' if info['mode'] == 'None' or cat_int == 0
                 else f"{name}\n{info['mode']} – {CAT_NAMES[cat_int]}")
        ax.text(site_lon + 0.08, site_lat, label, ha='left', va='center', fontsize=6.5,
                fontweight='bold', color='#1a3a5c', zorder=9, transform=ccrs.PlateCarree(),
                path_effects=[pe.withStroke(linewidth=2, foreground='white')])

    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=4)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor='#555544', zorder=5)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='#aaaaaa', alpha=0.8,
                      linestyle='--', zorder=2)
    gl.top_labels = gl.right_labels = False
    _draw_gauge(fig.add_axes([0.53, 0.61, 0.28, 0.28]))

    if title is None:
        title = 'MHW / MCS Flag Map'
    ax.set_title(f"{title}  ·  {depth_name}\n"
                 f"{dates[0].strftime('%d %b')} – {dates[-1].strftime('%d %b %Y')}",
                 fontsize=12, color='#1a3a5c', pad=8)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  wrote {out_path}')


def _finish_timeseries_axes(fig, ax, dates, today, title, ylabel, ncol,
                            legend_order=None):
    """
    Shared styling for the per-site time-series figures.

    legend_order optionally fixes the order of the legend entries, so that the
    data lines lead and the reference lines follow, rather than taking whatever
    order they happened to be drawn in.
    """
    ax.axvline(today, color='black', lw=1.0, zorder=6)
    ax.text(today + pd.Timedelta(hours=4), ax.get_ylim()[1], 'Today',
            va='bottom', ha='left', fontsize=10)
    ax.set_title(title, fontsize=14, fontweight='bold', pad=10, color='#1a3a5c')
    ax.set_ylabel(ylabel, fontsize=11, fontweight='bold', color='#1a3a5c')
    ax.set_xlim(dates[0], dates[-1])
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d'))
    # the full dates run into each other on a narrow figure, so slant them
    fig.autofmt_xdate(rotation=30, ha='right')
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_color('#1a3a5c')
    ax.spines['bottom'].set_color('#1a3a5c')

    handles, labels = ax.get_legend_handles_labels()
    if legend_order is not None:
        ordered = [(h, l) for want in legend_order
                   for h, l in zip(handles, labels) if l == want]
        handles, labels = [h for h, _ in ordered], [l for _, l in ordered]
    ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.18),
              ncol=ncol, fontsize=9, frameon=False)


def plot_timeseries_mhw(ts_file, clim_file, thresh_file, out_dir, today,
                        level=-1, depth_name='Surface'):
    """
    Per-site temperature time-series against the day-of-year climatology and
    the MHW/MCS percentile thresholds, split into hindcast and forecast at
    'today'.

    Parameters
    ----------
    ts_file     : time-series file from get_ts_multivar, containing 'temp'
    clim_file   : day-of-year climatology file
    thresh_file : day-of-year percentile threshold file
    out_dir     : directory to write one png per site into
    today       : hindcast/forecast transition date
    level       : s_rho index; -1 (default) is the surface, 0 the bottom
    depth_name  : label for that level, used in the titles and file names
    """
    ds = open_ts(ts_file)
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    today = pd.Timestamp(today).normalize()

    ds_clim = mhw.load_and_harmonize_baselines(clim_file, thresh_file)
    # the climatology is on the model grid, so the sites index straight into it
    js, is_ = site_grid_indices(ds)
    dates = pd.to_datetime(ds['time'].values)
    # positional day-of-year lookup, so a climatology without day 366 still works
    doy_idx, _ = mhw.build_doy_alignment_index(ds['time'].values,
                                               ds_clim['dayofyear'].values)

    obs_m, fct_m = dates <= today, dates >= today

    for k, site in enumerate(ds['site'].values):
        name = str(site)
        temp = ds['temp'].sel(site=site).isel(s_rho=level).values
        at_site = dict(s_rho=level, eta_rho=js[k], xi_rho=is_[k])
        seas  = ds_clim['climatology'].isel(**at_site).values[doy_idx]
        h_thr = ds_clim['threshold_90'].isel(**at_site).values[doy_idx]
        c_thr = ds_clim['threshold_10'].isel(**at_site).values[doy_idx]

        fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
        ax.yaxis.grid(True, color='#cccccc', linewidth=0.7, zorder=0)
        ax.set_facecolor('white'); fig.patch.set_facecolor('white')

        # shade wherever the temperature crosses a threshold
        ax.fill_between(dates, temp, h_thr, where=(temp > h_thr), interpolate=True,
                        color='#f0ad4e', alpha=0.85, zorder=1)
        ax.fill_between(dates, temp, c_thr, where=(temp < c_thr), interpolate=True,
                        color='#5bc0de', alpha=0.85, zorder=1)

        ax.plot(dates, seas,  ':',  color='gray',    label='Climatology',   lw=1.5, zorder=2)
        ax.plot(dates, h_thr, '--', color='#d9534f', label='MHW threshold', lw=1.2, zorder=2)
        ax.plot(dates, c_thr, '--', color='#337ab7', label='MCS threshold', lw=1.2, zorder=2)
        if obs_m.any():
            ax.plot(dates[obs_m], temp[obs_m], color='#777777', lw=2.5,
                    label='SST observed', zorder=5)
        if fct_m.any():
            ax.plot(dates[fct_m], temp[fct_m], color='black', lw=2.5,
                    label='SST forecast', zorder=5)

        lat_used = float(ds['lat_rho'].sel(site=site))
        lon_used = float(ds['lon_rho'].sel(site=site))
        _finish_timeseries_axes(
            fig, ax, dates, today,
            f'{name}  ({abs(lat_used):.3f}°S, {lon_used:.3f}°E)',
            'Temperature [°C]', ncol=3,
            legend_order=['SST observed', 'MHW threshold', 'Climatology',
                          'SST forecast', 'MCS threshold'])

        plt.savefig(out_dir / f"{name.replace(' ', '_')}_{depth_name}_{today.strftime('%Y%m%d')}.png",
                    dpi=150, bbox_inches='tight')
        plt.close()

    ds_clim.close()
    print(f'  wrote {len(ds["site"])} {depth_name} time-series to {out_dir}')


def plot_timeseries_stratification(ts_file, out_dir, today, target_depth=None):
    """
    Per-site stratification time-series (bottom minus target-depth density),
    split into hindcast and forecast at 'today'.

    target_depth is only used to label the y axis. The depth the stratification
    was actually computed against is set by add_stratification, so pass the same
    value here if you want it named rather than described generically.
    """
    ds = open_ts(ts_file)
    if 'stratification' not in ds:
        print('  no stratification in the time-series file - skipping those figures')
        return
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    today = pd.Timestamp(today).normalize()

    dates = pd.to_datetime(ds['time'].values)
    obs_m, fct_m = dates <= today, dates >= today
    against = f'{target_depth:g}m' if target_depth is not None else 'target depth'
    ylabel = f'Stratification (bottom − {against} density) [kg m$^{{-3}}$]'

    for site in ds['site'].values:
        name = str(site)
        strat = ds['stratification'].sel(site=site).values

        fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
        ax.yaxis.grid(True, color='#cccccc', linewidth=0.7, zorder=0)
        ax.set_facecolor('white'); fig.patch.set_facecolor('white')
        ax.axhline(0, color='#999999', lw=1.0, ls=':', zorder=1)

        if obs_m.any():
            ax.plot(dates[obs_m], strat[obs_m], color='#777777', lw=2.5,
                    label='Stratification (observed)', zorder=5)
        if fct_m.any():
            ax.plot(dates[fct_m], strat[fct_m], color='black', lw=2.5,
                    label='Stratification (forecast)', zorder=5)

        lat_used = float(ds['lat_rho'].sel(site=site))
        lon_used = float(ds['lon_rho'].sel(site=site))
        _finish_timeseries_axes(
            fig, ax, dates, today,
            f'{name}  ({abs(lat_used):.3f}°S, {lon_used:.3f}°E)',
            ylabel, ncol=2)

        plt.savefig(out_dir / f"{name.replace(' ', '_')}_Stratification_{today.strftime('%Y%m%d')}.png",
                    dpi=150, bbox_inches='tight')
        plt.close()

    print(f'  wrote {len(ds["site"])} stratification time-series to {out_dir}')
