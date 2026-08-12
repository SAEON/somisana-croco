"""
Figures and animations for the operational MHW/MCS products.

What is here is the plotting that crocotools_py.plotting cannot do generically:
per-site time-series against the day-of-year climatology, and the coastal flag
maps. The spatial animations are plain map animations, so they are made by
calling crocplot on the products file from the workflow rather than by bespoke
code in here.

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

import crocotools_py.postprocess as post
import crocotools_py.marineheatwaves as mhw


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

TARGETS = {"Kleinsee":       (17.030382, -29.680623), "Hondeklipbaai":  (17.252461, -30.315292), "Doringbaai":     (18.213554, -31.814509),
    "Elandsbaai":     (18.30165,  -32.312317), "Laaiplek":       (18.125354, -32.742041), "Paternoster":    (17.870305, -32.777566),
    "Saldanha":       (17.929861, -33.074807), "Yzerfontein":    (18.13382,  -33.361876), "Bloubergstrand": (18.443896, -33.803906),
    "Oudekraal":      (18.342541, -33.980098), "Cape Point":     (18.46024,  -34.358313), "Simonstown":     (18.442294, -34.176514),
    "Strand":         (18.810174, -34.120553), "Hangklip":       (18.803882, -34.374716), "Kleinmond":      (19.026591, -34.355882),
    "Hermanus":       (19.256989, -34.425957), "Gansbaai":       (19.323381, -34.576985),}

WINDOW_DAYS = 10

FILL_MOD   = "#ffc73e";  FILL_STR   = "#f77819"

FILL_SEV   = "#bf460c";  FILL_EXT   = "#4e1909"

FILL_C_MOD = "#a6d3e8";  FILL_C_STR = "#5da6c9"

FILL_C_SEV = "#2074a3";  FILL_C_EXT = "#103c68"

MHW_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_MOD,   2: FILL_STR,   3: FILL_SEV,   4: FILL_EXT}

MCS_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_C_MOD, 2: FILL_C_STR, 3: FILL_C_SEV, 4: FILL_C_EXT}

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

def plot_timeseries_stratification(ds_cat, today, out_dir):
    """ds_cat is the same already-open products dataset -
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

def plot_operational_mhw_mcs(forecast_file, cat_file, clim_file, thresh_file, out_dir, today, Yorig=2000):
    """
    Render the operational MHW/MCS figures and animations from the output of
    the products file.

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

    # The spatial animations (categories, anomalies, thermal fronts and
    # stratification) are ordinary map animations of a single variable, so the
    # workflow makes them by calling crocplot on the products file rather than
    # duplicating the animation machinery here. See the products job in
    # .github/workflows/postprocess_croco.yml.

    # Stratification plots - only if stratification is present in the products file
    if 'stratification' in ds_cat:
        print("\nStratification Time Series")
        plot_timeseries_stratification(ds_cat, today, out_dir)

    ds_fcst.close(); ds_clim.close(); ds_cat.close()
    print(f"\nAll operational visuals saved cleanly to: {out_dir}")
