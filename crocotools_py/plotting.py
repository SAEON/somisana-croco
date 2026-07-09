import os
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import sys
import crocotools_py.postprocess as post

import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Wedge
from pathlib import Path
from matplotlib.collections import LineCollection

class LandmaskFeature(cfeature.GSHHSFeature):
    """from the OpenDrift code"""
    def __init__(self, scale='auto', globe=None, **kwargs):
        super().__init__(scale, **kwargs)

        if globe is not None:
            self._crs = ccrs.PlateCarree(globe=globe)

    def geometries(self):
        self.intersecting_geometries(extent=None)

    def intersecting_geometries(self, extent):
        global __polys__

        if self._scale == 'auto':
            scale = self._scale_from_extent(extent)
        else:
            scale = self._scale[0]
        return super().intersecting_geometries(extent)

def plot_land(ax, ocean_color = 'white', land_color = cfeature.COLORS['land'], lscale = 'auto', globe=None):
    """
    Plot the landmask or the shapes from GSHHG.
    """
    land = LandmaskFeature(scale=lscale, facecolor=land_color, globe=globe)
    ax.add_feature(land, zorder=2, facecolor=land_color, edgecolor='black')

def setup_plot(ax, fname, extents=None, add_land=True, land_color=('k', 0), lscale='h'):
    if extents is None:
        lon = post.get_grd_var(fname,'lon_rho')
        lat = post.get_grd_var(fname,'lat_rho')
        lon_min = min(np.ravel(lon))
        lon_max = max(np.ravel(lon))
        lat_min = min(np.ravel(lat))
        lat_max = max(np.ravel(lat))
        factor=0.05 
        dl = 0.5 * (lon_max - lon_min + lat_max - lat_min) * factor
        extents=[lon_min-dl,lon_max+dl,lat_min-dl,lat_max+dl]
    
    ax.set_extent(extents)
    if add_land:
        plot_land(ax,land_color=land_color,lscale=lscale)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='dimgrey', alpha=0.5, linestyle=':')
    gl.right_labels = False
    gl.top_labels = False
    
    return extents

def plot_var(ax,var_data,lon,lat, ticks = [], cmap = 'Spectral_r'):
    if len(ticks)==0:
        ticks = np.linspace(min(np.ravel(var_data)),max(np.ravel(var_data)),num=20)
    levs = np.array(ticks)
    cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)
    
    var_plt = ax.pcolormesh(lon, lat, var_data, cmap=cmap, norm=cmap_norm, transform=ccrs.PlateCarree())
    return var_plt

def plot_cbar(ax,var_plt, ticks=[], tick_font = 13, label='values', label_font=15, loc=None, aspect_ratio=1, orientation='vertical'):
    if loc is None:
        loc = [0.8, 0.2, 0.015, 0.6]
    cbarax = plt.gcf().add_axes(loc) 
    cbar_plt = plt.colorbar(var_plt, cbarax, ticks=ticks, orientation=orientation)
    cbar_plt.set_label(label, fontsize=label_font)
    cbar_plt.ax.tick_params(labelsize=tick_font)
    return cbar_plt

def plot_time(ax,time, loc=[0.5,1.01], time_fmt = '%Y-%m-%d %H:%M', time_font=15):
    time_plt = ax.text(loc[0], loc[1],  pd.Timestamp(time).strftime(time_fmt), ha='center', fontsize=time_font, transform=ax.transAxes)
    return time_plt

def plot_uv(ax,u,v,lon,lat, extents = None, skip_uv = 5, scale = 10, ref_vector = None, col = 'k'):
    width=0.0035  
    uv_plt = ax.quiver(lon[::skip_uv, ::skip_uv], lat[::skip_uv, ::skip_uv], u[::skip_uv, ::skip_uv], v[::skip_uv, ::skip_uv],
                      scale=scale, color=col, width=width, transform=ccrs.PlateCarree(), zorder=1)
    
    if ref_vector is not None:
        d_y=extents[3]-extents[2]
        d_x=extents[1]-extents[0]
        loc_x=extents[0]+d_x*0.04
        loc_y=extents[3]-d_y*0.04
        ax.quiver(np.array([loc_x]), np.array([loc_y]), np.array([ref_vector]), np.array([0]), 
                  scale=scale, color=col, width=width, transform=ccrs.PlateCarree(), zorder=1)
        loc_y_txt=loc_y-d_y*0.04
        ax.text(loc_x, loc_y_txt,  str(ref_vector)+' m s$^{-1}$', ha='left', transform=ccrs.PlateCarree())
    return uv_plt

def plot_isobaths(ax,fname,isobaths):
    h = post.get_grd_var(fname,'h')
    lon=post.get_grd_var(fname,'lon_rho')
    lat=post.get_grd_var(fname,'lat_rho')
    contour = ax.contour(lon,lat,h,isobaths,colors='k', linestyles='dashed', linewidths=1, transform=ccrs.PlateCarree())
    ax.clabel(contour, fmt='%1.0f', inline=True)

def get_uv_params(spd,scale_uv,ref_vector,skip_uv,num_vectors,aspect_ratio):
    if skip_uv is None:
        [Ny,Nx]=np.shape(spd)
        skip_uv=int(Ny/num_vectors*aspect_ratio)
    
    max_spd=np.nanmax(spd.ravel())
    if max_spd < 0.25:
        scale = 2.5; ref = 0.1
    elif max_spd < 0.5:
        scale = 5; ref = 0.25
    elif max_spd < 0.75:
        scale = 7.5; ref = 0.5
    elif max_spd < 1:
        scale = 10; ref = 0.5
    else:
        scale = 15; ref = 1
    
    if scale_uv is None: scale_uv = scale
    if ref_vector is None: ref_vector = ref
    return scale_uv, skip_uv, ref_vector

def plot(fname, ax=None, var='temp', grdname=None, time=slice(None), level=None, ticks = None, cmap = 'Spectral_r', extents = None, Yorig = None, add_cbar = True, cbar_loc = None, cbar_label = None, add_vectors = True, scale_uv = None, ref_vector = None, skip_uv = None, num_vectors=25, skip_time = 1, add_time_label = True, add_land = True, isobaths = None, jpg_out=None, gif_out=None, mp4_out=None):
    if grdname is None: grdname = fname
    if level is None:
        if isinstance(fname, xr.Dataset) or isinstance(fname, xr.DataArray): ds = fname.copy()
        else: ds = post.get_ds(fname,var)
        level = len(ds.s_rho) - 1
        ds.close()
        
    ds = post.get_var(fname,var,grdname=grdname,time=time,level=level,Yorig=Yorig)
    da_var=ds[var]
    time_var=np.atleast_1d(ds.time.values)
    lon = post.get_grd_var(grdname,'lon_rho').values
    lat = post.get_grd_var(grdname,'lat_rho').values
    
    data_plt = da_var.values if len(time_var)==1 else da_var.isel(time=0).values
    
    if ticks is None:
        is_anomaly = var.endswith('_anom')
        num_ticks = 10
        if is_anomaly:
            vmax = np.nanpercentile(abs(da_var), 99)
            vmax = round(vmax, 2 - int(np.floor(np.log10(abs(vmax)))) - 1)
            vmin = -vmax; cmap = 'bwr'; cbar_label = 'temperature anomaly ($\degree$C)'
        else:
            vmin = np.nanpercentile(da_var, 1); vmax = np.nanpercentile(da_var, 99)
            vmin = round(vmin, 2 - int(np.floor(np.log10(abs(vmin)))) - 1)
            vmax = round(vmax, 2 - int(np.floor(np.log10(abs(vmax)))) - 1)
    
        step = (vmax - vmin) / num_ticks
        step = round(step, 2 - int(np.floor(np.log10(abs(step)))) - 1)
        vmax = vmin + num_ticks * step
        ticks = np.arange(vmin, vmax + step / num_ticks, step)
    
    if extents is None:
        lon_min = min(np.ravel(lon)); lon_max = max(np.ravel(lon))
        lat_min = min(np.ravel(lat)); lat_max = max(np.ravel(lat))
        factor=0.05; dl = 0.5 * (lon_max - lon_min + lat_max - lat_min) * factor
        extents=[lon_min-dl,lon_max+dl,lat_min-dl,lat_max+dl]
    
    proj = ccrs.Mercator()
    x0, y0 = proj.transform_point(extents[0], extents[2], ccrs.PlateCarree())
    x1, y1 = proj.transform_point(extents[1], extents[3], ccrs.PlateCarree())
    aspect_ratio = abs(x1 - x0) / abs(y1 - y0)
    
    if ax is None:        
        fig_width = 6 * aspect_ratio if aspect_ratio > 1 else 6
        fig_height = 6 if aspect_ratio > 1 else 6 / aspect_ratio
        figsize = (0.1 * fig_width + fig_width + (0.2 * fig_width if add_cbar else 0.1 * fig_width), fig_height)
        fig = plt.figure(figsize=figsize) 
        ax = fig.add_axes([0.1, 0.1, 0.7 if add_cbar else 0.8, 0.8], projection=proj)
    
    extents = setup_plot(ax, grdname, extents = extents, add_land=add_land)
    var_plt = plot_var(ax,data_plt,lon,lat, ticks=ticks, cmap=cmap)
    if add_time_label: time_plt = plot_time(ax,time_var[0])
    if add_cbar: plot_cbar(ax,var_plt,label=cbar_label,ticks=ticks,aspect_ratio=aspect_ratio)
    
    if add_vectors:
        ds_uv = post.get_uv(fname,grdname=grdname,time=time,level=level,Yorig=Yorig)
        da_u=ds_uv.u; da_v=ds_uv.v
        u_plt = da_u.values if len(time_var)==1 else da_u.isel(time=0).values
        v_plt = da_v.values if len(time_var)==1 else da_v.isel(time=0).values
        spd_plt = np.sqrt(u_plt**2+v_plt**2)
        scale_uv, skip_uv, ref_vector = get_uv_params(spd_plt,scale_uv,ref_vector,skip_uv,num_vectors,aspect_ratio)
        uv_plt = plot_uv(ax,u_plt,v_plt,lon,lat, scale = scale_uv, skip_uv = skip_uv, ref_vector = ref_vector, extents=extents)
        
    if isobaths is not None: plot_isobaths(ax,grdname,isobaths)
    if len(time_var)==1:
        if jpg_out is not None: plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')
    else:
        def plot_tstep(i):
            var_i=da_var.isel(time=i).values
            if add_time_label: time_plt.set_text(pd.Timestamp(time_var[i]).strftime('%Y-%m-%d %H:%M'))
            var_plt.set_array(var_i.ravel())
            if add_vectors: uv_plt.set_UVC(da_u.isel(time=i).values[::skip_uv, ::skip_uv], da_v.isel(time=i).values[::skip_uv, ::skip_uv])
        anim = FuncAnimation(fig, plot_tstep, frames=range(0,len(time_var),skip_time)) 
        if gif_out is not None: anim.save(gif_out, writer='imagemagick')
        if mp4_out is not None: anim.save(mp4_out, writer="ffmpeg")

def plot_blk(croco_grd, croco_blk_file, var='wspd', figsize=(6,6), tstep=0, Yorig = 1993, ticks = [], cmap = 'Spectral_r', extents = None, cbar_loc = None, add_vectors = True, scale_uv = 150, skip_uv = 10, skip_time = 1, jpg_out=None, write_jpg=False, gif_out=None, write_gif=False, tstep_end=None):
    ds_croco_grd=xr.open_dataset(croco_grd)
    lon_rho=ds_croco_grd.lon_rho.values; lat_rho=ds_croco_grd.lat_rho.values
    angle=ds_croco_grd.angle.values; cos_a = np.cos(angle); sin_a = np.sin(angle)
    ds_croco_grd.close()

    ds_blk = xr.open_dataset(croco_blk_file, decode_times=False)
    ds_blk_days = np.float64(ds_blk.bulk_time.values)
    blk_time = datetime(Yorig,1,1) + timedelta(days = ds_blk_days[tstep])
    ds_blk_t = ds_blk.isel(bulk_time=tstep)
    var_data=ds_blk_t[var].values; u = post.u2rho(ds_blk_t.uwnd.values); v = post.v2rho(ds_blk_t.vwnd.values)
    u_rot = u * cos_a - v * sin_a; v_rot = v * cos_a + u * sin_a

    fig = plt.figure(); ax = plt.axes(projection=ccrs.Mercator()); setup_plot(ax, croco_grd)
    if len(ticks)==0: ticks = np.linspace(min(np.ravel(var_data)),max(np.ravel(var_data)),num=20)
    cmap_norm = mplc.BoundaryNorm(boundaries=np.array(ticks), ncolors=256)
    var_plt = ax.pcolormesh(lon_rho, lat_rho, var_data, cmap=cmap, norm=cmap_norm, transform=ccrs.PlateCarree())
    plot_cbar(var_plt,label=var,ticks=ticks,loc=cbar_loc)

    if add_vectors: uv_plt = ax.quiver(lon_rho[::skip_uv, ::skip_uv], lat_rho[::skip_uv, ::skip_uv], u_rot[::skip_uv, ::skip_uv], v_rot[::skip_uv, ::skip_uv], scale=scale_uv, transform=ccrs.PlateCarree(), zorder=1)
    time_plt = ax.text(0.5, 1.01, 'tstep = ' + str(tstep)+': '+datetime.strftime(blk_time, '%Y-%m-%d %H:%M:%S'), ha='center', transform=ax.transAxes)
    if write_jpg: plt.savefig(jpg_out or croco_blk_file.split('.nc')[0]+'_'+var+'_tstep'+str(tstep)+'.jpg',dpi=500,bbox_inches = 'tight')
    if write_gif:
        def plot_tstep(i):
            ds_blk_i = ds_blk.isel(bulk_time=i)
            u_i = post.u2rho(ds_blk_i.uwnd.values) * cos_a - post.v2rho(ds_blk_i.vwnd.values) * sin_a
            v_i = post.v2rho(ds_blk_i.vwnd.values) * cos_a + post.u2rho(ds_blk_i.uwnd.values) * sin_a
            time_plt.set_text('tstep = '+str(i)+': '+datetime.strftime(datetime(Yorig,1,1) + timedelta(days = ds_blk_days[i]), '%Y-%m-%d %H:%M:%S')) 
            var_plt.set_array(ds_blk_i[var].values.ravel())
            uv_plt.set_UVC(u_i[::skip_uv, ::skip_uv], v_i[::skip_uv, ::skip_uv])
        anim = FuncAnimation(fig, plot_tstep, frames=range(tstep, tstep_end or len(ds_blk.bulk_time), skip_time)) 
        anim.save(gif_out or croco_blk_file.split('.nc')[0]+'_'+var+'.gif', writer='imagemagick')

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

CMAP_9 = mplc.ListedColormap([FILL_C_EXT, FILL_C_SEV, FILL_C_STR, FILL_C_MOD,"#ffffff", FILL_MOD, FILL_STR, FILL_SEV, FILL_EXT])
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
        # Mask out land/fill values (-127) before looking for peaks
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

        # Concatenate entire timeline window cleanly
        all_dates = np.concatenate([obs_dates, fct_dates])
        all_temp  = np.concatenate([obs_temp, fct_temp])
        all_seas  = np.concatenate([obs_seas, fct_seas])
        all_h_thr = np.concatenate([obs_h_thr, fct_h_thr])
        all_c_thr = np.concatenate([obs_c_thr, fct_c_thr])

        fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
        ax.yaxis.grid(True, color="#cccccc", linewidth=0.7, zorder=0)
        ax.set_facecolor("white"); fig.patch.set_facecolor("white")

        # Plot climatology bounds
        ax.plot(all_dates, all_seas, ":", color="gray", label="Climatology", lw=1.5, zorder=2)
        ax.plot(all_dates, all_h_thr, "--", color="#d9534f", label="MHW threshold", lw=1.2, zorder=2)
        ax.plot(all_dates, all_c_thr, "--", color="#337ab7", label="MCS threshold", lw=1.2, zorder=2)

        # Draw historical vs forecast metrics
        if len(obs_dates) > 0:
            ax.plot(obs_dates, obs_temp, color="#777777", lw=2.5, label="SST observed", zorder=5)
        if len(fct_dates) > 0:
            ax.plot(fct_dates, fct_temp, color="black", lw=2.5, label="SST forecast", zorder=5)

        ax.axvline(today, color="black", lw=1.2, zorder=6)
        ax.text(today + pd.Timedelta(hours=4), ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1]-ax.get_ylim()[0]), "Today", va="top", ha="left", fontsize=9, fontweight="bold")

        ax.set_title(f"{site_name}  ({abs(data['lat']):.3f}°S, {data['lon']:.3f}°E)", fontsize=14, fontweight="bold", pad=10, color="#1a3a5c")
        ax.set_ylabel(f"Temperature [°C]", fontsize=11, fontweight="bold", color="#1a3a5c")
        ax.set_xlim(all_dates[0], all_dates[-1])
        
        # Format calendar days along X axis cleanly
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
    # Mask un-computed or masked grid locations safely (-127 back to NaN)
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

# --- Master Wrapper Function ---
def plot_operational_mhw_mcs(forecast_file, cat_file, clim_file, thresh_file, out_dir, start_date, end_date, Yorig=2000):
    print("Rendering Operational MHW/MCS Visuals")
    out_dir = Path(out_dir)
    
    # Standard NumPy-backed loading with explicit dropping to completely prevent memory overload
    vars_to_drop = ['u', 'v', 'salt', 'ubar', 'vbar']
    ds_clim_raw = xr.open_dataset(clim_file, drop_variables=vars_to_drop)
    ds_thresh_raw = xr.open_dataset(thresh_file)
    
    ds_clim = xr.Dataset(coords=ds_clim_raw.coords)
    
    if 'dayofyear' in ds_clim_raw.dims:
        ds_clim_raw = ds_clim_raw.rename_dims({'dayofyear': 'day_of_year'}).rename({'dayofyear': 'day_of_year'})

    ds_clim['climatology'] = ds_clim_raw['temp'] if 'temp' in ds_clim_raw.data_vars else ds_clim_raw['climatology']
    if 'zeta' in ds_clim_raw.data_vars: 
        ds_clim['zeta'] = ds_clim_raw['zeta']
        
    # Explicitly assign raw data pointers to enforce coordinate compatibility
    ds_clim['threshold_90'] = (('day_of_year', 's_rho', 'eta_rho', 'xi_rho'), ds_thresh_raw['threshold_90'].variable.data)
    ds_clim['threshold_10'] = (('day_of_year', 's_rho', 'eta_rho', 'xi_rho'), ds_thresh_raw['threshold_10'].variable.data)
        
    for v in ['lon_rho', 'lat_rho', 'day_of_year']:
        if v in ds_clim_raw.coords: 
            ds_clim = ds_clim.assign_coords({v: ds_clim_raw[v]})

    ds_cat  = xr.open_dataset(cat_file)
    ds_fcst = post.handle_time(post.get_ds(forecast_file, "temp"), Yorig=Yorig)
    
    lat = ds_fcst.lat_rho.values if "lat_rho" in ds_fcst else ds_fcst.lat.values
    lon = ds_fcst.lon_rho.values if "lon_rho" in ds_fcst else ds_fcst.lon.values
    if lat.ndim > 2: lat, lon = lat[0], lon[0]
    h = ds_fcst.h.values if "h" in ds_fcst else np.zeros_like(lat)
    if h.ndim > 2: h = h[0]
    nlev = len(ds_fcst.s_rho) if "s_rho" in ds_fcst else ds_fcst.dims.get("s_rho", 32)
    today = pd.Timestamp(ds_fcst.time.values[4]).normalize()

    # Stripped out '100m' varying level to streamline operational performance
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