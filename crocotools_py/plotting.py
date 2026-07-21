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