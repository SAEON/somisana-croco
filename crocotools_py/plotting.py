#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import os
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import sys
import crocotools_py.postprocess as post

def setup_plot(ax, fname, extents=[]):
    '''
    generic stuff applicable to all 2D plots
    '''
    # extents = [lon_min, lon_max, lat_min, lat_max]
    #
    # first need to get the domain extents if it's not set autmatically
    if len(extents) == 0:
        lon,lat,_ = post.get_lonlatmask(fname)
        lon_min = min(np.ravel(lon))
        lon_max = max(np.ravel(lon))
        lat_min = min(np.ravel(lat))
        lat_max = max(np.ravel(lat))
        factor=0.05 # factor of domain size used to get dl
        dl = 0.5 * (lon_max - lon_min + lat_max - lat_min) * factor
        extents=[lon_min-dl,lon_max+dl,lat_min-dl,lat_max+dl]
    
    ax.set_extent(extents)
    ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='dimgrey', alpha=0.5, linestyle=':')
    gl.right_labels = False
    gl.top_labels = False

def plot_var(ax,fname,var,
             tstep=0, # index or a datetime (in which case specify ref_date)
             level=0,
             ref_date = None, 
             ticks = [], # the ticks to plot
             cmap = 'Spectral_r',
             extents = []
             ):
    '''
    Add a variable to a 2D plot
    
    see post.get_var() for how inputs must be defined
    '''
    
    # get the data
    var_data = post.get_var(fname,var,tstep=tstep,level=level,ref_date=ref_date)
    lon = var_data.lon_rho.values
    lat = var_data.lat_rho.values
    var_data=var_data.values   
    
    # set up the plot
    setup_plot(ax,fname,extents = extents)
    
    # set up the cmap to handle non-uniform input ticks
    if len(ticks)==0:
        ticks = np.linspace(min(np.ravel(var_data)),max(np.ravel(var_data)),num=20)
    n_levels = len(ticks)
    #vmin = min(ticks) 
    #vmax = max(ticks)
    levs = np.array(ticks)
    cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)
    
    # plot the data
    var_plt = ax.pcolormesh(lon,
                              lat,
                              var_data,
                              cmap=cmap,
                              norm=cmap_norm,
                              transform=ccrs.PlateCarree())
    
    return var_plt

def plot_cbar(var_plt,
             ticks=[],
             tick_font = 12,
             label='values',
             label_font=14,
             loc=[1., 0.2, 0.02, 0.6], # [left, bottom, width, height]
             orientation='vertical'):
    
    '''
    Add a colorbar to a plot
    '''
    
    cbarax = plt.gcf().add_axes(loc) 
    cbar_plt = plt.colorbar(var_plt, cbarax,
                        ticks=ticks,
                        orientation=orientation)
    cbar_plt.set_label(label, fontsize=label_font)
    cbar_plt.ax.tick_params(labelsize=tick_font)
    
    return cbar_plt

def plot_time(ax,fname,
             loc=[0.5,1.01],
             tstep=0,
             ref_date = datetime(2000, 1, 1, 0, 0, 0),
             time_fmt = '%Y-%m-%d %H:%M',
             time_font=12):
    '''
    Add time text to a 2D plot
    '''
    
    time = post.get_time(fname, ref_date=ref_date,time_lims=tstep)[0]
    time_plt = ax.text(loc[0], loc[1],  datetime.strftime(time, time_fmt),
        ha='center', fontsize=time_font,
        transform=ax.transAxes)
    
    return time_plt

def plot_uv(ax,fname,
              tstep=0,
              level=0,
              ref_date=None,
              extents = [],
              scale = 10,
              skip = 1,
              col = 'k'
              ):
    '''
    Add vectors to a 2D plot
    '''
    
    # get the data
    u, v = post.get_uv(fname,tstep=tstep,level=level,ref_date=ref_date)
    lon = u.lon_rho.values
    lat = u.lat_rho.values
    u=u.values
    v=v.values
     
    # set up the plot
    setup_plot(ax, fname, extents = extents)
    
    # plot the data
    uv_plt = ax.quiver(lon[::skip, ::skip],
                      lat[::skip, ::skip],
                      u[::skip, ::skip],
                      v[::skip, ::skip],
                      scale=scale,
                      transform=ccrs.PlateCarree(), zorder=1)
    
    return uv_plt

def plot_isobaths(ax,fname,isobaths):
    h = post.get_var(fname,'h')
    lon=h.lon_rho.values
    lat=h.lat_rho.values
    ax.contour(lon,lat,h,isobaths,colors='k',transform=ccrs.PlateCarree())
    
def plot(fname,
        var='temp',
        figsize=(6,6), # (hz,vt)
        tstep=0, # the step to plot, or the first step to animate. Can be an integer or a datetime object
        level=None, # this gets set to surface level if not provided
        ticks = np.linspace(10,22,num=13), # the ticks to plot
        cmap = 'Spectral_r',
        extents = [],
        ref_date = None, # datetime, from CROCO model setup
        cbar_loc = [0.9, 0.2, 0.02, 0.6],
        cbar_label = 'temperature ($\degree$C)',
        add_vectors = True,
        scale_uv = 7,
        skip_uv = 3,
        skip_time = 1, # every nth time-step will be animated (if provided)
        isobaths = None, # optional list of isobaths to overlay over plot
        jpg_out=None,
        write_jpg=False,
        gif_out=None,
        write_gif=False,
        tstep_end=None, # The last timestep to animate. Only used if write_gif = True.
        ):
    '''
    this is a convenience function for doing a quick 2D plot with minimal coding.
    this might also be used as example code for doing your own plots 
    there's also an option to turn the plot into an animation
    '''
    
    fig = plt.figure(figsize=figsize) 
    ax = plt.axes(projection=ccrs.Mercator())
    
    # if no level provided then default to the surface level
    if level is None:
        print('no vertical level provided - defaulting to surface layer')
        with xr.open_dataset(fname) as ds: 
            level = len(ds.s_rho) - 1
    
    var_plt = plot_var(ax,fname,var,
             tstep=tstep,
             level=level,
             ref_date=ref_date, 
             ticks=ticks)
    
    time_plt = plot_time(ax,fname,tstep=tstep,ref_date=ref_date)
    
    plot_cbar(var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc)
    if add_vectors:
        uv_plt = plot_uv(ax,fname,
                      tstep=tstep,
                      level=level,
                      ref_date=ref_date,
                      scale = scale_uv,
                      skip = skip_uv,
                      )
        
    if isobaths is not None:
        plot_isobaths(ax,fname,isobaths)
    
    # write a jpg if specified
    if write_jpg:
        if jpg_out is None:
            jpg_out = fname.split('.nc')[0]+'_'+var+'_'+datetime.strftime(post.get_time(fname, ref_date)[tstep], '%Y%m%d_%H')+'.jpg'
        plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')
    
    # and/or write a gif if specified
    if write_gif: # do the animation
        # this only works if you stecify tstep as an integer, not a datetime
        def plot_tstep(i):
            # get the data for this time-step
            time_i = post.get_time(fname, ref_date)[i]
            var_i = post.get_var(fname,var,tstep=i,level=level,ref_date=ref_date)
            var_i = var_i.values
            u_i, v_i = post.get_uv(fname,tstep=i,level=level,ref_date=ref_date)
            u_i = u_i.values
            v_i = v_i.values
            
            # update the figure for this time-step
            time_plt.set_text(datetime.strftime(time_i, '%Y-%m-%d %H:%M:%S'))    
            var_plt.set_array(var_i.ravel())
            uv_plt.set_UVC(u_i[::skip_uv, ::skip_uv],
                                        v_i[::skip_uv, ::skip_uv])
        
        # animate
        if tstep_end is None: # if not defined then animate to end of file
            with xr.open_dataset(fname) as ds: 
                tstep_end = len(ds.time) - 1
        
        anim = FuncAnimation(
            fig, plot_tstep, frames=range(tstep,tstep_end,skip_time)) 
        
        if gif_out is None:
            gif_out = fname.split('.nc')[0]+'_'+var+'.gif'
        anim.save(gif_out, writer='imagemagick')
        
if __name__ == "__main__":
    
    fname='../output/avg_Y2011M2.nc'
    plot(fname,
          write_gif=True,
          level = -100,
          ticks = np.linspace(10,24,num=15),
          tstep_end=20,
          cbar_loc = [0.85, 0.2, 0.02, 0.6])
    
        
