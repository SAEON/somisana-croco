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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import sys
import crocotools_py.postprocess as post

def setup_plot(ax, fname, extents=None):
    '''
    generic stuff applicable to all 2D plots
    '''
    # extents = [lon_min, lon_max, lat_min, lat_max]
    #
    # first need to get the domain extents if it's not set autmatically
    if extents is None:
        lon = post.get_grd_var(fname,'lon_rho')
        lat = post.get_grd_var(fname,'lat_rho')
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
             grdname=None,
             tstep=0, # index or a datetime (in which case specify ref_date)
             level=0,
             ref_date = None, 
             ticks = [], # the ticks to plot
             cmap = 'Spectral_r',
             extents = None
             ):
    '''
    Add a variable to a 2D plot
    
    see post.get_var() for how inputs must be defined
    '''
    
    # get the data
    if var=='spd':
        u, v = post.get_uv(fname,grdname=grdname,tstep=tstep,level=level,ref_date=ref_date)
        var_data = np.sqrt(u**2+v**2)
    else:
        var_data = post.get_var(fname,var,grdname=grdname,tstep=tstep,level=level,ref_date=ref_date)
    lon = post.get_grd_var(grdname,'lon_rho').values
    lat = post.get_grd_var(grdname,'lat_rho').values
    var_data=var_data.values   
    
    # set up the plot
    setup_plot(ax,grdname,extents = extents)
    
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

def plot_cbar(ax,var_plt,
             ticks=[],
             tick_font = 12,
             label='values',
             label_font=14,
             loc=None, # [left, bottom, width, height]
             aspect_ratio=1,
             orientation='vertical'):
    
    '''
    Add a colorbar to a plot
    '''
    if loc is None:
        # dynamically position to colorbar
        colorbar_offset = 0.05  # Fixed offset from the right of the plot
        x_position = 0.8 + colorbar_offset * (1 / aspect_ratio)# adjust `x` based on aspect ratio, to no avail
        x_thickness = 0.02/aspect_ratio
        loc = [x_position, 0.2, x_thickness, 0.6]
    
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
              grdname=None,
              tstep=0,
              level=0,
              ref_date=None,
              extents = None,
              skip_uv = 5,
              scale = 10,
              ref_vector = None,
              col = 'k'
              ):
    '''
    Add vectors to a 2D plot
    '''
    
    # get the data
    u, v = post.get_uv(fname,grdname=grdname,tstep=tstep,level=level,ref_date=ref_date)
    lon = post.get_grd_var(grdname,'lon_rho').values
    lat = post.get_grd_var(grdname,'lat_rho').values
    u=u.values
    v=v.values
     
    # set up the plot
    setup_plot(ax, grdname, extents = extents)
        
    # plot the data
    width=0.0035  # Explicitly set the vector width for consistency between data and reference vector
    uv_plt = ax.quiver(lon[::skip_uv, ::skip_uv],
                      lat[::skip_uv, ::skip_uv],
                      u[::skip_uv, ::skip_uv],
                      v[::skip_uv, ::skip_uv],
                      scale=scale,
                      color=col,
                      width=width,
                      transform=ccrs.PlateCarree(), zorder=1)
    
    if ref_vector is not None:
        # add a reference vector
        loc_x=min(lon.ravel())
        loc_y=max(lat.ravel())
        ax.quiver(np.array([loc_x]), np.array([loc_y]),
                  np.array([ref_vector]), np.array([0]), 
                  scale=scale,
                  color=col,
                  width=width,
                  transform=ccrs.PlateCarree(), zorder=1,
                  )
        # add the label for it
        # try to define the y location of the text by shifting it down by a defined fraction of the y-range in the data
        d_y=max(lat.ravel())-min(lat.ravel())
        loc_y_txt=loc_y-d_y*0.04
        ax.text(loc_x, loc_y_txt,  str(ref_vector)+' m s$^{-1}$',
            ha='left', 
            #fontsize=time_font,
            transform=ccrs.PlateCarree())
    return uv_plt

def plot_isobaths(ax,fname,isobaths):
    h = post.get_grd_var(fname,'h')
    lon=h.lon_rho.values
    lat=h.lat_rho.values
    contour = ax.contour(lon,lat,h,isobaths,colors='k', 
               linestyles='dashed',
               linewidths=1,
               transform=ccrs.PlateCarree())
    # Label the contours with their respective isobaths values
    ax.clabel(contour, fmt='%1.0f', inline=True)

def get_uv_params(spd,scale_uv,ref_vector,skip_uv,num_vectors,aspect_ratio):
    # utility function to dynamically define params for scaling vectors    
    
    if skip_uv is None:
        # dynamically define how many vectors we want to plot
        [Ny,Nx]=np.shape(spd)
        skip_uv=int(Ny/num_vectors*aspect_ratio)
    
    max_spd=np.nanmax(spd.values.ravel())
    # Define thresholds for scaling and reference vector size
    # based on the maximum speed being plotted
    if max_spd < 0.25:
        scale = 2.5
        ref = 0.1
    elif max_spd < 0.5:
        scale = 5
        ref = 0.25
    elif max_spd < 0.75:
        scale = 7.5
        ref = 0.5
    elif max_spd < 1:
        scale = 10
        ref = 0.5
    else:
        scale = 15
        ref = 1
    
    # only update if they haven't already been specified by the user
    if scale_uv is None:
        scale_uv =scale
    if ref_vector is None:
        ref_vector =ref
        
    return scale_uv, skip_uv, ref_vector

def plot(fname,
        var='temp',
        grdname=None,
        tstep=0, # the step to plot, or the first step to animate. Can be an integer or a datetime object
        level=None, # this gets set to surface level if not provided
        ticks = np.linspace(12,22,num=11), # the ticks to plot
        cmap = 'Spectral_r',
        extents = None,
        ref_date = None, # datetime, from CROCO model setup
        cbar_loc = None, # [left, bottom, width, height]
        cbar_label = 'temperature ($\degree$C)',
        add_vectors = True,
        scale_uv = None,
        ref_vector = None,
        skip_uv = None,
        num_vectors=25, # baseline number of vectors in each direction, for a aspect ration of 1
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
    
    if grdname is None:
        grdname = fname
    
    # automatically set the figure size based on the grid aspect ratio
    lon = post.get_grd_var(grdname,'lon_rho')
    lat = post.get_grd_var(grdname,'lat_rho')
    d_lon = max(np.ravel(lon)) - min(np.ravel(lon))
    d_lat = max(np.ravel(lat)) - min(np.ravel(lat))
    aspect_ratio = d_lon/d_lat
    # set figsize according to the model aspect ratio
    if aspect_ratio>1:
        figsize = (6*aspect_ratio,6)
    else:
        figsize = (6,6/aspect_ratio)
    
    fig = plt.figure(figsize=figsize) 
    ax = plt.axes(projection=ccrs.Mercator())
    
    # if no level provided then default to the surface level
    if level is None:
        print('no vertical level provided - defaulting to surface layer')
        if isinstance(fname, xr.Dataset) or isinstance(fname, xr.DataArray):
            ds = fname.copy()
        else:
            ds = post.get_ds(fname,var)
        level = len(ds.s_rho) - 1
        ds.close()
    
    var_plt = plot_var(ax,fname,var,
             grdname=grdname,
             tstep=tstep,
             level=level,
             ref_date=ref_date, 
             ticks=ticks,
             extents=extents)
    
    time_plt = plot_time(ax,fname,tstep=tstep,ref_date=ref_date)
    
    plot_cbar(ax,var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc,aspect_ratio=aspect_ratio)
    
    if add_vectors:
        
        # dynamically update the vector scaling paramseters based on the actual veclocity data for this time-step
        u, v = post.get_uv(fname,grdname=grdname,tstep=tstep,level=level,ref_date=ref_date)
        spd = np.sqrt(u**2+v**2)
        scale_uv, skip_uv, ref_vector = get_uv_params(spd,scale_uv,ref_vector,skip_uv,num_vectors,aspect_ratio)
        
        uv_plt = plot_uv(ax,fname,
                      grdname=grdname,
                      tstep=tstep,
                      level=level,
                      ref_date=ref_date,
                      scale = scale_uv,
                      skip_uv = skip_uv,
                      ref_vector = ref_vector,
                      extents=extents
                      )
        
    if isobaths is not None:
        plot_isobaths(ax,grdname,isobaths)
    
    # write a jpg if specified
    if write_jpg:
        if jpg_out is None:
            jpg_out = fname.split('.nc')[0]+'_'+var+'_'+datetime.strftime(post.get_time(fname, ref_date)[tstep], '%Y%m%d_%H')+'.jpg'
        plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')
    
    # and/or write a gif if specified
    if write_gif: # do the animation
        # this only works if you specify tstep as an integer, not a datetime
        def plot_tstep(i):
            # get the data for this time-step
            time_i = post.get_time(fname, ref_date)[i]
            var_i = post.get_var(fname,var,grdname=grdname,tstep=i,level=level,ref_date=ref_date)
            var_i = var_i.values
            u_i, v_i = post.get_uv(fname,grdname=grdname,tstep=i,level=level,ref_date=ref_date)
            u_i = u_i.values
            v_i = v_i.values
            
            # update the figure for this time-step
            time_plt.set_text(datetime.strftime(time_i, '%Y-%m-%d %H:%M'))    
            var_plt.set_array(var_i.ravel())
            uv_plt.set_UVC(u_i[::skip_uv, ::skip_uv],
                                        v_i[::skip_uv, ::skip_uv])
        
        # animate
        if tstep_end is None: # if not defined then animate to end of file
            tstep_end=len(post.get_time(fname))
        
        anim = FuncAnimation(
            fig, plot_tstep, frames=range(tstep,tstep_end,skip_time)) 
        
        if gif_out is None:
            gif_out = fname.split('.nc')[0]+'_'+var+'.gif'
        anim.save(gif_out, writer='imagemagick')

def plot_blk(croco_grd, # the croco grid file - needed as not saved in the blk file
        croco_blk_file, # the croco blk file
        var='wspd',
        figsize=(6,6), # (hz,vt)
        tstep=0, # the step to plot (not going to worry about decoding the actual dates here)
        ref_date = datetime(1993,1,1), # datetime, from CROCO model setup
        ticks = [], # the ticks to plot
        cmap = 'Spectral_r',
        extents = None,#  [lon_min, lon_max, lat_min, lat_max]
        cbar_loc = None, #[0.7, 0.2, 0.02, 0.6],
        add_vectors = True,
        scale_uv = 150,
        skip_uv = 10,
        skip_time = 1, # every nth time-step will be animated (if provided)
        jpg_out=None,
        write_jpg=False,
        gif_out=None,
        write_gif=False,
        tstep_end=None, # The last timestep to animate. Only used if write_gif = True.
        ):
    '''
    this is a convenience function for doing a quick 2D plot 
    of a croco input blk file, just to make sure the forcing looks ok
    a quick plot on winds at the first time-step would be just
    plot_blk(croco_grd, croco_blk_file)
    '''
    # get the croco grid variables
    ds_croco_grd=xr.open_dataset(croco_grd)
    lon_rho=ds_croco_grd.lon_rho.values
    lat_rho=ds_croco_grd.lat_rho.values
    angle=ds_croco_grd.angle.values
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    ds_croco_grd.close()

    ds_blk = xr.open_dataset(croco_blk_file, decode_times=False)
    # get an array of datetimes for this file
    ds_blk_days = np.float64(ds_blk.bulk_time.values)
    blk_time = ref_date + timedelta(days = ds_blk_days[tstep])
    ds_blk_t = ds_blk.isel(bulk_time=tstep)

    var_data=ds_blk_t[var].values
    u = ds_blk_t.uwnd.values
    v = ds_blk_t.vwnd.values
    # get u,v onto rho grid
    u = post.u2rho(u)
    v = post.v2rho(v)
    # rotate to be east,north components
    u_rot = u * cos_a - v * sin_a
    v_rot = v * cos_a + u * sin_a

    fig = plt.figure() 
    ax = plt.axes(projection=ccrs.Mercator())

    setup_plot(ax, croco_grd)

    # set up the cmap to handle non-uniform input ticks
    if len(ticks)==0:
        ticks = np.linspace(min(np.ravel(var_data)),max(np.ravel(var_data)),num=20)
    n_levels = len(ticks)
    #vmin = min(ticks) 
    #vmax = max(ticks)
    levs = np.array(ticks)
    cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)

    # plot the data
    var_plt = ax.pcolormesh(lon_rho,
                              lat_rho,
                              var_data,
                              cmap=cmap,
                              norm=cmap_norm,
                              transform=ccrs.PlateCarree())

    plot_cbar(var_plt,label=var,ticks=ticks,loc=cbar_loc)

    if add_vectors:
        uv_plt = ax.quiver(lon_rho[::skip_uv, ::skip_uv],
                          lat_rho[::skip_uv, ::skip_uv],
                          u_rot[::skip_uv, ::skip_uv],
                          v_rot[::skip_uv, ::skip_uv],
                          scale=scale_uv,
                          transform=ccrs.PlateCarree(), zorder=1)
    
    # show tstep
    time_plt = ax.text(0.5, 1.01, 'tstep = ' + str(tstep)+': '+datetime.strftime(blk_time, '%Y-%m-%d %H:%M:%S'),
        ha='center', # fontsize=12,
        transform=ax.transAxes)
    
    # write a jpg if specified
    if write_jpg:
        if jpg_out is None:
            jpg_out = croco_blk_file.split('.nc')[0]+'_'+var+'_tstep'+str(tstep)+'.jpg'
        plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')
    
    # and/or write a gif if specified
    if write_gif: # do the animation
        # this only works if you stecify tstep as an integer, not a datetime
        def plot_tstep(i):
            # get the data for this time-step
            ds_blk_i = ds_blk.isel(bulk_time=i)
            var_i=ds_blk_i[var].values
            u = ds_blk_i.uwnd.values
            v = ds_blk_i.vwnd.values
            # get u,v onto rho grid
            u = post.u2rho(u)
            v = post.v2rho(v)
            # rotate to be east,north components
            u_i = u * cos_a - v * sin_a
            v_i = v * cos_a + u * sin_a
            
            # update the figure for this time-step
            blk_time = ref_date + timedelta(days = ds_blk_days[i])
            time_plt.set_text('tstep = '+str(i)+': '+datetime.strftime(blk_time, '%Y-%m-%d %H:%M:%S')) 
            var_plt.set_array(var_i.ravel())
            uv_plt.set_UVC(u_i[::skip_uv, ::skip_uv],
                                        v_i[::skip_uv, ::skip_uv])
        
        # animate
        if tstep_end is None: # if not defined then animate to end of file
            with xr.open_dataset(croco_blk_file) as ds: 
                tstep_end = len(ds.bulk_time)
        
        anim = FuncAnimation(
            fig, plot_tstep, frames=range(tstep,tstep_end,skip_time)) 
        
        if gif_out is None:
            gif_out = croco_blk_file.split('.nc')[0]+'_'+var+'.gif'
        anim.save(gif_out, writer='imagemagick')
        
if __name__ == "__main__":
    
    fname='../output/avg_Y2011M2.nc'
    plot(fname,
          write_gif=True,
          level = -100,
          ticks = np.linspace(10,24,num=15),
          tstep_end=20,
          cbar_loc = [0.85, 0.2, 0.02, 0.6])
    
        
