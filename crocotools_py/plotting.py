import os
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
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
    (from the OpenDrift code)
    lscale = resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
    """
    land = LandmaskFeature(scale=lscale, facecolor=land_color, globe=globe)

    ax.add_feature(land, zorder=2,
                   facecolor=land_color,
                   edgecolor='black')

def setup_plot(ax, fname, extents=None, add_land=True, land_color=('k', 0), lscale='h'):
    '''
    generic stuff applicable to all 2D plots
    
    lscale = resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
    '''
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
    if add_land:
        plot_land(ax,land_color=land_color,lscale=lscale)
        # ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='dimgrey', alpha=0.5, linestyle=':')
    gl.right_labels = False
    gl.top_labels = False
    
    return extents

def plot_var(ax,var_data,lon,lat,
             ticks = [], # the ticks to plot
             cmap = 'Spectral_r'
             ):
    '''
    Add a variable to a 2D plot
    '''
    
    # set up the cmap to handle non-uniform input ticks
    if len(ticks)==0:
        ticks = np.linspace(min(np.ravel(var_data)),max(np.ravel(var_data)),num=20)
    # n_levels = len(ticks)
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
             tick_font = 13,
             label='values',
             label_font=15,
             loc=None, # [left, bottom, width, height]
             aspect_ratio=1,
             orientation='vertical'):
    
    '''
    Add a colorbar to a plot
    '''
    if loc is None:
        # this can be hard coded because we took care to set up the figsize and
        # axis location to accomodate the colorbar
        x_position = 0.8
        x_thickness = 0.015
        loc = [x_position, 0.2, x_thickness, 0.6]
    
    cbarax = plt.gcf().add_axes(loc) 
    
    cbar_plt = plt.colorbar(var_plt, cbarax,
                        ticks=ticks,
                        orientation=orientation)
    cbar_plt.set_label(label, fontsize=label_font)
    cbar_plt.ax.tick_params(labelsize=tick_font)
    
    return cbar_plt

def plot_time(ax,time,
             loc=[0.5,1.01],
             time_fmt = '%Y-%m-%d %H:%M',
             time_font=15):
    '''
    Add time text to a 2D plot
    '''
    
    time_plt = ax.text(loc[0], loc[1],  pd.Timestamp(time).strftime(time_fmt),#datetime.strftime(time, time_fmt),
        ha='center', fontsize=time_font,
        transform=ax.transAxes)
    
    return time_plt

def plot_uv(ax,u,v,lon,lat,
              extents = None,
              skip_uv = 5,
              scale = 10,
              ref_vector = None,
              col = 'k'
              ):
    '''
    Add vectors to a 2D plot
    '''
        
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
        d_y=extents[3]-extents[2]
        d_x=extents[1]-extents[0]
        loc_x=extents[0]+d_x*0.04
        loc_y=extents[3]-d_y*0.04
        ax.quiver(np.array([loc_x]), np.array([loc_y]),
                  np.array([ref_vector]), np.array([0]), 
                  scale=scale,
                  color=col,
                  width=width,
                  transform=ccrs.PlateCarree(), zorder=1,
                  )
        # add the label for it
        # try to define the y location of the text by shifting it down by a defined fraction of the y-range in the data
        loc_y_txt=loc_y-d_y*0.04
        ax.text(loc_x, loc_y_txt,  str(ref_vector)+' m s$^{-1}$',
            ha='left', 
            #fontsize=time_font,
            transform=ccrs.PlateCarree())
    return uv_plt

def plot_isobaths(ax,fname,isobaths):
    h = post.get_grd_var(fname,'h')
    lon=post.get_grd_var(fname,'lon_rho')
    lat=post.get_grd_var(fname,'lat_rho')
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
    
    max_spd=np.nanmax(spd.ravel())
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
        ax=None, # allowing for adding to an existing axis
        var='temp', # croco variable to plot
        grdname=None, # option croco grid file (if grid variables arem't in the croco output file)
        time=slice(None), # see post.get_var() for 'time' format. If a single value, then a plot is made, if a slice then an animation between the define slice limits is made
        level=None, # see post.get_var() for 'level' format. Has to be a single value for this function to do a plot
        ticks = None, #np.linspace(12,22,num=11), (gets set automatically if None)
        cmap = 'Spectral_r',
        extents = None, # [lon0,lon1,lat0,lat1] whole domain plotted if None
        Yorig = None, # Origin year used in setting up CROCO model time
        add_cbar = True, # add a colorbar?
        cbar_loc = None, # [left, bottom, width, height] (gets set automatically if None)
        cbar_label = None, # 'temperature ($\degree$C)', we just use 'var' is None
        add_vectors = True, # add horizontal vectors?
        scale_uv = None, # define the horizontal vector scaling (gets set automatically if None)
        ref_vector = None, # value of reference vector (gets set automatically if None)
        skip_uv = None, # only every nth vector will get plotted (automatically defined in None)
        num_vectors=25, # baseline number of vectors in each direction, given an aspect ratio of 1
        skip_time = 1, # every nth time-step will be animated (if provided)
        add_time_label = True,
        add_land = True,
        isobaths = None, # optional list of isobaths to overlay over plot
        jpg_out=None, # full path to jpg output
        gif_out=None, # full path to gif output
        mp4_out=None, # option to rather write an mp4
        ):
    '''
    this is a convenience function for doing a quick 2D plot with minimal coding.
    this might also be used as example code for doing your own plots 
    there's also an option to turn the plot into an animation
    '''
    
    if grdname is None:
        grdname = fname
    
    # if no level provided then default to the surface level for the plot/animation
    if level is None:
        print('no vertical level provided - defaulting to surface layer')
        if isinstance(fname, xr.Dataset) or isinstance(fname, xr.DataArray):
            ds = fname.copy()
        else:
            ds = post.get_ds(fname,var)
        level = len(ds.s_rho) - 1
        ds.close()
        
    # get the data we want to 
    print('extracting the data to plot')
    ds = post.get_var(fname,var,grdname=grdname,time=time,level=level,Yorig=Yorig)
    da_var=ds[var]
    time_var=np.atleast_1d(ds.time.values)
    lon = post.get_grd_var(grdname,'lon_rho').values
    lat = post.get_grd_var(grdname,'lat_rho').values
    
    
    if len(time_var)==1:
        data_plt=da_var.values
    else:
        # this will be an animation, starting with the first time-step
        data_plt=da_var.isel(time=0).values
    
    if ticks is None:
    
        is_anomaly = var.endswith('_anom')
        num_ticks = 10
    
        if is_anomaly:
            vmax = np.nanpercentile(abs(da_var), 99)
            # Symmetric color scale around 0
            vmax = round(vmax, 2 - int(np.floor(np.log10(abs(vmax)))) - 1)
            vmin = -vmax
            cmap = 'bwr'
            #isobaths=[200,500],
            cbar_label = 'temperature anomaly ($\degree$C)'
        else:
            vmin = np.nanpercentile(da_var, 1)
            vmax = np.nanpercentile(da_var, 99)
            # Round to two significant figures
            vmin = round(vmin, 2 - int(np.floor(np.log10(abs(vmin)))) - 1)
            vmax = round(vmax, 2 - int(np.floor(np.log10(abs(vmax)))) - 1)
    
        # Shared logic: step rounding and tick generation
        step = (vmax - vmin) / num_ticks
        step = round(step, 2 - int(np.floor(np.log10(abs(step)))) - 1)
        vmax = vmin + num_ticks * step
        ticks = np.arange(vmin, vmax + step / num_ticks, step)
    
    # compute the extents from the grid if not explicitly defined
    if extents is None:
        lon_min = min(np.ravel(lon))
        lon_max = max(np.ravel(lon))
        lat_min = min(np.ravel(lat))
        lat_max = max(np.ravel(lat))
        factor=0.05 # factor of domain size used to get dl
        dl = 0.5 * (lon_max - lon_min + lat_max - lat_min) * factor
        extents=[lon_min-dl,lon_max+dl,lat_min-dl,lat_max+dl]
    
    # Create a Mercator projection
    proj = ccrs.Mercator()
    # Convert the corner points of the extents to Mercator projected coordinates
    # this is needed to compute the plot aspect ratio properly
    x0, y0 = proj.transform_point(extents[0], extents[2], ccrs.PlateCarree())  # lon0, lat0
    x1, y1 = proj.transform_point(extents[1], extents[3], ccrs.PlateCarree())  # lon1, lat1
    # Calculate the width and height of the domain in projected coordinates
    width = abs(x1 - x0)
    height = abs(y1 - y0)
    aspect_ratio = width/height
    
    if ax is None:        
        # set figsize according to the plot aspect ratio
        if aspect_ratio>1:
            fig_width = 6*aspect_ratio
            fig_height = 6
        else:
            fig_width = 6
            fig_height = 6/aspect_ratio
        # cbar_ax_width = 0.2 * fig_width if add_cbar else 0
        buffer_left=0.1 * fig_width
        buffer_right=0.2 * fig_width if add_cbar else 0.1 * fig_width
        figsize = (buffer_left + fig_width + buffer_right, fig_height)
        
        fig = plt.figure(figsize=figsize) 
        # [left, bottom, width, height] in fractions of figure dimensions
        width=0.7 if add_cbar else 0.8
        ax = fig.add_axes([0.1, 0.1, width, 0.8], projection=proj)
    
    # set up the plot
    extents = setup_plot(ax, grdname, extents = extents, add_land=add_land)
    
    # plot the data
    var_plt = plot_var(ax,data_plt,lon,lat, 
             ticks=ticks,
             cmap=cmap
             )
    
    # add a time label
    if add_time_label:
        time_plt = plot_time(ax,time_var[0])
    
    if add_cbar:
        if cbar_label is None:
            cbar_label = var
        plot_cbar(ax,var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc,aspect_ratio=aspect_ratio)
    
    if add_vectors:
        
        print('getting the u/v vectors')
        
        # dynamically update the vector scaling paramseters based on the actual veclocity data
        ds_uv = post.get_uv(fname,grdname=grdname,time=time,level=level,Yorig=Yorig)
        
        da_u=ds_uv.u
        da_v=ds_uv.v
        
        if len(time_var)==1:
            u_plt=da_u.values
            v_plt=da_v.values
        else:
            # this will be an animation, starting with the first time-step
            u_plt=da_u.isel(time=0).values
            v_plt=da_v.isel(time=0).values
        
        spd_plt = np.sqrt(u_plt**2+v_plt**2)
        
        # subset spd based on the plot extents before running get_uv_params()
        # I'm commenting this as it doesn't work for very curvilinear grids
        # the user always has the option of specifying the size of the ref_vector
        # eta_rho,_,xi_rho,_=post.domain_to_slice(slice(None),slice(None),slice(None),slice(None),extents,grdname,'temp')
        # spd=spd[eta_rho,xi_rho]
        
        scale_uv, skip_uv, ref_vector = get_uv_params(spd_plt,scale_uv,ref_vector,skip_uv,num_vectors,aspect_ratio)
        
        uv_plt = plot_uv(ax,u_plt,v_plt,lon,lat,
                      scale = scale_uv,
                      skip_uv = skip_uv,
                      ref_vector = ref_vector,
                      extents=extents
                      )
        
    if isobaths is not None:
        plot_isobaths(ax,grdname,isobaths)
    
    # write a jpg if specified
    if len(time_var)==1: # single time-step specified, so a plot, not an animation
        if jpg_out is not None:
            print('writing '+jpg_out)
            plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')
    
    else: # do the animation
        def plot_tstep(i):
            # get the data for this time-step
            var_i=da_var.isel(time=i).values
            
            # update the figure for this time-step
            if add_time_label:
                time_plt.set_text(pd.Timestamp(time_var[i]).strftime('%Y-%m-%d %H:%M'))
                var_plt.set_array(var_i.ravel())
            
            if add_vectors:
                u_i=da_u.isel(time=i).values
                v_i=da_v.isel(time=i).values
                uv_plt.set_UVC(u_i[::skip_uv, ::skip_uv],
                                        v_i[::skip_uv, ::skip_uv])
        
        # animate
        print('making animation')
        tstep_end=len(time_var)
        
        anim = FuncAnimation(
            fig, plot_tstep, frames=range(0,tstep_end,skip_time)) 
        
        if gif_out is not None:
            print('writing '+gif_out)
            anim.save(gif_out, writer='imagemagick')
        if mp4_out is not None:
            print('writing '+mp4_out)
            anim.save(mp4_out, writer="ffmpeg")

def plot_blk(croco_grd, # the croco grid file - needed as not saved in the blk file
        croco_blk_file, # the croco blk file
        var='wspd',
        figsize=(6,6), # (hz,vt)
        tstep=0, # the step to plot (not going to worry about decoding the actual dates here)
        Yorig = 1993, # from CROCO model setup
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
    blk_time = datetime(Yorig,1,1) + timedelta(days = ds_blk_days[tstep])
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
            blk_time = datetime(Yorig,1,1) + timedelta(days = ds_blk_days[i])
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
        

# MHW/MCS PLOTTING (TIME-SERIES, FLAG MAPS, INTENSITY GIF)

# sites of interest
TARGETS = {
    "Kleinsee":       (17.030382, -29.680623),
    "Hondeklipbaai":  (17.252461, -30.315292),
    "Doringbaai":     (18.213554, -31.814509),
    "Elandsbaai":     (18.30165,  -32.312317),
    "Laaiplek":       (18.125354, -32.742041),
    "Paternoster":    (17.870305, -32.777566),
    "Saldanha":       (17.929861, -33.074807),
    "Yzerfontein":    (18.13382,  -33.361876),
    "Bloubergstrand": (18.443896, -33.803906),
    "Oudekraal":      (18.342541, -33.980098),
    "Cape Point":     (18.46024,  -34.358313),
    "Simonstown":     (18.442294, -34.176514),
    "Strand":         (18.810174, -34.120553),
    "Hangklip":       (18.803882, -34.374716),
    "Kleinmond":      (19.026591, -34.355882),
    "Hermanus":       (19.256989, -34.425957),
    "Gansbaai":       (19.323381, -34.576985),
}

WINDOW_DAYS = 10

# Hobday category Colors

FILL_MOD   = "#ffc73e";  FILL_STR   = "#f77819"
FILL_SEV   = "#bf460c";  FILL_EXT   = "#4e1909"
FILL_C_MOD = "#a6d3e8";  FILL_C_STR = "#5da6c9"
FILL_C_SEV = "#2074a3";  FILL_C_EXT = "#103c68"

MHW_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_MOD,   2: FILL_STR,   3: FILL_SEV,   4: FILL_EXT}
MCS_FLAG_COLOURS = {0: "#4CAF7D", 1: FILL_C_MOD, 2: FILL_C_STR, 3: FILL_C_SEV, 4: FILL_C_EXT}

CMAP_9 = mplc.ListedColormap([FILL_C_EXT, FILL_C_SEV, FILL_C_STR, FILL_C_MOD,"#ffffff",
    FILL_MOD, FILL_STR, FILL_SEV, FILL_EXT
])
CMAP_9.set_bad("white")
BNORM_9 = mplc.BoundaryNorm([-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5], CMAP_9.N)

# --- Helpers ---
def nearest(lon2d, lat2d, lon0, lat0):
    d2 = (lon2d - lon0)**2 + (lat2d - lat0)**2
    return np.unravel_index(np.argmin(d2), d2.shape)

def doy_index(times):
    return np.clip(pd.to_datetime(times).dayofyear.values - 1, 0, 365)

def compute_site_flag_data(sites, cat_ds, lev):
    site_data = {}
    for site_name, data in sites.items():
        pj, pi = data["pj"], data["pi"]
        cat = (cat_ds["category"]
               .isel(s_rho=lev, eta_rho=pj, xi_rho=pi)
               .load().values.astype(float))

        mhw_days = cat[cat > 0]
        mcs_days = np.abs(cat[cat < 0])

        avg_mhw   = float(np.mean(mhw_days)) if len(mhw_days) > 0 else 0.0
        avg_mcs   = float(np.mean(mcs_days)) if len(mcs_days) > 0 else 0.0
        score_mhw = avg_mhw * len(mhw_days)
        score_mcs = avg_mcs * len(mcs_days)

        if score_mhw >= score_mcs:
            site_data[site_name] = {"mode": "MHW", "avg_cat": avg_mhw}
        else:
            site_data[site_name] = {"mode": "MCS", "avg_cat": avg_mcs}
    return site_data


def plot_timeseries_multisite(sites, today, output_dir, depth_name):
    out_dir = Path(output_dir)
    today   = pd.Timestamp(today)

    for site_name, data in sites.items():
        fct_dates = pd.DatetimeIndex(data["fct_dates"])
        fct_temp  = np.asarray(data["fct_temp"])
        fct_seas  = np.asarray(data["fct_seas"])
        fct_h_thr = np.asarray(data["fct_h_thr"])
        fct_c_thr = np.asarray(data["fct_c_thr"])

        obs_dates = data.get("obs_dates")
        obs_temp  = data.get("obs_temp")
        obs_seas  = data.get("obs_seas")
        obs_h_thr = data.get("obs_h_thr")
        obs_c_thr = data.get("obs_c_thr")

        if obs_dates is not None and obs_temp is not None:
            cut   = np.asarray(obs_dates) < np.asarray(fct_dates[0])
            n_cut = int(cut.sum())
            def _cat(obs_arr, fct_arr):
                base = np.asarray(obs_arr)[cut] if obs_arr is not None else np.full(n_cut, np.nan)
                return np.concatenate([base, np.asarray(fct_arr)])
            all_dates = pd.DatetimeIndex(list(np.asarray(obs_dates)[cut]) + list(fct_dates))
            all_temp  = _cat(obs_temp,  fct_temp)
            all_seas  = _cat(obs_seas,  fct_seas)
            all_h_thr = _cat(obs_h_thr, fct_h_thr)
            all_c_thr = _cat(obs_c_thr, fct_c_thr)
        else:
            all_dates, all_temp, all_seas, all_h_thr, all_c_thr = fct_dates, fct_temp, fct_seas, fct_h_thr, fct_c_thr

        if not np.isfinite(all_temp).any():
            continue

        t_min, t_max = today - pd.Timedelta(days=int(WINDOW_DAYS * 0.5)), today + pd.Timedelta(days=int(WINDOW_DAYS * 0.5))
        mask  = (all_dates >= t_min) & (all_dates <= t_max)
        d, T, S, H, C = all_dates[mask], all_temp[mask], all_seas[mask], all_h_thr[mask], all_c_thr[mask]
        obs_m, fct_m = d <= today, d >= today
        dh, dc = np.maximum(H - S, 1e-6), np.maximum(S - C, 1e-6)

        fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
        ax.yaxis.grid(True, color="#cccccc", linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)

        mhw_fills = [(T > H, H, np.minimum(T, H+dh), FILL_MOD), (T > H+dh, H+dh, np.minimum(T, H+2*dh), FILL_STR),
                     (T > H+2*dh, H+2*dh, np.minimum(T, H+3*dh), FILL_SEV), (T > H+3*dh, H+3*dh, T, FILL_EXT)]
        mcs_fills = [(T < C, np.maximum(T, C-dc), C, FILL_C_MOD), (T < C-dc, np.maximum(T, C-2*dc), C-dc, FILL_C_STR),
                     (T < C-2*dc, np.maximum(T, C-3*dc), C-2*dc, FILL_C_SEV), (T < C-3*dc, T, C-3*dc, FILL_C_EXT)]
        
        for where, lo, hi, col in mhw_fills + mcs_fills:
            ax.fill_between(d, lo, hi, where=where, color=col, alpha=0.85, interpolate=True, zorder=1)

        ax.plot(d, H, color="#d62728", ls="--", lw=1.6, zorder=4)
        ax.plot(d, C, color="#1f77b4", ls="--", lw=1.6, zorder=4)
        ax.plot(d, S, color="#7f7f7f", ls=":",  lw=1.4, zorder=4)
        ax.plot(d[obs_m], T[obs_m], color="#555555", lw=2.2, zorder=5, solid_capstyle="round")
        if fct_m.any(): ax.plot(d[fct_m], T[fct_m], color="black", lw=2.2, zorder=5, solid_capstyle="round")

        ax.axvline(today, color="black", lw=1.2, zorder=6)
        ax.text(today + pd.Timedelta(hours=6), ax.get_ylim()[1], " Today", va="top", ha="left", fontsize=9)

        lon_v, lat_v = data["lon"], data["lat"]
        ax.set_title(f"{site_name}  ({abs(lat_v):.3f}°{'S' if lat_v < 0 else 'N'}, {lon_v:.3f}°{'W' if lon_v < 0 else 'E'})", fontsize=14, fontweight="bold", pad=10)
        ax.set_ylabel(f"Temperature [°C]", fontsize=11)
        ax.set_xlim(t_min, t_max)
        for spine in ("top", "right"): ax.spines[spine].set_visible(False)

        handles = [Line2D([0],[0], color="#555555", lw=2.2, label="SST observed"), Line2D([0],[0], color="black", lw=2.2, label="SST forecast"),
                   Line2D([0],[0], color="#d62728", lw=1.6, ls="--", label="MHW threshold"), Line2D([0],[0], color="#1f77b4", lw=1.6, ls="--", label="MCS threshold"),
                   Line2D([0],[0], color="#7f7f7f", lw=1.4, ls=":", label="Climatology")]
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=9, frameon=False)
        
        plt.tight_layout()
        plt.savefig(out_dir / f"{site_name.replace(' ', '_')}_{depth_name}_{today.strftime('%Y%m%d')}.png", dpi=150, bbox_inches="tight")
        plt.close()

def plot_flag_map(site_data, today, start_date, end_date, out_path, lat, lon, depth_name="Surface"):
    out_path = Path(out_path)
    def _flag_col(mode, cat):
        c = max(0, min(4, int(round(cat if pd.notna(cat) else 0))))
        return (MHW_FLAG_COLOURS if mode == "MHW" else MCS_FLAG_COLOURS)[c]

    coast_order = [
    "Kleinsee",
    "Hondeklipbaai",
    "Doringbaai",
    "Elandsbaai",
    "Laaiplek",
    "Paternoster",
    "Saldanha",
    "Yzerfontein",
    "Bloubergstrand",
    "Oudekraal",
    "Cape Point",
    "Simonstown",
    "Strand",
    "Hangklip",
    "Kleinmond",
    "Hermanus",
    "Gansbaai"]
    
    BOX_SIZE, BOX_STEP, OFFSHORE, all_boxes = 0.75, 0.50, -0.10, []
    
    for k in range(len(coast_order) - 1):
        lon0, lat0 = TARGETS[coast_order[k]]; lon1, lat1 = TARGETS[coast_order[k + 1]]
        info = site_data.get(coast_order[k], {"mode": "MHW", "avg_cat": 0})
        seg = np.hypot(lon1 - lon0, lat1 - lat0)
        px, py = -(lat1 - lat0) / seg, (lon1 - lon0) / seg 
        for j in range(max(2, int(seg / BOX_STEP))):
            t = j / max(2, int(seg / BOX_STEP))
            all_boxes.append((lon0 + t*(lon1-lon0) + px*OFFSHORE, lat0 + t*(lat1-lat0) + py*OFFSHORE, _flag_col(info["mode"], info["avg_cat"])))
            
    lon_l, lat_l = TARGETS[coast_order[-1]]; lon_p, lat_p = TARGETS[coast_order[-2]]
    seg = np.hypot(lon_l - lon_p, lat_l - lat_p); px, py = -(lat_l - lat_p)/seg, (lon_l - lon_p)/seg
    all_boxes.append((lon_l + px*OFFSHORE, lat_l + py*OFFSHORE, _flag_col(site_data.get(coast_order[-1], {"mode": "MHW", "avg_cat": 0})["mode"], site_data.get(coast_order[-1], {"mode": "MHW", "avg_cat": 0})["avg_cat"])))

    fig = plt.figure(figsize=(10, 13), dpi=150)
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.set_extent([15.8, 20.5, -36.0, -28.0], crs=ccrs.PlateCarree())

    for cx, cy, col in all_boxes:
        ax.add_patch(FancyBboxPatch((cx - BOX_SIZE/2, cy - BOX_SIZE/2), BOX_SIZE, BOX_SIZE, boxstyle="round,pad=0.04", facecolor=col, edgecolor="white", linewidth=0.6, zorder=3, transform=ccrs.PlateCarree()))

    for site_name, (site_lon, site_lat) in TARGETS.items():
        info = site_data.get(site_name, {"mode": "MHW", "avg_cat": 0})
        cat_int = max(0, min(4, int(round(info["avg_cat"]))))
        ax.plot(site_lon, site_lat, "o", ms=4, color="white", zorder=8, mec="black", mew=0.8, transform=ccrs.PlateCarree())
        ax.text(site_lon + 0.08, site_lat, f"{site_name}\n{info['mode']} – {['None', 'Moderate', 'Strong', 'Severe', 'Extreme'][cat_int]}", ha="left", va="center", fontsize=6.5, fontweight="bold", color="#1a3a5c", zorder=9, transform=ccrs.PlateCarree(), path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=4); ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor="#555544", zorder=5)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="#aaaaaa", alpha=0.8, linestyle="--", zorder=2)
    gl.top_labels = gl.right_labels = False

    ax.set_title(f"SA West Coast  ·  MHW / MCS Flag Map  ·  {depth_name}\nForecast: {pd.to_datetime(start_date).strftime('%d %b')} – {pd.to_datetime(end_date).strftime('%d %b %Y')}", fontsize=12, color="#1a3a5c", pad=8)
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()

def _update_spatial_frame(frame, cat_data, time_data, mesh_obj, title_obj, d_name):
    mesh_obj.set_array(cat_data[frame].ravel())
    title_obj.set_text(f"MHW & MCS Categories ({d_name})\nDate: {str(time_data[frame])[:10]}")
    return mesh_obj, title_obj

def animate_spatial_categories(cat_ds, ds_fcst, lat, lon, depth_name, lev, is_varying, idx_2d, out_path):
    out_path = Path(out_path)
    times = pd.to_datetime(cat_ds.time.values)

    if is_varying:
        nj, ni = idx_2d.shape
        jj, ii = np.meshgrid(np.arange(nj), np.arange(ni), indexing="ij")
        valid  = idx_2d >= 0
        idx_cl = np.where(valid, idx_2d, 0)
        cat = cat_ds["category"].values[:, idx_cl, jj, ii].astype(float)
        cat[:, ~valid] = np.nan
    else:
        cat = cat_ds["category"].isel(s_rho=lev).values.astype(float)

    mask = ds_fcst["mask_rho"].values if "mask_rho" in ds_fcst else np.ones_like(lat)
    if mask.ndim > 2: mask = mask[0]
    cat = np.where(mask[np.newaxis, :, :] == 1, cat, np.nan)

    fig = plt.figure(figsize=(9, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(lon, lat, cat[0], transform=ccrs.PlateCarree(), cmap=CMAP_9, norm=BNORM_9, shading="auto")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.LAND, facecolor="white", edgecolor='black', zorder=2)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=":")

    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.04, ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4])
    cbar.ax.set_yticklabels(["Ext", "Sev", "Str", "Mod", "Neut", "Mod", "Str", "Sev", "Ext"])
    cbar.set_label("MCS (Cold)  ←  Intensity  →  MHW (Heat)")

    title = ax.set_title(f"MHW & MCS Categories ({depth_name})\nDate: {str(times[0])[:10]}")
    ani = FuncAnimation(fig, _update_spatial_frame, frames=len(times), fargs=(cat, times, mesh, title, depth_name), blit=False)
    ani.save(out_path, writer='ffmpeg', fps=5, dpi=120)
    plt.close(fig)

# --- Master Wrapper Function ---
def plot_operational_mhw_mcs(forecast_file, cat_file, clim_file, out_dir, start_date, end_date, Yorig=2000):
    """
    Operational entry point to run time series, flag maps, and GIF animations.
    """
    print("\n=== Rendering Operational MHW/MCS Visuals ===")
    out_dir = Path(out_dir)
    ds_clim = xr.open_dataset(clim_file)
    ds_cat  = xr.open_dataset(cat_file)
    ds_fcst = post.handle_time(post.get_ds(forecast_file, "temp"), Yorig=Yorig)
    
    lat = ds_fcst.lat_rho.values if "lat_rho" in ds_fcst else ds_fcst.lat.values
    lon = ds_fcst.lon_rho.values if "lon_rho" in ds_fcst else ds_fcst.lon.values
    if lat.ndim > 2: lat, lon = lat[0], lon[0]
    h = ds_fcst.h.values if "h" in ds_fcst else np.zeros_like(lat)
    if h.ndim > 2: h = h[0]
    nlev = len(ds_fcst.s_rho) if "s_rho" in ds_fcst else ds_fcst.dims.get("s_rho", 32)
    # today = pd.Timestamp(ds_fcst.time.values[0]).normalize()
    today = pd.Timestamp.now().normalize()

    Cs_r = ds_fcst.Cs_r.values if "Cs_r" in ds_fcst else np.linspace(-1, 0, nlev)
    sc_r = ds_fcst.sc_r.values if "sc_r" in ds_fcst else np.linspace(-1, 0, nlev)
    hc   = float(ds_fcst.hc.values) if "hc" in ds_fcst else 10.0
    z = np.zeros((nlev, *h.shape))
    for k in range(nlev): z[k] = ((hc * sc_r[k] + h * Cs_r[k]) / (hc + h)) * h
    idx_100m = np.argmin(np.abs(z - (-100.0)), axis=0)
    idx_100m[h < 100] = -1

    depth_levels = {
        "Surface": {"type": "fixed",   "lev": nlev - 1},
        "100m":    {"type": "varying", "lev": idx_100m},
        "Bottom":  {"type": "fixed",   "lev": 0},
    }

    ds_fcst_single = post.handle_time(xr.open_dataset(forecast_file, decode_times=False), Yorig=Yorig)

    for depth_name, depth_info in depth_levels.items():
        print(f"\nProcessing Depth: {depth_name}...")
        sites = {}
        for site_name, (site_lon, site_lat) in TARGETS.items():
            pj, pi = nearest(lon, lat, site_lon, site_lat)
            lev_site = depth_info["lev"] if depth_info["type"] == "fixed" else int(depth_info["lev"][pj, pi])
            if lev_site < 0: continue

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
        if depth_info["type"] == "fixed":
            lev_flag = depth_info["lev"] 
        else:
            # Get all depth values for sites
            depth_values = [depth_info["lev"][d["pj"], d["pi"]] for d in sites.values()]
            # Filter out NaNs and negative values
            valid_depths = [int(v) for v in depth_values if not np.isnan(v) and v >= 0]
    
            if valid_depths:
                lev_flag = int(np.median(valid_depths))
            else:
                # Fallback if NO sites are deep enough (e.g., all sites are < 100m)
                lev_flag = 0
                
        plot_flag_map(compute_site_flag_data(sites, ds_cat, lev_flag), today, start_date, end_date, out_dir / f"FlagMap_{depth_name}_{today.strftime('%Y%m%d')}.png", lat, lon, depth_name)

        print("  -> Spatial GIF...")
        varying = depth_info["type"] == "varying"
        animate_spatial_categories(ds_cat, ds_fcst, lat, lon, depth_name, depth_info["lev"], varying, depth_info["lev"] if varying else None, out_dir / f"Categories_Animation_{depth_name}.mp4")

    ds_fcst_single.close(); ds_fcst.close(); ds_clim.close(); ds_cat.close()
    print(f"\nAll visuals saved to: {out_dir}")
        
if __name__ == "__main__":
    
    fname='../output/avg_Y2010M12.nc'
    plot(fname)
    
        
