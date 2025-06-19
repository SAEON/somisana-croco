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
import matplotlib.colors as mcolors

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

    # fallback in case cfeature.COLORS is broken
    if land_color is None or not mcolors.is_color_like(land_color):
        land_color = 'lightgray'  # or 'black', or any valid fallback

    land = LandmaskFeature(scale=lscale, facecolor=land_color, globe=globe)

    ax.add_feature(land, zorder=2,
                   facecolor=land_color,
                   edgecolor='black')

def setup_plot(ax, fname, extents=None, land_color=('k', 0), lscale='h'):
    '''
    generic stuff applicable to all 2D plots
    
    lscale = resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
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
             tstep=0,
             ref_date = datetime(2000, 1, 1, 0, 0, 0),
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
        time=slice(None), # see post.get_var() for 'time' format. If a single value, then a plot is made, if two values, then an animation between those times is made
        level=None, # see post.get_var() for 'level' format. Has to be a single value for this function to do a plot
        ticks = None, #np.linspace(12,22,num=11), (gets set automatically if None)
        cmap = 'Spectral_r',
        extents = None, # [lon0,lon1,lat0,lat1] whole domain plotted if None
        ref_date = None, # datetime, from CROCO model setup
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
    ds = post.get_var(fname,var,grdname=grdname,time=time,level=level,ref_date=ref_date)
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
        # get the range of the data to plot (using 5th and 95th percentiles)
        vmin=np.nanpercentile(da_var, 1)
        vmax=np.nanpercentile(da_var, 99)
        # round these to two significant figures
        vmin=round(vmin, 2 - int(np.floor(np.log10(abs(vmin)))) - 1)
        vmax=round(vmax, 2 - int(np.floor(np.log10(abs(vmax)))) - 1)
        num_ticks = 10
        step = (vmax - vmin) / num_ticks
        step = round(step, 2 - int(np.floor(np.log10(abs(step)))) - 1)
        # update vmax based on the rounded step
        vmax = vmin + num_ticks * step
        # Generate the ticks using the rounded step size
        ticks = np.arange(vmin, vmax + step/10, step) # Add a small value to ensure new_vmax is included
        
        # === NEW CODE START ===
        # Handle anomaly plots with symmetric color scale
        is_anomaly = var.endswith('_anom')
        if is_anomaly:
            cmap = 'bwr'
            abs_vals = np.abs(da_var.values).flatten()
            abs_vals = abs_vals[~np.isnan(abs_vals)]
            threshold = np.percentile(abs_vals, 95)
            vmin, vmax = -threshold, threshold
            ticks = np.linspace(vmin, vmax, 11)
        # === NEW CODE END ===
    
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
    extents = setup_plot(ax, grdname, extents = extents)
    
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
        ds_uv = post.get_uv(fname,grdname=grdname,time=time,level=level,ref_date=ref_date)
        
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
        
        # === NEW CODE START ===    
        # 🔹 NEW BLOCK: Auto-generate anomaly plot if plotting temp
        if var == 'temp' and 'temp_anom' in ds and gif_out is not None:
            gif_out_anom = gif_out.replace(".gif", "_anom.gif")  # ✅ more robust
            print(f'🔁 Also plotting anomaly variable to {gif_out_anom}')
            plot(
                fname=fname,
                var='temp_anom',
                level=level,
                cmap='bwr',
                ref_date=ref_date,
                gif_out=gif_out_anom,
                add_vectors=add_vectors,
                time=time,
                skip_time=skip_time,
                grdname=grdname,
                add_time_label=add_time_label,
                isobaths=isobaths
            )
            # === NEW CODE END ===
        
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
    
    fname='../output/avg_Y2010M12.nc'
    plot(fname)
    
        
