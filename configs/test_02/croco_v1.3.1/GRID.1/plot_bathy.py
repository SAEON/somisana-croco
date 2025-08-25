import crocotools_py.postprocess as post
import crocotools_py.plotting as crocplot
import crocotools_py.validation as val
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo

fname='../GRID/croco_grd.nc'

h_p = post.get_grd_var(fname, 'h')
lon_rho_p, lat_rho_p, mask_rho_p = post.get_lonlatmask(fname)
h_p = h_p * mask_rho_p

fname='croco_grd.nc.1'

h = post.get_grd_var(fname, 'h')
lon_rho, lat_rho, mask_rho = post.get_lonlatmask(fname)
h = h * mask_rho

figsize=(6,6) # (hz,vt)
# ticks = np.linspace(0,100,num=11) # the ticks to plot
ticks = [0,5,10,15,20,25,30,35,40,45,50]
cmap = cmo.deep
extents = [52.3,54.6,23.7,25.5]
cbar_loc = [0.92, 0.2, 0.02, 0.6]
cbar_label = 'bathymetry (m)'

fig = plt.figure(figsize=figsize) 
ax = plt.axes(projection=ccrs.Mercator())

# set up the plot
crocplot.setup_plot(ax,fname,extents = extents)

levs = np.array(ticks)
cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)

# plot the data
var_plt_p = ax.pcolormesh(lon_rho_p,
                          lat_rho_p,
                          h_p,
                          cmap=cmap,
                          norm=cmap_norm,
                          transform=ccrs.PlateCarree())

var_plt = ax.pcolormesh(lon_rho,
                          lat_rho,
                          h,
                          cmap=cmap,
                          norm=cmap_norm,
                          transform=ccrs.PlateCarree())

# crocplot.plot_cbar(var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc)
crocplot.plot_cbar(ax,var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc)

# Coords for Mirfa
lat_ts = 24.122461
lon_ts = 53.444186
ax.scatter(lon_ts,lat_ts,20,color='k',transform=ccrs.PlateCarree())
time_plt = ax.text(lon_ts-0.1, lat_ts-0.1, 'Salmiya',
    ha='right', va='top', fontsize=10,
    transform=ccrs.PlateCarree())

jpg_out = 'plot_bathy.jpg'
plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')

