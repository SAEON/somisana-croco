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

fname='croco_grd.nc'

h = post.get_grd_var(fname, 'h')
lon_rho, lat_rho, mask_rho = post.get_lonlatmask(fname)
h = h * mask_rho

figsize=(6,6) # (hz,vt)
# ticks = np.linspace(0,100,num=11) # the ticks to plot
ticks = [0,10,20,30,40,50,60,70,80,90,100,200,300]
cmap = cmo.deep
extents = []
cbar_loc = [0.92, 0.2, 0.02, 0.6]
cbar_label = 'bathymetry (m)'

fig = plt.figure(figsize=figsize) 
ax = plt.axes(projection=ccrs.Mercator())

# set up the plot
crocplot.setup_plot(ax,fname,extents = extents)

levs = np.array(ticks)
cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)

# plot the data
var_plt = ax.pcolormesh(lon_rho,
                          lat_rho,
                          h,
                          cmap=cmap,
                          norm=cmap_norm,
                          transform=ccrs.PlateCarree())

crocplot.plot_cbar(var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc)

# Coords for Salmiya
lat_ts = 29.356870
lon_ts = 48.111210
ax.scatter(lon_ts,lat_ts,20,color='k',transform=ccrs.PlateCarree())
time_plt = ax.text(lon_ts-0.1, lat_ts-0.1, 'Salmiya',
    ha='right', va='top', fontsize=10,
    transform=ccrs.PlateCarree())

jpg_out = 'plot_bathy.jpg'
plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')

