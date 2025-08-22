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

figsize=(8,5) # (hz,vt)
# ticks = np.linspace(0,100,num=11) # the ticks to plot
ticks = [0,50,100,200,300,500,1000,2000,3000,4000,5000]
cmap = cmo.deep
extents = None
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

crocplot.plot_cbar(ax,var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc)

jpg_out = 'plot_bathy.jpg'
plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')

