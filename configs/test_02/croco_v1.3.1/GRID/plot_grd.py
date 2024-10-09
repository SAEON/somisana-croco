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

# h = post.get_grd_var(fname, 'h')
lon_rho, lat_rho, mask_rho = post.get_lonlatmask(fname)
# h = h * mask_rho

figsize=(6,6) # (hz,vt)
extents = None

fig = plt.figure(figsize=figsize) 
ax = plt.axes(projection=ccrs.Mercator())

# Plot grid lines 
for i in range(lon_rho.shape[0]):
    plt.plot(lon_rho[i, :], lat_rho[i, :], color='grey', transform=ccrs.PlateCarree(), linewidth = 0.2)

for j in range(lon_rho.shape[1]):
    plt.plot(lon_rho[:, j], lat_rho[:, j], color='grey', transform=ccrs.PlateCarree(), linewidth = 0.2)

# set up the plot
crocplot.setup_plot(ax,fname,extents = extents)

# Coords for Salmiya
lat_ts = 29.356870
lon_ts = 48.111210
ax.scatter(lon_ts,lat_ts,20,color='k',transform=ccrs.PlateCarree())
time_plt = ax.text(lon_ts-0.1, lat_ts-0.1, 'Salmiya',
    ha='right', va='top', fontsize=10,
    transform=ccrs.PlateCarree())

jpg_out = 'plot_grd.jpg'
plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')

