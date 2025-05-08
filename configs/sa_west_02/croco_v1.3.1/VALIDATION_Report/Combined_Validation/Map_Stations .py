#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:05:34 2024
@author: nkululeko
"""

import os
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# Directory containing the NetCDF files
netcdf_dir = '/media/disk02/DATA-20231010T133411Z-003/DATA/AJSMIT_UTR/Proccessed/NetCDF_Stations'
savepath = '/home/nc.memela/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/'

# Initialize a list to store station data
stations = []

# Iterate over all NetCDF files in the directory
for filename in os.listdir(netcdf_dir):
    if filename.endswith('.nc'):
        filepath = os.path.join(netcdf_dir, filename)
        # Open the NetCDF file
        ds = xr.open_dataset(filepath)
        # Extract the station, lat, and lon
        station = ds['station'].values[0]
        lat = ds['lat'].values[0]
        lon = ds['lon'].values[0]
        # Append the data to the list
        stations.append({'station': station, 'lat': lat, 'lon': lon})

# Create a DataFrame from the collected station data
stations_df = pd.DataFrame(stations)

# Define map extent (adjust as needed)
extent = [17, 19.5, -34.5, -31.5]  # [min_lon, max_lon, min_lat, max_lat]

# Define zoomed-in extent
# [min_lon, max_lon, min_lat, max_lat]
zoom_extent = [18.25, 18.55, -34.4, -33.9]

# Define unique markers
markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'H',
           'd', '<', '>', 'P', 'X', 'x']  # Changed '|' to 's'

# Variable for marker size (adjust as needed)
marker_size = 11

# Plot the main map
fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(extent, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.RIVERS)
ax.add_feature(cfeature.LAKES)

# Plot each station as a scatter point with unique markers on the main map
for station, marker in zip(stations_df['station'].unique(), markers):
    station_data = stations_df[stations_df['station'] == station]
    # Make '*' and '^' markers red
    color = 'red' if marker == '*' or marker == '^' else 'k'
    plt.scatter(station_data['lon'], station_data['lat'], marker=marker,
                color=color, s=marker_size**2, transform=ccrs.PlateCarree(), label=station)

# Set plot title and labels
plt.title('Map of Validation UTR and surface Insitu Stations')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Add gridlines and labels to main map
gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
gl.top_labels = False
gl.right_labels = False
gl.xlines = True
gl.ylines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([17, 18, 19, 20])
gl.ylocator = mticker.FixedLocator([-34.5, -34, -33.5, -33, -32.5, -32, -31.5])

# Create legend with unique markers for each station
legend_elements = [Line2D([0], [0], marker=marker, color='k' if marker != '*' and marker != '^' else 'red', label=station, markersize=marker_size, linestyle='None')
                   for station, marker in zip(stations_df['station'].unique(), markers)]
plt.legend(handles=legend_elements, loc='upper right', fontsize='small')

# Create inset map for zoomed-in region
ax_inset = fig.add_axes([0.55, 0.28, 0.3, 0.33], projection=ccrs.PlateCarree())
ax_inset.set_extent(zoom_extent, crs=ccrs.PlateCarree())
ax_inset.add_feature(cfeature.COASTLINE)
ax_inset.add_feature(cfeature.BORDERS, linestyle=':')
ax_inset.add_feature(cfeature.LAND)
ax_inset.add_feature(cfeature.OCEAN)
ax_inset.add_feature(cfeature.RIVERS)
ax_inset.add_feature(cfeature.LAKES)

# Plot stations in the zoomed-in region
for station, marker in zip(stations_df['station'].unique(), markers):
    station_data = stations_df[stations_df['station'] == station]
    # Make '*' and '^' markers red
    color = 'red' if marker == '*' or marker == '^' else 'k'
    ax_inset.scatter(station_data['lon'], station_data['lat'], marker=marker,
                     color=color, s=marker_size**2, transform=ccrs.PlateCarree())

# Add gridlines and labels to inset map
gl_inset = ax_inset.gridlines(draw_labels=True, linestyle='--', color='gray')
gl_inset.top_labels = False
gl_inset.right_labels = False
gl_inset.xlines = True
gl_inset.ylines = True
gl_inset.xformatter = LONGITUDE_FORMATTER
gl_inset.yformatter = LATITUDE_FORMATTER
gl_inset.xlocator = mticker.FixedLocator([18.3, 18.4, 18.5, 18.6, 18.7])
gl_inset.ylocator = mticker.FixedLocator(
    [-34.8, -34.7, -34.6, -34.5, -34.4, -34.3, -34.2, -34.1, -34.0])

# Add white solid box around Saint Helena Bay on the parent map (ax)
rect_saint_helena = Rectangle((17.9, -32.8), 0.5, 1.1, linewidth=2, edgecolor='red',
                              linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_saint_helena)
ax.annotate('St Helena Bay', xy=(17.7, -32.5), xytext=(17.7, -32.5),
            textcoords='data', fontsize=12, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

# Add red solid box around Saldanha Bay on the parent map (ax)
rect_saldanha = Rectangle((17.8, -33.4), 0.5, 0.5, linewidth=2, edgecolor='red',
                          linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_saldanha)
ax.annotate('Saldanha Bay', xy=(17.8, -33.4), xytext=(17.6, -33.3),
            textcoords='data', fontsize=12, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

# Add white dashed box around Atlantic Seaboard (Cape Town area) on the parent map (ax)
rect_atlantic = Rectangle((18.0, -34.2), 0.4, 0.5, linewidth=2, edgecolor='red',
                          linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_atlantic)
ax.annotate('Atlantic Seaboard', xy=(17.7, -34.2), xytext=(17.7, -34.1),
            textcoords='data', fontsize=12, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

# Add red dashed box around False Bay on the parent map (ax)
rect_false_bay = Rectangle((18.4, -34.4), 0.6, 0.4, linewidth=2, edgecolor='red',
                           linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_false_bay)
ax.annotate('False Bay', xy=(18.5, -34.45), xytext=(18.5, -34.45),
            textcoords='data', fontsize=12, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))


# Add red dashed box on the parent map (ax)
rect = Rectangle((18.25, -34.4), 0.3, 0.5, linewidth=2, edgecolor='black',
                 linestyle=':', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect)


# Save the plot as a .png file in the same directory as the script
plt.savefig(savepath + 'UTR_&_Surface_stations.jpg')

# Show the plot
plt.show()


# %% perfect one with other stations

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:05:34 2024
@author: nkululeko
"""


# Directory containing the NetCDF files for AJSMIT_UTR
netcdf_dir = '/media/disk02/DATA-20231010T133411Z-003/DATA/AJSMIT_UTR/Proccessed/NetCDF_Stations'

# Directory for GPITCHER_SHB data
gpitcher_shb_dir = '/media/disk02/Nkululeko_Work_Laptop/Linux_Ubuntu_Nkululeko_WSL/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/GPITCHER_SHB/'

# Corrected GPITCHER_SHB files to be plotted
gpitcher_files = ["Validated_WQM_20m_TS.nc", "Validated_WQM_70m_TS.nc"]

# Initialize a list to store station data
stations = []

# Iterate over all NetCDF files in the AJSMIT_UTR directory
for filename in os.listdir(netcdf_dir):
    if filename.endswith('.nc'):
        filepath = os.path.join(netcdf_dir, filename)
        # Open the NetCDF file
        ds = xr.open_dataset(filepath)
        # Extract the station, lat, and lon
        station = ds['station'].values[0]
        lat = ds['lat'].values[0]
        lon = ds['lon'].values[0]
        # Append the data to the list
        stations.append({'station': station, 'lat': lat, 'lon': lon})

# Create a DataFrame from the collected station data
stations_df = pd.DataFrame(stations)

# Define map extent (adjust as needed)
extent = [17, 19.5, -34.5, -31.5]  # [min_lon, max_lon, min_lat, max_lat]

# Define unique markers
markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'H',
           'd', '<', '>', 'P', 'X', 'x']  # Changed '|' to 's'

# Variable for marker size (adjust as needed)
marker_size = 11

# Plot the main map
fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(extent, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.RIVERS)
ax.add_feature(cfeature.LAKES)

# Plot each AJSMIT_UTR station in red with unique markers
for station, marker in zip(stations_df['station'].unique(), markers):
    station_data = stations_df[stations_df['station'] == station]
    # Make '*' and '^' markers red
    color = 'red' if marker == '*' or marker == '^' else 'k'
    plt.scatter(station_data['lon'], station_data['lat'], marker=marker,
                color=color, s=marker_size**2, transform=ccrs.PlateCarree(), label=station)

# Now plot the GPITCHER_SHB files in green
for file, marker in zip(gpitcher_files, markers):
    filepath = os.path.join(gpitcher_shb_dir, file)
    ds = xr.open_dataset(filepath)
    lat = ds['latitude'].values[0]
    lon = ds['longitude'].values[0]
    plt.scatter(lon, lat, marker=marker, color='green',
                s=marker_size**2, transform=ccrs.PlateCarree(), label=file)

# Set plot title and labels
plt.title('Map of Validation UTR and GPITCHER_SHB Stations')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Add gridlines and labels to main map
gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
gl.top_labels = False
gl.right_labels = False
gl.xlines = True
gl.ylines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([17, 18, 19, 20])
gl.ylocator = mticker.FixedLocator([-34.5, -34, -33.5, -33, -32.5, -32, -31.5])

# Create legend with unique markers for each station
legend_elements = [Line2D([0], [0], marker=marker, color='k' if marker != '*' and marker != '^' else 'red', label=station, markersize=marker_size, linestyle='None')
                   for station, marker in zip(stations_df['station'].unique(), markers)]
legend_elements += [Line2D([0], [0], marker=marker, color='green', label=file, markersize=marker_size, linestyle='None')
                    for file, marker in zip(gpitcher_files, markers)]

plt.legend(handles=legend_elements, loc='upper right', fontsize='small')

# Add white solid box around Saint Helena Bay on the parent map (ax)
rect_saint_helena = Rectangle((17.9, -32.8), 0.5, 1.1, linewidth=2, edgecolor='red',
                              linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_saint_helena)

# Add red solid box around Saint Helena Bay on the parent map (ax)
rect_saint_helena = Rectangle((17.8, -33.4), 0.5, 0.5, linewidth=2, edgecolor='red',
                              linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_saint_helena)

# Add white dashed box around Cape Town area on the parent map (ax)
rect_cape_town = Rectangle((18.0, -34.2), 0.4, 0.5, linewidth=2, edgecolor='red',
                           linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_cape_town)

# Add red dashed box on the parent map (ax)
rect = Rectangle((18.4, -34.4), 0.6, 0.4, linewidth=2, edgecolor='red',
                 linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect)

# Save the plot as a .png file in the same directory as the script
plt.savefig('UTR_and_GPITCHER_SHB_stations.jpg')

# Show the plot
plt.show()


# %% Merge the above 2: perfect one with other stations

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:05:34 2024
@author: nkululeko
"""


# Directory containing the NetCDF files for AJSMIT_UTR
netcdf_dir = '/media/disk02/DATA-20231010T133411Z-003/DATA/AJSMIT_UTR/Proccessed/NetCDF_Stations'

# Directory for GPITCHER_SHB data
gpitcher_shb_dir = '/media/disk02/Nkululeko_Work_Laptop/Linux_Ubuntu_Nkululeko_WSL/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/GPITCHER_SHB/'

ADCP_shb_dir = '/media/disk02/Nkululeko_Work_Laptop/Linux_Ubuntu_Nkululeko_WSL/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/adcp_mooring/'

adcp_files = ["Validated_adcp_mooring_1.nc"]
adcp_leg = ["ADCP & Wirewalker_mooring"]

# Corrected GPITCHER_SHB files to be plotted
gpitcher_files = ["Validated_WQM_20m_TS.nc", "Validated_WQM_70m_TS.nc"]
gpitcher_leg = ["WQM_20m_TS", "WQM_70m_TS"]

# Initialize a list to store station data
stations = []

# Iterate over all NetCDF files in the AJSMIT_UTR directory
for filename in os.listdir(netcdf_dir):
    if filename.endswith('.nc'):
        filepath = os.path.join(netcdf_dir, filename)
        # Open the NetCDF file
        ds = xr.open_dataset(filepath)
        # Extract the station, lat, and lon
        station = ds['station'].values[0]
        lat = ds['lat'].values[0]
        lon = ds['lon'].values[0]
        # Append the data to the list
        stations.append({'station': station, 'lat': lat, 'lon': lon})

# Create a DataFrame from the collected station data
stations_df = pd.DataFrame(stations)

# Define map extent (adjust as needed)
extent = [17, 19.5, -34.5, -31.5]  # [min_lon, max_lon, min_lat, max_lat]

# Define zoomed-in extent
# [min_lon, max_lon, min_lat, max_lat]
zoom_extent = [18.25, 18.55, -34.4, -33.9]


# Define unique markers
markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'H',
           'd', '<', '>', 'P', 'X', 'x']  # Changed '|' to 's'

# Variable for marker size (adjust as needed)
marker_size = 11

# Plot the main map
fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(extent, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.RIVERS)
ax.add_feature(cfeature.LAKES)

# Plot each AJSMIT_UTR station in red with unique markers
for station, marker in zip(stations_df['station'].unique(), markers):
    station_data = stations_df[stations_df['station'] == station]
    # Make '*' and '^' markers red
    color = 'red' if marker == '*' or marker == '^' else 'k'
    plt.scatter(station_data['lon'], station_data['lat'], marker=marker,
                color=color, s=marker_size**2, transform=ccrs.PlateCarree(), label=station)

# Now plot the GPITCHER_SHB files in green
for file, marker in zip(gpitcher_files, markers):
    filepath = os.path.join(gpitcher_shb_dir, file)
    ds = xr.open_dataset(filepath)
    lat = ds['latitude'].values[0]
    lon = ds['longitude'].values[0]
    plt.scatter(lon, lat, marker=marker, color='green',
                s=marker_size**2, transform=ccrs.PlateCarree(), label=file)

# Now plot the GPITCHER_SHB files in green
for file, marker in zip(adcp_files, markers):
    filepath = os.path.join(ADCP_shb_dir, file)
    ds = xr.open_dataset(filepath)
    lat = ds['latitude'].values[0]
    lon = ds['longitude'].values[0]
    plt.scatter(lon, lat, marker=marker, color='red', s=marker_size **
                2, transform=ccrs.PlateCarree(), label=file)


# Set plot title and labels
plt.title('Map of UTR and surface Insitu Stations used for validation', fontsize=16)
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Add gridlines and labels to main map
gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
gl.top_labels = False
gl.right_labels = False
gl.xlines = True
gl.ylines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([17, 18, 19, 20])
gl.ylocator = mticker.FixedLocator([-34.5, -34, -33.5, -33, -32.5, -32, -31.5])

# Create legend with unique markers for each station
# Create legend with unique markers for each station
legend_elements = [Line2D([0], [0], marker=marker, color='k' if marker != '*' and marker != '^' else 'red', label=station, markersize=marker_size, linestyle='None')
                   for station, marker in zip(stations_df['station'].unique(), markers)]
legend_elements += [Line2D([0], [0], marker=marker, color='green', label=file, markersize=marker_size, linestyle='None')
                    for file, marker in zip(gpitcher_leg, markers)]
legend_elements += [Line2D([0], [0], marker=marker, color='red', label=file, markersize=marker_size, linestyle='None')
                    for file, marker in zip(adcp_leg, markers)]

plt.legend(handles=legend_elements, loc='upper right', fontsize='small')

# Create inset map for zoomed-in region
ax_inset = fig.add_axes([0.55, 0.275, 0.3, 0.32],
                        projection=ccrs.PlateCarree())
ax_inset.set_extent(zoom_extent, crs=ccrs.PlateCarree())
ax_inset.add_feature(cfeature.COASTLINE)
ax_inset.add_feature(cfeature.BORDERS, linestyle=':')
ax_inset.add_feature(cfeature.LAND)
ax_inset.add_feature(cfeature.OCEAN)
ax_inset.add_feature(cfeature.RIVERS)
ax_inset.add_feature(cfeature.LAKES)

# Plot stations in the zoomed-in region
for station, marker in zip(stations_df['station'].unique(), markers):
    station_data = stations_df[stations_df['station'] == station]
    # Make '*' and '^' markers red
    color = 'red' if marker == '*' or marker == '^' else 'k'
    ax_inset.scatter(station_data['lon'], station_data['lat'], marker=marker,
                     color=color, s=marker_size**2, transform=ccrs.PlateCarree())

# Add gridlines and labels to inset map
gl_inset = ax_inset.gridlines(draw_labels=True, linestyle='--', color='gray')
gl_inset.top_labels = False
gl_inset.right_labels = False
gl_inset.xlines = True
gl_inset.ylines = True
gl_inset.xformatter = LONGITUDE_FORMATTER
gl_inset.yformatter = LATITUDE_FORMATTER
gl_inset.xlocator = mticker.FixedLocator([18.3, 18.4, 18.5, 18.6, 18.7])
gl_inset.ylocator = mticker.FixedLocator(
    [-34.8, -34.7, -34.6, -34.5, -34.4, -34.3, -34.2, -34.1, -34.0])


# Add white solid box around Saint Helena Bay on the parent map (ax)
rect_saint_helena = Rectangle((17.9, -32.8), 0.5, 1.1, linewidth=2, edgecolor='red',
                              linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_saint_helena)
ax.annotate('St Helena Bay', xy=(17.7, -32.5), xytext=(17.7, -32.5),
            textcoords='data', fontsize=12, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

# Add red solid box around Saldanha Bay on the parent map (ax)
rect_saldanha = Rectangle((17.8, -33.4), 0.5, 0.5, linewidth=2, edgecolor='red',
                          linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_saldanha)
ax.annotate('Saldanha Bay', xy=(17.8, -33.4), xytext=(17.6, -33.3),
            textcoords='data', fontsize=12, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

# Add white dashed box around Atlantic Seaboard (Cape Town area) on the parent map (ax)
rect_atlantic = Rectangle((18.0, -34.2), 0.4, 0.5, linewidth=2, edgecolor='red',
                          linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_atlantic)
ax.annotate('Atlantic Seaboard', xy=(17.7, -34.2), xytext=(17.7, -34.1),
            textcoords='data', fontsize=12, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

# Add red dashed box around False Bay on the parent map (ax)
rect_false_bay = Rectangle((18.4, -34.4), 0.6, 0.4, linewidth=2, edgecolor='red',
                           linestyle='-', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect_false_bay)
ax.annotate('False Bay', xy=(18.5, -34.45), xytext=(18.5, -34.45),
            textcoords='data', fontsize=12, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))


# Add red dashed box on the parent map (ax)
rect = Rectangle((18.25, -34.4), 0.3, 0.5, linewidth=2, edgecolor='black',
                 linestyle=':', facecolor='none', transform=ccrs.PlateCarree())
ax.add_patch(rect)


# Save the plot as a .png file in the same directory as the script
plt.savefig(savepath + 'UTR_&_Surface_stations.jpg')

# Show the plot
plt.show()


# %%

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:05:34 2024
@author: nkululeko
"""

import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Define the sub-regions
sub_regions = {
    # [min_lon, max_lon, min_lat, max_lat]
    "St Helena Bay": [17.9, 18.7, -32.8, -31.7],
    "Saldanha Bay": [17.8, 18.3, -33.4, -32.9],
    "Atlantic Seaboard": [18.0, 18.4, -34.2, -33.7],
    "False Bay": [18.4, 19.0, -34.4, -34.0]
}

# Directory containing the NetCDF files for AJSMIT_UTR
netcdf_dir = '/media/disk02/DATA-20231010T133411Z-003/DATA/AJSMIT_UTR/Proccessed/NetCDF_Stations'
savepath = '/home/nc.memela/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/VALIDATION_Report/Combined_Validation/'

# Directory for GPITCHER_SHB data
gpitcher_shb_dir = '/media/disk02/Nkululeko_Work_Laptop/Linux_Ubuntu_Nkululeko_WSL/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/GPITCHER_SHB/'
# Directory for ADCP data
adcp_shb_dir = '/media/disk02/Nkululeko_Work_Laptop/Linux_Ubuntu_Nkululeko_WSL/somisana-croco/configs/swcape_02/croco_v1.3.1/C01_I01_GLORYS_ERA5/Validation/adcp_mooring/'

# List of ADCP files
adcp_files = ["Validated_adcp_mooring_1.nc"]
# Corrected GPITCHER_SHB files to be plotted
gpitcher_files = ["Validated_WQM_20m_TS.nc", "Validated_WQM_70m_TS.nc"]

# Initialize a list to store station data
stations = []

# Define markers for different datasets
ajmsit_markers = ['o', 's', '^', 'D', 'v', 'p', 'h']  # Red markers for AJSMIT_UTR
gpitcher_markers = ['o', 's', '^', 'D']  # Green markers for GPITCHER_SHB
adcp_marker = 'o'  # Pink markers for ADCP

# Read AJSMIT_UTR data (red markers)
for filename in os.listdir(netcdf_dir):
    if filename.endswith('.nc'):
        filepath = os.path.join(netcdf_dir, filename)
        # Open the NetCDF file
        ds = xr.open_dataset(filepath)
        # Extract the station, lat, and lon
        station = ds['station'].values[0]
        lat = ds['lat'].values[0]
        lon = ds['lon'].values[0]
        # Append the data to the list with red marker designation and unique markers
        stations.append({'station': station, 'lat': lat, 'lon': lon, 'color': 'red', 'marker': ajmsit_markers[len(stations) % len(ajmsit_markers)]})

# Read GPITCHER_SHB data (green markers)
for filename in gpitcher_files:
    filepath = os.path.join(gpitcher_shb_dir, filename)
    ds = xr.open_dataset(filepath)
    # Extract latitude and longitude
    lat = ds['latitude'].values
    lon = ds['longitude'].values
    # Append the data to the list for all stations in the GPITCHER_SHB dataset
    for station_lat, station_lon in zip(lat, lon):
        stations.append({'station': filename[10:].replace('_TS.nc', ''),
                         'lat': station_lat, 'lon': station_lon, 'color': 'green', 'marker': gpitcher_markers[len(stations) % len(gpitcher_markers)]})

# Read ADCP data (pink markers)
for filename in adcp_files:
    filepath = os.path.join(adcp_shb_dir, filename)
    ds = xr.open_dataset(filepath)
    # Extract latitude and longitude
    lat = ds['latitude'].values
    lon = ds['longitude'].values
    # Append the data to the list for all stations in the ADCP dataset
    for station_lat, station_lon in zip(lat, lon):
        stations.append({'station': filename[10:].replace('_1.nc', ''),
                         'lat': station_lat, 'lon': station_lon, 'color': 'pink', 'marker': adcp_marker})

# Convert the station data to a list of dictionaries
stations_df = stations

# Function to create a map for each region and save it to a file
def plot_region(region_name, extent, annotate_text):
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    # Add features to the map
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)

    # Plot the stations within the sub-region
    for station in stations_df:
        lon = station['lon']
        lat = station['lat']
        color = station['color']  # Get the color based on the source
        marker = station['marker']  # Get the marker for the station
        if extent[0] <= lon <= extent[1] and extent[2] <= lat <= extent[3]:
            # Plot stations inside this sub-region with specified markers and black outline
            plt.scatter(lon, lat, marker=marker, s=200, facecolor=color,
                        edgecolor='black', linewidth=1.5, label=station['station'], transform=ccrs.PlateCarree())

    # Add title and annotation
    plt.title(f"Map of {region_name}", fontsize=16)
    ax.text(extent[0] + 0.05, extent[2] + 0.15, annotate_text, fontsize=12,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

    # Add legend
    plt.legend(loc='upper right', fontsize='small')

    # Save the figure to the specified save path
    plt.savefig(os.path.join(savepath, f"{region_name.replace(' ', '_')}_map.png"), bbox_inches='tight')
    plt.close(fig)  # Close the figure to avoid displaying it


# Plot maps for each sub-region and save them
plot_region("St Helena Bay", sub_regions["St Helena Bay"], "St Helena Bay")
plot_region("Saldanha Bay", sub_regions["Saldanha Bay"], "Saldanha Bay")
plot_region("Atlantic Seaboard",
            sub_regions["Atlantic Seaboard"], "Atlantic Seaboard")
plot_region("False Bay", sub_regions["False Bay"], "False Bay")
