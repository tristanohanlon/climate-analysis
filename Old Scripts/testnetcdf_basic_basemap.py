# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 12:55:59 2019

@author: toha006
"""
import datetime as dt  # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

dataset = Dataset('Data/19820101.atmos_month-001.nc', 'r') #define dataset file in the same directory as read only


lons = dataset.variables['lon'][:] #Extract longitude data
lats = dataset.variables['lat'][:] #Extract latitude data
t_surf = dataset.variables['t_surf'][:] #Extract surface tempertaure data, keyed to time, lat and long

t_surf_units = dataset.variables['t_surf'].units #Extract surface tempertaure units

t_surfav = np.mean(t_surf, axis=0) #Get mean surface tempertaure over time (first key array column)

#lon_0 = lons.mean() #centres on the mean longitude
#lat_0 = lats.mean() #centres on the mean latitude

m = Basemap(width=5000000,height=3500000, 
            resolution='l',projection='stere',
            lat_ts=0,lat_0=-40,lon_0=180)

lon, lat = np.meshgrid(lons, lats) # create a 2D array
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(t_surfav))

# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawcountries()

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(t_surf_units)

# Add Title
plt.title('Average Surface Temperature')

plt.show()