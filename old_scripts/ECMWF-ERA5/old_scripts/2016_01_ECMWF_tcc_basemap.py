# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: tristan
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm


dataset = Dataset('ERA5.2016.01.nc', 'r')

lons = dataset.variables['longitude'][:] #Extract longitude data
lats = dataset.variables['latitude'][:] #Extract latitude data
tcc = dataset.variables['tcc'][:] #Extract surface tempertaure data, keyed to time, lat and long

tcloud_av = np.mean(tcc, axis=0) #Get mean surface tempertaure over time (first key array column)

#lon_0 = lons.mean() #centres on the mean longitude
#lat_0 = lats.mean() #centres on the mean latitude
plt.figure(figsize = (12,8)) #set size of figure
m = Basemap(llcrnrlon=0.,llcrnrlat=-70.,urcrnrlon=360.,urcrnrlat=-50.,
            rsphere=(6378137, 6356752),
            resolution='l',projection='merc',
            lat_0=-50.,lon_0=180.,lat_ts=10)

lon, lat = np.meshgrid(lons, lats) # create a 2D array
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(tcloud_av))

# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 20.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawcountries()

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="20%")
#cbar.set_label(tcloud_units)

# Add Title
plt.title('Total Cloud Cover over Southern Ocean')

#m.plot()
plt.show()