
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
tcloud = dataset.variables['tot_cld_amt'][:] #Extract surface tempertaure data, keyed to time, lat and long

tcloud_units = dataset.variables['tot_cld_amt'].units #Extract surface tempertaure units

tcloud_av = np.mean(tcloud, axis=0) #Get mean surface tempertaure over time (first key array column)

#lon_0 = lons.mean() #centres on the mean longitude
#lat_0 = lats.mean() #centres on the mean latitude
plt.figure(figsize = (12,8)) #set size of figure
m = Basemap(llcrnrlon=0.,llcrnrlat=-70.,urcrnrlon=360.,urcrnrlat=-50.,
            rsphere=(6378137, 6356752),
            resolution='l',projection='merc',
            lat_0=-60.,lon_0=180.,lat_ts=5)

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
cbar = m.colorbar(cs, location='bottom', pad="30%")
cbar.set_label('Cloud Fraction')

# Add Title
plt.title('Total Cloud Cover over Southern Ocean - ECMWF (ERA5 Reanalysis)')

#m.plot()
plt.show()



#fig,ax=plt.subplots(figsize=(10,10),dpi=300)
#ax.imshow()
#ax.set_xticks(np.arange(0,180,10))
#ax.set_xticklabels(np.arange(-90,90,10))

