# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:13:24 2019

@author: Tristan
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from pyhdf import SD
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from netCDF4 import MFDataset
import xarray as xr
import matplotlib.pyplot as plt

# run get data xarray first

elon = dataset.variables['longitude'][:] #Extract longitude data
elat = dataset.variables['latitude'][:] #Extract latitude data
etcc = dataset.variables['tcc'][:] #Extract surface tempertaure data, keyed to time, lat and long

tcloud_av = np.mean(etcc, axis=0) #Get mean total cloud cover over time (first key array column)
tcloud_av1 = np.array(tcloud_av) # convert xarray back to numpy array

lat_bins=np.arange(-90,90,0.5)
lat_hist=np.digitize(elat,bins=lat_bins)
list1=[]
for i in range(len(lat_bins)):
    location=lat_hist==i
    cf_lat=tcloud_av1[location]
    nw_data=np.nanmean(cf_lat[(cf_lat<1.0)&(cf_lat>-0.01)])
    list1.append(nw_data)
    
f=SD.SD('CERES/2010/CERES_2010_tcc.hdf')
#view datasets
f.datasets()

tc=f.select('tcc') 
lat=f.select('latitude (dimension)')
lon=f.select('longitude (dimension)')
tc=tc.get()
lon=lon.get()
lat=lat.get()

a = np.reshape(tc,180)

tcc=a/100
    
plt.figure()

fig, ax = plt.subplots()
ax.plot(lat,tcc, '-b', label='CERES Total Cloud Cover')
ax.plot(lat_bins,list1, '-r', label='ECMWF Total Cloud Cover')
#ax.axis('equal')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

#plt.plot(lat[140:161],tcc[140:161])
plt.xlabel('Latitude')
plt.ylabel('Cloud Fraction')
plt.title('Total Cloud Fraction vs Latitude - CERES and ECMWF Reanalysis 2010')
plt.grid(True)
plt.show()





















