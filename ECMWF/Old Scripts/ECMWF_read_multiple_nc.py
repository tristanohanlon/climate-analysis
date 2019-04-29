# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: toha006
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from netCDF4 import MFDataset
import xarray as xr
import matplotlib.pyplot as plt

dataset = xr.open_mfdataset('2010/*.nc', chunks={'time': 200})
#dataset = MFDataset('2010/*.nc', 'r')

lon = dataset.variables['longitude'][:] #Extract longitude data
lat = dataset.variables['latitude'][:] #Extract latitude data
tcc = dataset.variables['tcc'][:] #Extract surface tempertaure data, keyed to time, lat and long
      
#data = (tcc - add_offset) * scale_factor

tcloud_av = np.mean(tcc, axis=0) #Get mean surface tempertaure over time (first key array column)

plt.figure()
plt.contourf(lon, lat, tcloud_av)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Total Cloud Fraction - ECMWF Reanalysis 2010')


lat_bins=np.arange(-90,90,0.5)
lat_hist=np.digitize(lat,bins=lat_bins)
list1=[]
for i in range(len(lat_bins)):
    location=lat_hist==i
    cf_lat=tcloud_av[location]
    nw_data=np.nanmean(cf_lat[(cf_lat<1.0)&(cf_lat>-0.01)])
    list1.append(nw_data)
plt.figure()
plt.plot(lat_bins,list1)
plt.xlabel('Latitude')
plt.ylabel('Total Cloud Fraction')
plt.title('Total Cloud Fraction vs Latitude - ECMWF Reanalysis 2010')
plt.grid(True)