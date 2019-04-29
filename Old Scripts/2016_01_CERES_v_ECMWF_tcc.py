# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:13:24 2019

@author: Tristan
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from pyhdf import SD
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt

dataset = Dataset('ECMWF/ERA5.2016.01.nc', 'r')

elon = dataset.variables['longitude'][:] #Extract longitude data
elat = dataset.variables['latitude'][:] #Extract latitude data
etcc = dataset.variables['tcc'][:] #Extract surface tempertaure data, keyed to time, lat and long

#add_offset = dataset.variables['tcc'].add_offset
#scale_factor = dataset.variables['tcc'].scale_factor
        
#data = (tcc - add_offset) * scale_factor

tcloud_av = np.mean(etcc, axis=0) #Get mean surface tempertaure over time (first key array column)

lat_bins=np.arange(-70,-50,0.5)
lat_hist=np.digitize(elat,bins=lat_bins)
list1=[]
for i in range(len(lat_bins)):
    location=lat_hist==i
    cf_lat=tcloud_av[location]
    nw_data=np.nanmean(cf_lat[(cf_lat<1.0)&(cf_lat>-0.01)])
    list1.append(nw_data)
    
f=SD.SD('CERES/CER_SSF1deg-Month_Terra-MODIS_Edition4A_401405.201601.hdf')
#view datasets
f.datasets()

tc=f.select('cld_amount_zon') 
tlc=f.select('cld_amount_liq_zon') 
tic=f.select('cld_amount_ice_zon') 
lat=f.select('latitude')
lon=f.select('longitude')
tc=tc.get()
tlc=tlc.get()
tic=tic.get()
lon=lon.get()
lat=lat.get()
tpcc=tc[-1,:] #selects the last index of the array (total cloud)
tcc=tpcc/100
tlpcc=tlc[-1,:]
tlcc=tlpcc/100
tipcc=tic[-1,:]
ticc=tipcc/100
    
plt.figure()

fig, ax = plt.subplots()
ax.plot(lat[140:161],tcc[140:161], '-b', label='CERES Total Cloud Cover')
ax.plot(lat[140:161],tlcc[140:161], '--b', label='CERES Liquid Cloud Cover')
ax.plot(lat[140:161],ticc[140:161], '--o', label='CERES Ice Cloud Cover')
ax.plot(lat_bins,list1, '-r', label='ECMWF Total Cloud Cover')
#ax.axis('equal')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

#plt.plot(lat[140:161],tcc[140:161])
plt.xlabel('Latitude')
plt.ylabel('Cloud Fraction')
plt.title('Total Cloud Fraction vs Latitude - CERES and ECMWF Reanalysis 01.2016')
plt.grid(True)
plt.show()





















