# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: tristan
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt

dataset = Dataset('ERA5.2016.01.nc', 'r')

lon = dataset.variables['longitude'][:] #Extract longitude data
lat = dataset.variables['latitude'][:] #Extract latitude data
tcc = dataset.variables['tcc'][:]
tciw = dataset.variables['tciw'][:]
tclw = dataset.variables['tclw'][:]


tcc_av = np.mean(tcc, axis=0)
tciw_av = np.mean(tciw, axis=0)
tclw_av = np.mean(tclw, axis=0)

lat_bins=np.arange(-70,-50,0.5)
lat_hist=np.digitize(lat,bins=lat_bins)
tcclist=[]
tciwlist=[]
tclwlist=[]
for i in range(len(lat_bins)):
    location=lat_hist==i
    tcc_lat=tcc_av[location]
    nw_tcc_data=np.nanmean(tcc_lat[(tcc_lat<1.0)&(tcc_lat>-0.01)])
    tcclist.append(nw_tcc_data)
    
    tciw_lat=tciw_av[location]
    nw_tciw_data=np.nanmean(tciw_lat[(tciw_lat<1.0)&(tciw_lat>-0.01)])
    tciwlist.append(nw_tciw_data)
    
    tclw_lat=tclw_av[location]
    nw_tclw_data=np.nanmean(tclw_lat[(tclw_lat<1.0)&(tclw_lat>-0.01)])
    tclwlist.append(nw_tclw_data)

plt.show()
plt.figure()

fig, ax = plt.subplots()
ax.plot(lat_bins,tcclist, '-b', label='Total Cloud Cover')
ax.plot(lat_bins,tciwlist, '--b', label='Ice Cloud Cover')
ax.plot(lat_bins,tclwlist, '-g', label='Liquid Cloud Cover')
#ax.axis('equal')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=3, fancybox=True, shadow=True);


plt.xlabel('Latitude')
plt.ylabel('Cloud Fraction')
plt.title('Cloud Fraction vs Latitude - ECMWF Reanalysis 01.2016')
plt.grid(True)