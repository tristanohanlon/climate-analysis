# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon

This will create a dataset of time averaged global total cloud cover with latitude.
Data is stored in the 2D array gfdl_tcc_lat.

[:,0] = latitude
[:,1] = cloud fraction
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt

dataset = Dataset('2010/clt_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')

lon = dataset.variables['lon'][:] #Extract longitude data
lat = dataset.variables['lat'][:] #Extract latitude data
gfdl_tcc = dataset.variables['clt'][:] #Extract total cloud cover data, keyed to time,lon and lat
        
#Select the months from 2010
#gfdl_tc = np.take(gfdl_tcc, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)


tc = np.mean(gfdl_tcc, axis=0) # average total cloud cover data over time
tcc = np.mean(tc, axis=1) # average total cloud cover data over longitude

gfdl_tcc_lat = np.vstack((lat, tcc)).T # Join the two lists as if they were two columns side by side, into a list of two elements each

plt.figure()
plt.plot(gfdl_tcc_lat[:,0],gfdl_tcc_lat[:,1])
plt.show()